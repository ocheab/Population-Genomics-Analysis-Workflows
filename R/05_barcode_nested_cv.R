source('R/01_helpers.R')
cat('WP5-WP6: Leakage-controlled candidate barcode evaluation (exploratory)\n')
Z<-readRDS(required(file.path(OUT,'rds','QC_cohorts.rds')))$primary2020
X<-t(Z$X$G); y<-Z$meta$state; sample_ids<-Z$meta$sample_id
assert(all(y %in% STATES),'Unexpected class labels')
assert(min(table(factor(y,STATES)))>=SETTINGS$outer_folds,'Too few samples for outer CV')
assert(length(unique(sample_ids))==length(sample_ids),'Duplicate sample IDs')
# Group-aware folds by recorded case identifiers, when duplicate case groups exist.
case<-Z$meta$case_id
case[is.na(case)|case=='']<-sample_ids[is.na(case)|case=='']
assert(!any(tapply(y,case,function(v)length(unique(v))>1L)),
  'Cross-state case groups detected: review before modelling')
make_folds<-function(y,k,groups=NULL,seed=1L) {
  set.seed(seed)
  if(is.null(groups)) groups<-as.character(seq_along(y))
  ug<-unique(groups)
  gs<-vapply(ug,function(g)unique(y[groups==g])[1],character(1))
  fold_by_g<-setNames(integer(length(ug)),ug)
  for (s in STATES) {
    g<-sample(ug[gs==s]);assert(length(g)>=k,paste('Too few independent',s,'groups for folds'))
    # Distribute group sizes greedily to reduce fold imbalance.
    loads<-integer(k)
    for (id in g) { f<-which.min(loads);fold_by_g[id]<-f;loads[f]<-loads[f]+sum(groups==id) }
  }
  as.integer(fold_by_g[groups])
}
# Training-fold-only site screening: biallelic, quality and MAC, THEN label ranking and LD pruning.
rank_training_sites<-function(X,y,max_k=max(SETTINGS$panel_sizes)) {
  n<-nrow(X); P<-ncol(X)
  called<-!is.na(X)
  invalid<-colSums(called & X!=0L & X!=1L)>0
  ni<-colSums(called)
  ai<-colSums(X==1L,na.rm=TRUE)
  state_n<-sapply(STATES,function(s) colSums(called[y==s,,drop=FALSE]))
  state_alt<-sapply(STATES,function(s) colSums(X[y==s,,drop=FALSE]==1L,na.rm=TRUE))
  denom<-vapply(STATES,function(s)sum(y==s),integer(1))
  state_rate<-sweep(state_n,2,pmax(1,denom),'/')
  eligible <- !invalid & ni/n >= SETTINGS$call_rate &
    rowSums(state_rate>=SETTINGS$state_call_rate)==length(STATES) &
    pmin(ai,ni-ai)>=SETTINGS$minor_allele_count & ai>0 & ai<ni
  ix<-which(eligible);if(!length(ix))return(integer())
  p<-ai[ix]/ni[ix]
  exp_alt<-sweep(state_n[ix,,drop=FALSE],1,p,'*')
  exp_ref<-state_n[ix,,drop=FALSE]-exp_alt
  obs_alt<-state_alt[ix,,drop=FALSE]
  chi<-rowSums((obs_alt-exp_alt)^2/pmax(exp_alt,1e-9)+
               ((state_n[ix,,drop=FALSE]-obs_alt)-exp_ref)^2/pmax(exp_ref,1e-9))
  ord<-ix[order(-chi,ix)]
  keep<-integer()
  for (j in ord) {
    redundant<-FALSE
    for (k in keep) {
      ok<-!is.na(X[,j]) & !is.na(X[,k])
      if(sum(ok)>=20 && length(unique(X[ok,j]))==2 && length(unique(X[ok,k]))==2) {
        r2<-suppressWarnings(cor(X[ok,j],X[ok,k])^2)
        if(is.finite(r2)&&r2>=SETTINGS$ld_prune_r2){redundant<-TRUE;break}
      }
    }
    if(!redundant) keep<-c(keep,j)
    if(length(keep)>=max_k)break
  }
  keep
}
# Equal state priors; score averages log likelihood across observed panel calls.
# Scores are NOT calibrated posterior probabilities, especially in one linked gene.
predict_panel<-function(Xtrain,ytrain,Xtest,feat,min_fraction=.75) {
  if(!length(feat))return(rep(NA_character_,nrow(Xtest)))
  tr<-Xtrain[,feat,drop=FALSE];te<-Xtest[,feat,drop=FALSE]
  pf<-sapply(STATES,function(s){
    z<-tr[ytrain==s,,drop=FALSE]
    (colSums(z==1L,na.rm=TRUE)+1)/(colSums(!is.na(z))+2)
  })
  if(is.null(dim(pf))) pf<-matrix(pf,nrow=length(feat),dimnames=list(NULL,STATES))
  score<-matrix(NA_real_,nrow(te),length(STATES),dimnames=list(NULL,STATES))
  for (i in seq_len(nrow(te))) {
    obs<-which(!is.na(te[i,]))
    if(length(obs)<max(1,ceiling(min_fraction*length(feat))))next
    for(s in seq_along(STATES)) {
      f<-pmin(pmax(pf[obs,s],1e-8),1-1e-8)
      score[i,s]<-mean(ifelse(te[i,obs]==1L,log(f),log(1-f)))
    }
  }
  apply(score,1,function(v)if(all(is.na(v)))NA_character_ else STATES[which.max(v)])
}
metrics<-function(truth,pred) {
  recall<-precision<-f1<-setNames(rep(0,length(STATES)),STATES)
  for(s in STATES) {
    tp<-sum(truth==s & !is.na(pred) & pred==s)
    fn<-sum(truth==s)-tp
    fp<-sum(truth!=s & !is.na(pred) & pred==s)
    recall[s]<-if(tp+fn>0)tp/(tp+fn) else NA_real_
    precision[s]<-if(tp+fp>0)tp/(tp+fp) else 0
    f1[s]<-if(precision[s]+recall[s]>0)2*precision[s]*recall[s]/(precision[s]+recall[s]) else 0
  }
  c(balanced_accuracy=mean(recall,na.rm=TRUE),macro_F1=mean(f1,na.rm=TRUE),
    overall_accuracy=mean(!is.na(pred)&truth==pred),coverage=mean(!is.na(pred)))
}
qc_predict<-function(trn,ytr,tst) {
  # nuisance-only baseline from Fws and genotype missingness, with train-estimated scaling.
  mu<-colMeans(trn);sig<-apply(trn,2,sd);sig[sig<1e-8]<-1
  trn<-sweep(sweep(trn,2,mu,'-'),2,sig,'/')
  tst<-sweep(sweep(tst,2,mu,'-'),2,sig,'/')
  centers<-sapply(STATES,function(s)colMeans(trn[ytr==s,,drop=FALSE]))
  if(is.null(dim(centers)))centers<-matrix(centers,nrow=2)
  apply(tst,1,function(z)STATES[which.min(colSums((centers-z)^2))])
}
cv_once<-function(X,y,nuisance,groups,repeat_id=1L,seed=1L,return_detail=TRUE) {
  folds<-make_folds(y,SETTINGS$outer_folds,groups,seed)
  pred<-pred_qc<-rep(NA_character_,length(y)); chosen<-integer(length(y))
  selected<-list()
  for(f in seq_len(SETTINGS$outer_folds)) {
    test<-which(folds==f);train<-which(folds!=f)
    xt<-X[train,,drop=FALSE];yt<-y[train]; nt<-nuisance[train,,drop=FALSE]
    inner<-make_folds(yt,SETTINGS$inner_folds,groups[train],seed+1000L+f)
    panel_score<-rep(NA_real_,length(SETTINGS$panel_sizes))
    for(kk in seq_along(SETTINGS$panel_sizes)) {
      ip<-rep(NA_character_,length(train));k<-SETTINGS$panel_sizes[kk]
      for(g in seq_len(SETTINGS$inner_folds)) {
        v<-which(inner==g);t<-which(inner!=g)
        ranked<-rank_training_sites(xt[t,,drop=FALSE],yt[t],max_k=k)
        ip[v]<-predict_panel(xt[t,,drop=FALSE],yt[t],xt[v,,drop=FALSE],head(ranked,k))
      }
      panel_score[kk]<-metrics(yt,ip)['balanced_accuracy']
    }
    k<-SETTINGS$panel_sizes[which.max(panel_score)] # smaller k breaks exact ties
    ranked<-rank_training_sites(xt,yt,max_k=k)
    sel<-head(ranked,k)
    pred[test]<-predict_panel(xt,yt,X[test,,drop=FALSE],sel)
    pred_qc[test]<-qc_predict(nt,yt,nuisance[test,,drop=FALSE])
    chosen[test]<-k
    selected[[f]]<-if(length(sel)) data.frame(repeat_id=rep(repeat_id,length(sel)),
      fold=rep(f,length(sel)),site_index=sel,pos=Z$site$pos[sel]) else
      data.frame(repeat_id=integer(),fold=integer(),site_index=integer(),pos=integer())
  }
  list(metrics=metrics(y,pred),nuisance_metrics=metrics(y,pred_qc),
       prediction=data.frame(repeat_id=repeat_id,sample_id=sample_ids,truth=y,
         predicted=pred,predicted_QC_only=pred_qc,outer_fold=folds,chosen_size=chosen),
       selected=do.call(rbind,selected))
}
Nuisance<-as.matrix(Z$meta[,c('Fws','pass_missingness_haploid')]);storage.mode(Nuisance)<-'double'
results<-list();scores<-list();markers<-list()
for(r in seq_len(SETTINGS$outer_repeats)) {
  rr<-cv_once(X,y,Nuisance,case,repeat_id=r,seed=SETTINGS$seed+100L*r)
  results[[r]]<-rr$prediction
  scores[[r]]<-data.frame(repeat_id=r,t(rr$metrics),t(rr$nuisance_metrics))
  markers[[r]]<-rr$selected
  cat(sprintf('Nested CV repeat %d/%d: balance accuracy %.3f, macro-F1 %.3f\n',
      r, SETTINGS$outer_repeats,rr$metrics['balanced_accuracy'],rr$metrics['macro_F1']))
}
P<-do.call(rbind,results)
S<-do.call(rbind,scores)
names(S)<-c('repeat_id','balanced_accuracy','macro_F1','overall_accuracy','coverage',
  'QC_only_balanced_accuracy','QC_only_macro_F1','QC_only_overall_accuracy','QC_only_coverage')
M<-do.call(rbind,markers)
write_tsv(P,'Table_S14_nested_CV_predictions.tsv','model')
write_tsv(S,'Table_5_nested_CV_metrics_by_repeat.tsv','model')
if(!is.null(M)&&nrow(M)>0) {
  stability<-sort(table(M$pos),decreasing=TRUE)
  stab<-data.frame(pos=as.integer(names(stability)),outer_fold_selections=as.integer(stability),
    n_outer_models=SETTINGS$outer_repeats*SETTINGS$outer_folds)
  stab$selection_fraction<-stab$outer_fold_selections/stab$n_outer_models
  write_tsv(stab,'Table_S15_candidate_marker_selection_stability.tsv','model')
}
# ONE full nested-CV repetition for observed vs null to keep permutation specifications identical.
obs_permutation<-S$balanced_accuracy[1]
null<-numeric(SETTINGS$barcode_perm_reps)
for(b in seq_along(null)) {
  set.seed(SETTINGS$seed+100000L+b)
  yp<-sample(y) # state counts retained; full rank/filter/tune procedure re-run
  gperm<-paste('perm',seq_along(yp),sep='_') # independence under null
  tmp<-cv_once(X,yp,Nuisance,gperm,repeat_id=b,seed=SETTINGS$seed+100L,
    return_detail=FALSE)
  null[b]<-tmp$metrics['balanced_accuracy']
  if(b%%10L==0L)cat('Whole-pipeline null permutations:',b,'/',length(null),'\n')
}
p_perm<-(1+sum(null>=obs_permutation))/(length(null)+1)
write_tsv(data.frame(null_balanced_accuracy=null),'Table_S16_whole_pipeline_permutation_null.tsv','model')
report<-data.frame(model='PI4K_SNP_equal_prior_centroid_probability',
  mean_repeated_balanced_accuracy=mean(S$balanced_accuracy),
  mean_repeated_macro_F1=mean(S$macro_F1),
  mean_QC_only_balanced_accuracy=mean(S$QC_only_balanced_accuracy),
  baseline_majority_accuracy=max(table(y))/length(y),
  baseline_majority_balanced_accuracy=1/length(STATES),
  observed_one_repeat_balanced_accuracy=obs_permutation,
  n_whole_pipeline_null_permutations=length(null),
  whole_pipeline_permutation_p=p_perm,
  independence_validation='NOT_DONE', external_validation='NOT_DONE')
write_tsv(report,'Table_5_barcode_internal_validation_summary.tsv','model')
# Confusion matrix only from first out-of-fold repetition; repeats are correlated.
first<-P[P$repeat_id==1L,,drop=FALSE]
conf<-as.data.frame.matrix(table(Observed=factor(first$truth,STATES),
    Predicted=factor(first$predicted,levels=STATES)))
conf<-cbind(observed_state=rownames(conf),conf)
write_tsv(conf,'Table_S17_first_repeat_confusion_matrix.tsv','model')
per_state<-do.call(rbind,lapply(STATES,function(s){
 n<-sum(first$truth==s);k<-sum(first$truth==s&!is.na(first$predicted)&first$predicted==s)
 ci<-wilson_ci(k,n)
 data.frame(state=s,n=n,correct=k,recall=k/n,recall_CI_low=ci[1],recall_CI_high=ci[2])
}))
write_tsv(per_state,'Table_S18_first_repeat_per_state_recall.tsv','model')
# Full-data marker ranking is for ASSAY FOLLOW-UP ONLY. Never use to compute CV accuracy.
ix<-rank_training_sites(X,y,max_k=max(SETTINGS$panel_sizes))
if(length(ix)) {
 cand<-data.frame(pos=Z$site$pos[ix],ref=Z$site$ref[ix],alt=Z$site$alt[ix],
   rank=seq_along(ix),note='EXPLORATORY_FULL_DATA_RANKING_NOT_INDEPENDENTLY_VALIDATED')
 write_tsv(cand,'Candidate_assay_SNPs_NOT_VALIDATED.tsv','assay')
}
if(requireNamespace('ggplot2',quietly=TRUE)) {
 library(ggplot2)
 cmp<-rbind(data.frame(repeat_id=S$repeat_id,model='PI4K SNP',BA=S$balanced_accuracy),
   data.frame(repeat_id=S$repeat_id,model='Fws + missingness',BA=S$QC_only_balanced_accuracy))
 p<-ggplot(cmp,aes(model,BA,color=model))+geom_boxplot(width=.45,outlier.shape=NA)+
    geom_jitter(width=.09,alpha=.8,size=2)+geom_hline(yintercept=.20,linetype='dashed')+
    scale_y_continuous(limits=c(0,1))+guides(color='none')+
    labs(title='Exploratory internal prediction of collection state',
      subtitle='Repeated nested CV; 20% balanced-accuracy chance reference; NOT external validation',
      x=NULL,y='Balanced accuracy')+clean_theme()
 export_plot(p,'Figure_8_nested_CV_performance',7,5)
 cf<-as.data.frame(table(factor(first$truth,STATES),factor(first$predicted,STATES)))
 names(cf)<-c('Observed','Predicted','n')
 p<-ggplot(cf,aes(Predicted,Observed,fill=n))+geom_tile(color='white')+
    geom_text(aes(label=n),size=3.5)+scale_fill_gradient(low='#FFFFFF',high='#376285')+
    labs(title='PI4K five-state confusion matrix (first held-out repeat)',
      subtitle='Counts; not pooled across correlated CV repeats',fill='Samples')+
    theme_minimal(base_size=11)+theme(panel.grid=element_blank())
 export_plot(p,'Figure_9_outer_CV_confusion_matrix',7,5.5)
}
cat('WP5-WP6 COMPLETE — internal-only barcode evaluation; not prospective validation.\n')
