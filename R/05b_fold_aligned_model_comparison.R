#!/usr/bin/env Rscript
# PI4K Step 29 — SOURCE-AUDITED, FOLD-ALIGNED INTERNAL MODEL COMPARISON
# Run from /mnt/c/PI4K_Nigeria_Study/PI4K_Southwest_R_Analysis_Suite
# Requires previously generated results/Southwest5/R_Analysis_v1/rds/QC_cohorts.rds
# Does NOT modify R/05_barcode_nested_cv.R, the underlying genotype data, or prior results.
# The fixed-four panel was selected after viewing this same cohort: DESCRIPTIVE ONLY.
# All SNP discovery for the nested-selected model is repeated within training folds.
source('R/01_helpers.R')
cat('STEP 29 CORRECTED: repeated fold-aligned 5-state PI4K modelling (offline)\n')
Z <- readRDS(required(file.path(OUT, 'rds', 'QC_cohorts.rds')))$primary2020
G <- t(Z$X$G) # independently sampled haploid-consensus isolates x PASS SNPs
Y <- as.character(Z$meta$state)
IDS <- as.character(Z$meta$sample_id)
assert(nrow(G) == length(Y) && ncol(G) == EXPECTED_PASS_SITES, 'Input genotype matrix shape unexpected')
assert(identical(rownames(G), IDS), 'Sample ID order does not match genotype matrix')
assert(identical(as.character(Z$meta$study), rep('1306-PF-NG-NGWA-SM', length(Y))), 'Unexpected study')
assert(all(Z$meta$year == 2020L), 'Unexpected study year')
assert(length(unique(IDS)) == length(IDS), 'Duplicate sample IDs')
assert(all(Y %in% STATES), 'Unexpected state labels')
assert(min(table(factor(Y, STATES))) >= 5L, 'Insufficient samples for five outer folds')
if (any(!is.na(G) & G < 0L)) stop('Invalid negative allele code')
N <- as.matrix(Z$meta[, c('Fws', 'pass_missingness_haploid'), drop=FALSE])
storage.mode(N) <- 'double'
assert(all(is.finite(N)), 'Missing or non-finite Fws/missingness: explicit handling required')
GROUP <- as.character(Z$meta$case_id)
GROUP[is.na(GROUP) | !nzchar(GROUP)] <- IDS[is.na(GROUP) | !nzchar(GROUP)]
assert(!any(vapply(split(Y, GROUP), function(v) length(unique(v)) > 1L, logical(1))),
       'A case group spans multiple states')

# Analysis settings are stored separately; original SETTINGS and scripts remain unchanged.
CFG <- list(seed=20261109L, outer_folds=5L, outer_repeats=5L,
            inner_folds=3L, panel_sizes=c(1L,3L,5L,8L,12L,16L),
            min_call=.90, min_state_call=.80, min_mac=3L,
            ld_r2=.90, min_panel_fraction=.75,
            combined_lambda=c(0,.25,.5,.75,1), tie_tolerance=1e-12)
stopifnot(length(CFG$panel_sizes)>0L, all(CFG$combined_lambda>=0 & CFG$combined_lambda<=1))
RESULT_DIR <- file.path(OUT, 'model_step29_corrected')
FIG_DIR <- file.path(OUT, 'figures_step29_corrected')
AUDIT_DIR <- file.path(OUT, 'audit_step29_corrected')
for(d in c(RESULT_DIR, FIG_DIR, AUDIT_DIR)) dir.create(d, recursive=TRUE, showWarnings=FALSE)
saveRDS(CFG, file.path(AUDIT_DIR,'frozen_step29_config.rds'))
writeLines(capture.output(str(CFG)), file.path(AUDIT_DIR,'frozen_step29_config.txt'))
write_table <- function(df, file, dir=RESULT_DIR) {
  path <- file.path(dir,file)
  write.table(df, path, sep='\t', quote=TRUE, row.names=FALSE, na='NA')
  cat('Saved:', path, '\n')
}

# Stratification is done by STATE; case groups remain intact across all folds.
make_group_folds <- function(y, groups, k, seed) {
  assert(length(y)==length(groups), 'Fold lengths mismatch')
  by_group <- split(seq_along(groups), groups)
  lab <- vapply(by_group, function(ix) {v <- unique(y[ix]); assert(length(v)==1L,'Mixed-state group');v}, character(1))
  assert(all(table(factor(lab,STATES))>=k), 'Too few independent case groups for k folds')
  set.seed(seed)
  fold_group <- setNames(integer(length(by_group)), names(by_group))
  for(s in STATES) {
    gg <- names(by_group)[lab==s]
    gg <- sample(gg, length(gg))
    load <- rep(0L,k)
    for(g in gg) {
      f <- which.min(load)
      fold_group[g] <- f
      load[f] <- load[f] + length(by_group[[g]])
    }
  }
  folds <- as.integer(fold_group[groups])
  assert(!anyNA(folds), 'Missing fold assignments')
  assert(all(table(factor(y,STATES),factor(folds,seq_len(k)))>0L), 'Some folds omit a state')
  folds
}

empty_panel <- function() data.frame(j=integer(), c0=integer(), c1=integer(), chi=numeric())

# Raw VCF allele codes are NOT assumed to be 0/1: 0/2 or 1/2 can be valid.
# Site eligibility, state-wise callability, MAC, ranking and LD pruning are
# computed using INNER/OUTER training samples ONLY.
rank_training <- function(x, y, max_k) {
  n <- nrow(x)
  candidates <- vector('list', ncol(x))
  cnt <- 0L
  for(j in seq_len(ncol(x))) {
    v <- x[,j]
    ok <- !is.na(v)
    nc <- sum(ok)
    if(nc/n < CFG$min_call) next
    codes <- sort(unique(v[ok]))
    if(length(codes)!=2L) next
    tab <- c(sum(v[ok]==codes[1]),sum(v[ok]==codes[2]))
    if(min(tab)<CFG$min_mac) next
    state_n <- vapply(STATES, function(s) sum(ok[y==s]), numeric(1))
    state_den <- vapply(STATES, function(s) sum(y==s), numeric(1))
    if(any(state_n/state_den < CFG$min_state_call)) next
    alt <- vapply(STATES, function(s) sum(v[y==s]==codes[2],na.rm=TRUE), numeric(1))
    p <- tab[2]/nc
    ea <- state_n*p
    er <- state_n*(1-p)
    chi <- sum((alt-ea)^2/pmax(ea,1e-12) + (state_n-alt-er)^2/pmax(er,1e-12))
    cnt <- cnt+1L
    candidates[[cnt]] <- data.frame(j=j,c0=as.integer(codes[1]),c1=as.integer(codes[2]),chi=chi)
  }
  if(!cnt) return(empty_panel())
  df <- do.call(rbind, candidates[seq_len(cnt)])
  df <- df[order(-df$chi,df$j),,drop=FALSE]
  keep <- logical(nrow(df))
  selected <- integer()
  for(i in seq_len(nrow(df))) {
    jj <- df$j[i]; c01 <- c(df$c0[i],df$c1[i])
    red <- FALSE
    for(k in selected) {
      kk <- df$j[k]; c02 <- c(df$c0[k],df$c1[k])
      a <- x[,jj]; b <- x[,kk]
      ok <- !is.na(a) & !is.na(b) & a %in% c01 & b %in% c02
      if(sum(ok) >= 20L) {
        a1 <- as.integer(a[ok]==c01[2]); b1 <- as.integer(b[ok]==c02[2])
        if(length(unique(a1))==2L && length(unique(b1))==2L) {
          r2 <- suppressWarnings(cor(a1,b1)^2)
          if(is.finite(r2) && r2>=CFG$ld_r2) {red <- TRUE; break}
        }
      }
    }
    if(!red) { keep[i] <- TRUE; selected <- c(selected,i) }
    if(length(selected)>=max_k) break
  }
  df[keep,,drop=FALSE]
}

# The fixed positions are POST-HOC selections from this dataset, NOT externally
# determined independent markers. This model is descriptive and cannot validate them.
FIXED4 <- c(409688L,410338L,410817L,413463L)
fixed_training <- function(x) {
  ii <- match(FIXED4,Z$site$pos)
  assert(!anyNA(ii), 'One or more fixed four sites not found in VCF')
  rows <- lapply(ii, function(j) {
    codes <- sort(unique(x[,j][!is.na(x[,j])]))
    if(length(codes)!=2L) return(NULL)
    data.frame(j=j,c0=as.integer(codes[1]),c1=as.integer(codes[2]),chi=NA_real_)
  })
  rows <- Filter(Negate(is.null),rows)
  if(!length(rows)) empty_panel() else do.call(rbind,rows)
}

# If a test sample shows an allele not observed in the training fold, its call
# is treated as unknown, never recoded using held-out information.
score_snp <- function(tr, ytr, te, panel) {
  ans <- matrix(NA_real_,nrow(te),length(STATES),dimnames=list(NULL,STATES))
  complete <- integer(nrow(te)); unexpected <- integer(nrow(te))
  if(!nrow(panel)) return(list(scores=ans,called=complete,unseen=unexpected,panel_n=0L))
  jx <- panel$j
  trcode <- matrix(NA_integer_, nrow(tr), length(jx))
  tecode <- matrix(NA_integer_, nrow(te), length(jx))
  for(h in seq_along(jx)) {
    a <- tr[,jx[h]]; b <- te[,jx[h]]
    trcode[a==panel$c0[h] & !is.na(a),h] <- 0L
    trcode[a==panel$c1[h] & !is.na(a),h] <- 1L
    tecode[b==panel$c0[h] & !is.na(b),h] <- 0L
    tecode[b==panel$c1[h] & !is.na(b),h] <- 1L
    unexpected <- unexpected + as.integer(!is.na(b) & !(b %in% c(panel$c0[h],panel$c1[h])))
  }
  complete <- rowSums(!is.na(tecode))
  for(s in seq_along(STATES)) {
    xx <- trcode[ytr==STATES[s],,drop=FALSE]
    f <- (colSums(xx==1L,na.rm=TRUE)+1)/(colSums(!is.na(xx))+2)
    for(i in seq_len(nrow(tecode))) {
      obs <- which(!is.na(tecode[i,]))
      if(length(obs)<max(1L,ceiling(CFG$min_panel_fraction*nrow(panel)))) next
      ans[i,s] <- mean(ifelse(tecode[i,obs]==1L,log(f[obs]),log(1-f[obs])))
    }
  }
  list(scores=ans,called=complete,unseen=unexpected,panel_n=nrow(panel))
}

score_qc <- function(tr, ytr, te) {
  means <- colMeans(tr); sds <- apply(tr,2,sd)
  sds[!is.finite(sds)|sds<1e-9] <- 1
  tr <- sweep(sweep(tr,2,means,'-'),2,sds,'/')
  te <- sweep(sweep(te,2,means,'-'),2,sds,'/')
  result <- matrix(NA_real_, nrow(te), length(STATES), dimnames=list(NULL,STATES))
  for (s in seq_along(STATES)) {
    ctr <- colMeans(tr[ytr==STATES[s],,drop=FALSE])
    result[,s] <- -rowMeans(sweep(te,2,ctr,'-')^2)/2
  }
  result
}

# Explicit handling of exact/near ties prevents accidental preference for Oyo,
# which happens to be STATES[1]. No calibrated posterior probability is claimed.
decide <- function(scores) {
  out <- rep(NA_character_, nrow(scores))
  tie <- rep(FALSE,nrow(scores))
  for(i in seq_len(nrow(scores))) {
    a <- scores[i,]
    if(any(!is.finite(a))) next
    ord <- order(a,decreasing=TRUE)
    if(a[ord[1]]-a[ord[2]] <= CFG$tie_tolerance) {tie[i] <- TRUE; next}
    out[i] <- STATES[ord[1]]
  }
  list(pred=out,ties=tie)
}

# Normalize the across-class variation in each score vector before blending;
# tunable lambda is selected exclusively by inner CV. These are model scores,
# NOT probabilities. A missing SNP prediction falls back to the QC score.
scale_scores <- function(s) {
  ans <- s
  for(i in seq_len(nrow(s))) {
    z <- s[i,]
    if(any(!is.finite(z))) next
    sdv <- sd(z)
    ans[i,] <- if(sdv<1e-12) rep(0,length(z)) else (z-mean(z))/sdv
  }
  ans
}
blend_scores <- function(snp,qc,lam) {
  sn <- scale_scores(snp); q <- scale_scores(qc)
  out <- q
  for(i in seq_len(nrow(q))) {
    if(all(is.finite(sn[i,]))) out[i,] <- lam*sn[i,] + (1-lam)*q[i,]
  }
  out
}
metrics <- function(truth, pred) {
  rec <- f1 <- numeric(length(STATES))
  for(i in seq_along(STATES)) {
    s <- STATES[i]
    tp <- sum(truth==s & !is.na(pred) & pred==s)
    n <- sum(truth==s)
    fp <- sum(truth!=s & !is.na(pred) & pred==s)
    rec[i] <- if(n)tp/n else NA_real_
    pr <- if(tp+fp)tp/(tp+fp) else 0
    f1[i] <- if(pr+rec[i])2*pr*rec[i]/(pr+rec[i]) else 0
  }
  c(balanced_accuracy=mean(rec),macro_F1=mean(f1),
    overall_accuracy=mean(!is.na(pred) & truth==pred),coverage=mean(!is.na(pred)))
}

# INNER CV hyperparameters: tune K for the SNP-only model and (K, lambda) for
# the combined model on the same inner validation observations.
tune_inner <- function(xt,yt,nt,groups,seed) {
  inner <- make_group_folds(yt,groups,CFG$inner_folds,seed)
  allk <- CFG$panel_sizes
  lam <- CFG$combined_lambda
  oof_snp <- matrix(NA_character_,length(yt),length(allk))
  oof_combo <- array(NA_character_,dim=c(length(yt),length(allk),length(lam)))
  for(f in seq_len(CFG$inner_folds)) {
    v <- which(inner==f);t <- which(inner!=f)
    ranked <- rank_training(xt[t,,drop=FALSE],yt[t],max(allk))
    qcs <- score_qc(nt[t,,drop=FALSE],yt[t],nt[v,,drop=FALSE])
    for(kk in seq_along(allk)) {
      sn <- score_snp(xt[t,,drop=FALSE],yt[t],xt[v,,drop=FALSE],head(ranked,allk[kk]))$scores
      oof_snp[v,kk] <- decide(sn)$pred
      for(ll in seq_along(lam)) {
        oof_combo[v,kk,ll] <- decide(blend_scores(sn,qcs,lam[ll]))$pred
      }
    }
  }
  snp_ba <- vapply(seq_along(allk),function(kk) metrics(yt,oof_snp[,kk])['balanced_accuracy'],numeric(1))
  combo_ba <- matrix(NA_real_, length(allk),length(lam))
  for(kk in seq_along(allk))for(ll in seq_along(lam)) {
    combo_ba[kk,ll] <- metrics(yt,oof_combo[,kk,ll])['balanced_accuracy']
  }
  # Exact ties broken by fewer SNPs; then by smaller SNP weight.
  si <- order(-snp_ba, allk)[1]
  grid <- expand.grid(k_idx=seq_along(allk),l_idx=seq_along(lam))
  grid$ba <- mapply(function(k,l) combo_ba[k,l],grid$k_idx,grid$l_idx)
  grid <- grid[order(-grid$ba, allk[grid$k_idx],lam[grid$l_idx]),,drop=FALSE]
  list(k_snp=allk[si],k_comb=allk[grid$k_idx[1]],lambda=lam[grid$l_idx[1]],
       snp_inner_BA=snp_ba[si],combo_inner_BA=grid$ba[1])
}

run_one <- function(y, groups, seed, repeat_id, produce_detail=TRUE) {
  fold <- make_group_folds(y,groups,CFG$outer_folds,seed)
  model_names <- c('majority','QC_only','SNP_nested','QC_plus_SNP','SNP_fixed4_posthoc')
  p <- matrix(NA_character_,length(y),length(model_names),dimnames=list(NULL,model_names))
  selected <- list(); params <- list();tie_counts <- integer(length(model_names))
  unknown_calls <- integer(length(y))
  for(f in seq_len(CFG$outer_folds)) {
    test <- which(fold==f); train <- which(fold!=f)
    xt <- G[train,,drop=FALSE]; xv <- G[test,,drop=FALSE]
    yt <- y[train]; nt <- N[train,,drop=FALSE]; nv <- N[test,,drop=FALSE]
    tuned <- tune_inner(xt,yt,nt,groups[train],seed+1000L+f)
    ranking <- rank_training(xt,yt,max(c(tuned$k_snp,tuned$k_comb)))
    snp_selected <- head(ranking,tuned$k_snp)
    combined_selected <- head(ranking,tuned$k_comb)
    q_scores <- score_qc(nt,yt,nv)
    s_scores <- score_snp(xt,yt,xv,snp_selected)
    c_scores <- score_snp(xt,yt,xv,combined_selected)
    f_scores <- score_snp(xt,yt,xv,fixed_training(xt))
    scores <- list(QC_only=q_scores, SNP_nested=s_scores$scores,
                   QC_plus_SNP=blend_scores(c_scores$scores,q_scores,tuned$lambda),
                   SNP_fixed4_posthoc=f_scores$scores)
    majority <- names(sort(table(factor(yt,STATES)),decreasing=TRUE))[1]
    p[test,'majority'] <- majority
    for(m in names(scores)) {
      dd <- decide(scores[[m]])
      p[test,m] <- dd$pred
      tie_counts[match(m,model_names)] <- tie_counts[match(m,model_names)]+sum(dd$ties)
    }
    unknown_calls[test] <- s_scores$unseen
    if(produce_detail) {
      selected[[f]] <- if(nrow(snp_selected)) {
        data.frame(repeat_id=rep(repeat_id,nrow(snp_selected)),
          fold=rep(f,nrow(snp_selected)),
          pos=Z$site$pos[snp_selected$j], original_code0=snp_selected$c0,
          original_code1=snp_selected$c1, stringsAsFactors=FALSE)
      } else data.frame(repeat_id=integer(),fold=integer(),pos=integer(),
                        original_code0=integer(),original_code1=integer())
      params[[f]] <- data.frame(repeat_id=repeat_id,fold=f,n_test=length(test),
          k_snp_requested=tuned$k_snp,k_snp_available=nrow(snp_selected),
          k_combo_requested=tuned$k_comb,k_combo_available=nrow(combined_selected),
          combo_weight_SNP=tuned$lambda,inner_snp_BA=tuned$snp_inner_BA,
          inner_combo_BA=tuned$combo_inner_BA,
          n_SNP_abstained=sum(is.na(p[test,'SNP_nested'])),
          n_combo_abstained=sum(is.na(p[test,'QC_plus_SNP'])),
          n_test_unseen_alleles=sum(s_scores$unseen),
          n_SNP_ties=sum(decide(s_scores$scores)$ties),
          stringsAsFactors=FALSE)
    }
  }
  if(!produce_detail) return(list(metrics=do.call(rbind,lapply(colnames(p),function(m)metrics(y,p[,m])))))
  predictions <- data.frame(repeat_id=repeat_id,sample_id=IDS,truth=y,outer_fold=fold,
                            p,unknown_allele_calls=unknown_calls,check.names=FALSE)
  res_metrics <- do.call(rbind,lapply(colnames(p),function(m){
    a <- metrics(y,p[,m]);data.frame(repeat_id=repeat_id,model=m,t(a))
  }))
  list(metrics=res_metrics,predictions=predictions,
       selected=do.call(rbind,selected),params=do.call(rbind,params),ties=tie_counts)
}

# The five repetitions represent re-partitions of the SAME 119 isolates.
all_results <- vector('list',CFG$outer_repeats)
for(r in seq_len(CFG$outer_repeats)) {
  all_results[[r]] <- run_one(Y,GROUP,seed=CFG$seed+100L*r,repeat_id=r)
  m <- all_results[[r]]$metrics
  sn <- m$balanced_accuracy[m$model=='SNP_nested']
  qc <- m$balanced_accuracy[m$model=='QC_only']
  co <- m$balanced_accuracy[m$model=='QC_plus_SNP']
  cat(sprintf('Repeat %d/%d: SNP %.3f | QC %.3f | combined %.3f\n',r,CFG$outer_repeats,sn,qc,co))
}
DF <- do.call(rbind,lapply(all_results,`[[`,'predictions'))
MET <- do.call(rbind,lapply(all_results,`[[`,'metrics'))
SEL <- do.call(rbind,lapply(all_results,`[[`,'selected'))
PARAM <- do.call(rbind,lapply(all_results,`[[`,'params'))
assert(nrow(DF)==length(Y)*CFG$outer_repeats, 'Unexpected saved prediction count')
assert(!anyDuplicated(DF[c('repeat_id','sample_id')]), 'Duplicate sample/repeat prediction')
assert(all(vapply(split(DF$sample_id,DF$repeat_id), function(x) setequal(x,IDS), logical(1))),
       'Not all isolates included in every repeat')
write_table(DF,'Table_29A_fold_aligned_out_of_fold_predictions.tsv')
write_table(MET,'Table_29B_metrics_by_repeat_and_model.tsv')
# Secondary, explicitly selection-conditional analysis: QC and SNP scored on
# exactly the subset to which the SNP model assigned a label. NOT the headline.
matched_evaluable <- do.call(rbind,lapply(sort(unique(DF$repeat_id)),function(r) {
  z <- DF[DF$repeat_id==r,,drop=FALSE]
  ok <- !is.na(z$SNP_nested) & !is.na(z$QC_only)
  a <- if(all(STATES %in% z$truth[ok])) metrics(z$truth[ok],z$SNP_nested[ok]) else
    setNames(rep(NA_real_,4), c('balanced_accuracy','macro_F1','overall_accuracy','coverage'))
  b <- if(all(STATES %in% z$truth[ok])) metrics(z$truth[ok],z$QC_only[ok]) else
    setNames(rep(NA_real_,4), c('balanced_accuracy','macro_F1','overall_accuracy','coverage'))
  data.frame(repeat_id=r,n_jointly_evaluable=sum(ok),
             subset_SNP_balanced_accuracy=unname(a['balanced_accuracy']),
             subset_QC_balanced_accuracy=unname(b['balanced_accuracy']),
             subset_delta_SNP_minus_QC=unname(a['balanced_accuracy']-b['balanced_accuracy']))
}))
write_table(matched_evaluable,'Table_29M_jointly_evaluable_secondary_comparison.tsv')
write_table(PARAM,'Table_29C_nested_hyperparameters_by_outer_fold.tsv')
write_table(SEL,'Table_29D_outer_training_marker_selections.tsv')
if(nrow(SEL)) {
  s <- as.data.frame(table(SEL$pos), stringsAsFactors=FALSE)
  names(s) <- c('pos','n_outer_fold_selections');s <- s[s$n_outer_fold_selections>0,,drop=FALSE]
  s$pos <- as.integer(as.character(s$pos))
  s$selection_fraction <- s$n_outer_fold_selections/(CFG$outer_repeats*CFG$outer_folds)
  write_table(s[order(-s$selection_fraction,s$pos),],'Table_29E_marker_selection_stability.tsv')
}

tab <- aggregate(MET[,c('balanced_accuracy','macro_F1','overall_accuracy','coverage')],
                 by=list(model=MET$model),FUN=mean)
write_table(tab,'Table_29F_repeated_mean_performance.tsv')
# Compare held-out test data, class membership, and abstentions on identical folds.
M <- reshape(MET[,c('repeat_id','model','balanced_accuracy','macro_F1','coverage')],
             idvar='repeat_id',timevar='model',direction='wide')
for(other in c('SNP_nested','QC_plus_SNP','SNP_fixed4_posthoc')) {
  for(metric in c('balanced_accuracy','macro_F1','coverage')) {
    M[[paste0('delta_',other,'_minus_QC_',metric)]] <-
      M[[paste0(metric,'.',other)]] - M[[paste0(metric,'.QC_only')]]
  }
}
write_table(M,'Table_29G_fold_aligned_paired_descriptive_deltas.tsv')
first <- DF[DF$repeat_id==1L,,drop=FALSE]
for(model in c('QC_only','SNP_nested','QC_plus_SNP','SNP_fixed4_posthoc')) {
  pred <- first[[model]]
  pred[is.na(pred)] <- 'ABSTAIN'
  cm <- as.data.frame.matrix(table(factor(first$truth,levels=STATES),
                                   factor(pred,levels=c(STATES,'ABSTAIN'))))
  cm <- cbind(truth=rownames(cm),cm)
  write_table(cm,paste0('Table_29H_first_repeat_confusion_',model,'.tsv'))
}
state_summary <- do.call(rbind,lapply(unique(MET$model),function(model) {
  do.call(rbind,lapply(STATES,function(state) {
    rr <- DF[DF$truth==state,,drop=FALSE]
    byrepeat <- split(rr,rr$repeat_id)
    recall <- vapply(byrepeat,function(z)mean(!is.na(z[[model]]) & z[[model]]==state),numeric(1))
    coverage <- vapply(byrepeat,function(z)mean(!is.na(z[[model]])),numeric(1))
    data.frame(model=model,state=state,n_unique_isolates=sum(Y==state),
      mean_repeated_recall=mean(recall),min_repeat_recall=min(recall),
      max_repeat_recall=max(recall),mean_coverage=mean(coverage))
  }))
}))
write_table(state_summary,'Table_29I_state_recall_and_coverage.tsv')

# Keep null hypothesis separate: label independence is NOT incremental predictive
# value beyond QC. Permutations optional, and re-fit entire nested workflow.
perm_n <- suppressWarnings(as.integer(Sys.getenv('PI4K_STEP29_PERM','0')))
assert(is.finite(perm_n) && perm_n>=0L, 'PI4K_STEP29_PERM must be a nonnegative integer')
if(perm_n>0L) {
  assert(length(unique(GROUP))==length(GROUP),
    'Label permutations require independent cases in this implementation')
  null <- numeric(perm_n)
  for(i in seq_len(perm_n)) {
    set.seed(CFG$seed+99999L+i)
    yp <- sample(Y)
    out <- run_one(yp,GROUP,seed=CFG$seed+100L,repeat_id=1L,produce_detail=FALSE)
    # Rows in return_detail=FALSE are in model order, so take SNP_nested (third).
    null[i] <- out$metrics[3,'balanced_accuracy']
    if(i%%10L==0L)cat(sprintf('Optional label-null rerun %d/%d\n',i,perm_n))
  }
  obs <- subset(MET,repeat_id==1L & model=='SNP_nested')$balanced_accuracy
  write_table(data.frame(null_balanced_accuracy=null),'Table_29J_label_independence_null.tsv')
  write_table(data.frame(observed_first_repeat_BA=obs,n_perm=perm_n,
    p_label_independence=(1+sum(null>=obs))/(perm_n+1),
    important_note='NOT a test of additional signal beyond Fws/quality'),
    'Table_29K_label_independence_permutation_summary.tsv')
}

writeLines(capture.output(sessionInfo()),file.path(AUDIT_DIR,'R_sessionInfo.txt'))
writeLines(c('Exploratory cross-validation only; no external data.',
    'Same-study 2020 cohort; geographic labels mean collection state, not infection origin.',
    'SNP ranking, MAC, per-state callability, LD pruning, allele-code mapping and panel tuning use training samples only.',
    'Outer folds are shared by majority, QC, SNP, combined and post-hoc fixed-four comparators.',
    'Abstentions are counted as failures in balanced accuracy, macro-F1 and overall accuracy.',
    'Exact score ties abstain; no state receives artificial preference from label order.',
    'Combined model may fall back to QC when the SNP panel cannot score a sample.',
    'The fixed four markers were selected from the same cohort: not independent validation.',
    'Repeated CV folds are correlated; no t tests on repeats and no independent confidence intervals.',
    'The optional shuffled-label test does NOT quantify genetic added value beyond QC.',
    'Unobserved held-out allele codes are recorded and treated as missing, not silently as REF.'),
  file.path(AUDIT_DIR,'interpretation_guardrails.txt'))

if(requireNamespace('ggplot2',quietly=TRUE)) {
  library(ggplot2)
  p <- ggplot(MET,aes(x=factor(repeat_id),y=balanced_accuracy,color=model,group=model))+
    geom_hline(yintercept=0.2,linetype='dashed',colour='grey65')+
    geom_line(linewidth=.65)+geom_point(size=2)+
    scale_y_continuous(limits=c(0,1),labels=function(x)sprintf('%d%%',round(x*100)))+
    labs(title='PI4K genetic and QC model performance, matched folds',
      subtitle='Five repeated nested CV partitions of the same 119 isolates',
      x='Outer fold repeat',y='Balanced accuracy',color='Model')+
    clean_theme()
  export_one <- function(p,stem,w=9,h=5) {
    ggplot2::ggsave(file.path(FIG_DIR,paste0(stem,'.pdf')),plot=p,width=w,height=h,units='in',device='pdf',useDingbats=FALSE)
    ggplot2::ggsave(file.path(FIG_DIR,paste0(stem,'.tiff')),plot=p,width=w,height=h,units='in',dpi=600,compression='lzw')
  }
  export_one(p,'Figure_29A_matched_fold_balanced_accuracy')
  dd <- data.frame(repeat_id=M$repeat_id,
                   snp=M$delta_SNP_nested_minus_QC_balanced_accuracy,
                   combined=M$delta_QC_plus_SNP_minus_QC_balanced_accuracy)
  long <- rbind(data.frame(repeat_id=dd$repeat_id,model='SNP minus QC',delta=dd$snp),
                data.frame(repeat_id=dd$repeat_id,model='Combined minus QC',delta=dd$combined))
  p2 <- ggplot(long,aes(x=factor(repeat_id),y=100*delta,group=model,color=model))+
    geom_hline(yintercept=0,linetype='dashed',colour='grey60')+
    geom_point(size=2.3)+geom_line(linewidth=.65)+
    labs(title='Paired descriptive change in held-out balanced accuracy',
      subtitle='Matched folds; not five independent studies; percentages are point differences',
      x='Outer fold repeat',y='Balanced-accuracy difference (percentage points)',color='Comparison')+
    clean_theme()
  export_one(p2,'Figure_29B_model_minus_QC_differences')
  p3 <- ggplot(state_summary,aes(x=state,y=mean_repeated_recall,fill=model))+
    geom_col(position='dodge')+
    scale_y_continuous(limits=c(0,1),labels=function(x)sprintf('%d%%',round(x*100)))+
    labs(title='Class-specific recall, repeated held-out predictions',
      subtitle='Abstentions count as incorrect; means over repeated splits of same isolates',
      x='Collection state',y='Mean recall',fill='Model')+
    clean_theme()
  export_one(p3,'Figure_29C_per_state_recall')
} else cat('ggplot2 not installed: numeric tables generated, figures skipped.\n')
cat('STEP 29 CORRECTED COMPLETE. Input files and earlier outputs untouched.\n')
