# STEP 25 - Audited PI4K state differentiation: five Southwest Nigeria states.
# Run from the PI4K_Southwest_R_Analysis_Suite root directory.
# Uses only previously verified, locally cached RDS objects. No network requests.
source('R/01_helpers.R')
cat('STEP 25: Audited PI4K geographic differentiation (OFFLINE)\n')
cohorts <- readRDS(required(file.path(OUT, 'rds', 'QC_cohorts.rds')))
assert(all(c('primary2020', 'all282') %in% names(cohorts)), 'Cohort RDS missing')
assert(ncol(cohorts$primary2020$X$G) == 119L, '2020 N mismatch')
assert(ncol(cohorts$all282$X$G) == 282L, 'Full cohort N mismatch')
set.seed(SETTINGS$seed + 25L)

# Recode the two distinct observed consensus allele codes at each selected SNP.
# This handles pooled allele pairs such as 0/2 or 1/2 (not just 0/1).
make_panel <- function(z, mode) {
  G <- z$X$G
  select <- if (mode == 'strict') which(z$site$informative) else which(z$site$allele_count == 2L)
  if (!length(select)) stop('No selected sites for ', mode)
  M <- matrix(NA_real_, nrow=length(select), ncol=ncol(G),
              dimnames=list(as.character(z$site$pos[select]), colnames(G)))
  allele_codes <- vector('list', length(select))
  for (j in seq_along(select)) {
    x <- G[select[j], ]
    codes <- sort(unique(as.integer(x[!is.na(x)])))
    assert(length(codes) == 2L, paste('Site is not consensus-biallelic:', z$site$pos[select[j]]))
    M[j, ] <- ifelse(is.na(x), NA_real_, as.numeric(x == codes[2]))
    allele_codes[[j]] <- codes
  }
  mapping <- data.frame(pos_1based=z$site$pos[select], ref=z$site$ref[select],
     vcf_ALT=z$site$alt[select],
     allele_code_0=vapply(allele_codes, `[[`, integer(1), 1L),
     allele_code_1=vapply(allele_codes, `[[`, integer(1), 2L),
     stringsAsFactors=FALSE)
  assert(all(is.na(M) | M %in% c(0,1)), 'Nonbinary consensus genotype after recoding')
  list(G=M, mapping=mapping)
}

# Ratio-of-sums Hudson-style estimator. Preserve negative estimates and NA values.
# SNPs in this one gene are linked: isolate bootstrap reflects conditional sampling
# uncertainty and MUST NOT be called genome-wide or locus-independent uncertainty.
hudson_ratio <- function(A, B) {
  nA <- rowSums(!is.na(A)); nB <- rowSums(!is.na(B))
  pA <- rowSums(A == 1, na.rm=TRUE)/pmax(nA, 1L)
  pB <- rowSums(B == 1, na.rm=TRUE)/pmax(nB, 1L)
  den <- pA*(1-pB) + pB*(1-pA)
  valid <- nA>=2L & nB>=2L & is.finite(den) & den>0
  if (!any(valid)) return(c(fst=NA_real_, sites=0))
  num <- (pA[valid]-pB[valid])^2 -
    pA[valid]*(1-pA[valid])/(nA[valid]-1L) -
    pB[valid]*(1-pB[valid])/(nB[valid]-1L)
  c(fst=sum(num)/sum(den[valid]), sites=sum(valid))
}

pairwise_audit <- function(panel, states, cohort_name, marker_policy,
                           nboot=499L, nperm=999L) {
  rows <- list()
  for (pair in combn(STATES, 2L, simplify=FALSE)) {
    iA <- which(states == pair[1L]); iB <- which(states == pair[2L])
    A <- panel[, iA, drop=FALSE]; B <- panel[, iB, drop=FALSE]
    est <- hudson_ratio(A, B)
    boots <- rep(NA_real_, nboot)
    for (r in seq_len(nboot)) {
      aa <- A[,sample(seq_len(ncol(A)), ncol(A), replace=TRUE),drop=FALSE]
      bb <- B[,sample(seq_len(ncol(B)), ncol(B), replace=TRUE),drop=FALSE]
      boots[r] <- unname(hudson_ratio(aa,bb)['fst'])
    }
    combined <- cbind(A,B); nA <- ncol(A); nTotal <- ncol(combined)
    perm <- rep(NA_real_, nperm)
    for (r in seq_len(nperm)) {
      ix <- sample.int(nTotal)
      perm[r] <- unname(hudson_ratio(
        combined[,ix[seq_len(nA)],drop=FALSE],
        combined[,ix[(nA+1L):nTotal],drop=FALSE])['fst'])
    }
    valid_boot <- boots[is.finite(boots)]
    valid_perm <- perm[is.finite(perm)]
    p_one_sided <- if (is.finite(est['fst']) && length(valid_perm)>0L)
      (1+sum(valid_perm>=est['fst']))/(1+length(valid_perm)) else NA_real_
    ci <- if (length(valid_boot)>=max(30L, as.integer(.75*nboot)))
      unname(quantile(valid_boot, c(.025,.975), na.rm=TRUE)) else c(NA_real_,NA_real_)
    rows[[length(rows)+1L]] <- data.frame(
      cohort=cohort_name, marker_policy=marker_policy,
      state_a=pair[1L], state_b=pair[2L], n_a=ncol(A), n_b=ncol(B),
      n_selected_sites=nrow(panel), n_pair_informative_sites=as.integer(est['sites']),
      hudson_Fst_signed=unname(est['fst']),
      isolate_bootstrap_CI_low=ci[1L], isolate_bootstrap_CI_high=ci[2L],
      bootstrap_valid_draws=length(valid_boot), permutation_valid_draws=length(valid_perm),
      p_one_sided_positive_differentiation=p_one_sided,
      uncertainty_note='isolates resampled; PI4K SNPs are linked',
      stringsAsFactors=FALSE)
  }
  out <- do.call(rbind, rows)
  out$BH_q_10_pairs <- p.adjust(out$p_one_sided_positive_differentiation, method='BH')
  out$interpretation_flag <- ifelse(out$n_pair_informative_sites<2L, 'INSUFFICIENT_PAIRWISE_SITES',
      ifelse(out$n_pair_informative_sites<5L,'VERY_FEW_LINKED_SITES','EXPLORATORY_ONE_GENE'))
  out
}

# Four candidate SNPs: 5x2 exact Monte Carlo test with 9999 replicates;
# report allele code mapping, state denominators and BH adjustment across four sites.
marker_associations <- function(panel, states, mapping) {
  rr <- vector('list', nrow(panel))
  for(i in seq_len(nrow(panel))) {
    x<-panel[i,]; ok<-!is.na(x)
    tab<-table(factor(states[ok],levels=STATES), factor(x[ok],levels=c(0,1)))
    assert(all(rowSums(tab)>0L), paste('No calls in a state at',mapping$pos_1based[i]))
    ft<-suppressWarnings(fisher.test(tab,simulate.p.value=TRUE,B=9999L))
    rr[[i]]<-data.frame(pos_1based=mapping$pos_1based[i],ref=mapping$ref[i],
      vcf_ALT=mapping$vcf_ALT[i],observed_lower_allele_code=mapping$allele_code_0[i],
      observed_higher_allele_code=mapping$allele_code_1[i],
      n_called=sum(ok), n_higher_allele=sum(x==1,na.rm=TRUE),
      p_Fisher_montecarlo_5x2=unname(ft$p.value),
      stringsAsFactors=FALSE)
  }
  out <- do.call(rbind,rr)
  out$BH_q_four_markers <- p.adjust(out$p_Fisher_montecarlo_5x2,method='BH')
  out
}

# Explicit primary and sensitivity definitions: do not silently replace four SNPs
# with 43 SNPs in the principal differentiation table.
analyses <- list(
  list(label='2020_same_study',policy='strict',nboot=499L,nperm=999L),
  list(label='2020_same_study',policy='consensus_biallelic_exploratory',nboot=299L,nperm=499L),
  list(label='all282_provenance_uncontrolled',policy='strict',nboot=299L,nperm=499L)
)
result <- list()
for (cfg in analyses) {
  z <- if (cfg$label=='2020_same_study') cohorts$primary2020 else cohorts$all282
  mode <- if (cfg$policy=='strict') 'strict' else 'exploratory'
  p <- make_panel(z,mode)
  cat(sprintf('%s / %s: %d selected sites, %d samples\n',cfg$label,cfg$policy,nrow(p$G),ncol(p$G)))
  if(cfg$label=='2020_same_study' && cfg$policy=='strict') {
    assert(nrow(p$G)==4L,'Expected four strict SNPs; stop and audit if changed')
    assoc<-marker_associations(p$G,z$meta$state,p$mapping)
    write_tsv(assoc,'Table_S28_four_marker_global_state_exact_tests.tsv')
    write_tsv(p$mapping,'Table_S29_audited_marker_encoding.tsv')
  }
  q <- pairwise_audit(p$G,z$meta$state,cfg$label,cfg$policy,
                      cfg$nboot,cfg$nperm)
  result[[length(result)+1L]] <- q
}
Q <- do.call(rbind,result)
primary <- Q[Q$cohort=='2020_same_study' & Q$marker_policy=='strict',,drop=FALSE]
assert(nrow(primary)==10L,'Expected 10 primary pairwise state contrasts')
write_tsv(primary,'Table_4_primary2020_four_marker_Hudson_Fst_AUDITED.tsv')
write_tsv(Q,'Table_S26_PI4K_Fst_primary_and_sensitivity_AUDITED.tsv')
write_tsv(Q[,c('cohort','marker_policy','state_a','state_b','n_a','n_b',
 'n_selected_sites','n_pair_informative_sites','interpretation_flag')],
 'Table_S27_Fst_pairwise_information_audit.tsv')

if (requireNamespace('ggplot2', quietly=TRUE)) {
  library(ggplot2)
  p <- primary
  p$contrast <- factor(paste(p$state_a,p$state_b,sep=' vs '),
       levels=rev(paste(p$state_a,p$state_b,sep=' vs ')))
  f <- ggplot(p,aes(hudson_Fst_signed,contrast)) +
    geom_vline(xintercept=0,linetype='dashed',color='grey45') +
    geom_segment(aes(x=isolate_bootstrap_CI_low,xend=isolate_bootstrap_CI_high,
                     y=contrast,yend=contrast),color='#52708B',linewidth=.65) +
    geom_point(size=2.2,color='#A45541')+
    labs(title='PI4K geographic differentiation: common-study 2020 cohort',
      subtitle='Four linked SNPs; isolate-bootstrap 95% intervals; exploratory estimates',
      x='Hudson-type FST (signed; no zero truncation)', y=NULL)+clean_theme()
  export_plot(f,'Figure_7a_primary2020_four_marker_Fst_AUDITED',9.5,6.2)
}
writeLines(c('STEP25 COMPLETE - no genomic download or genotype modification',
  'Primary: same-study 2020 cohort / four prespecified high-quality SNPs.',
  'Sensitivity: 2020 consensus-biallelic sites, full-282 strict sites.',
  'SNPs come from a single short linked locus, so confidence intervals are conditional on this locus.',
  'Permutation p-values are exploratory; study-year confounding remains for the 282 cohort.',
  'No validated five-state diagnostic barcode has been established.'),
  file.path(OUT,'audit','Step25_interpretation_readme.txt'))
cat('STEP 25 COMPLETE: audited FST tables, four-marker Fisher tests, optional Figure 7a.\n')
