# General statistical methodology

This document describes the kinds of methods implemented in the scripts; it contains no empirical estimates, state associations, allele positions or study results.

## Genotype quality and allele encoding

Calls may be represented in more than one encoded slot. For an effectively haploid biological target, concordant calls may be represented as one observed allele; discordant calls require an explicit missing/ambiguity policy. Site and sample callability must be measured separately. A variant site's PASS flag does not establish individual call quality.

## Allele frequencies and uncertainty

Within a group, the frequency of a named allele is `x / n`, where `x` is the number of isolates bearing the allele among `n` callable isolates. Retain monomorphic groups in the output, and state the denominator for every comparison. Wilson intervals describe sampling uncertainty for a binomial frequency under their assumptions.

## Within-population diversity

For a biallelic haploid locus with allele frequency `p`, the standard expected gene diversity is `2p(1-p)`. Unbiased sample corrections are available when the callable sample size is sufficient. Diversity over a selected-SNP panel is not equivalent to genome-wide heterozygosity or nucleotide diversity across a fully callable gene sequence.

## Genetic differentiation

Pairwise Hudson-type differentiation statistics compare within- and between-population diversity. Numerical outcomes depend on the weighting and site policy. Retain signed estimates, document pair-informative-site counts, and distinguish estimates based on one linked genomic region from genome-wide results. Bootstrap intervals and permutation tests require their own clear resampling units.

## LD and multilocus configurations

Pairwise squared correlation (`r^2`) estimates association between encoded variants where joint calls exist. Low or high `r^2` within one gene must not be interpreted as independent genomic replication. Genotype combinations derived from haploid consensus do not establish phased haplotypes in mixed infection data.

## Prediction and validation

Evaluate feature selection, imputation and hyperparameter tuning within the training partitions of nested CV. Use matched test folds for model comparison. Report class-specific recall, macro-F1, balanced accuracy, overall accuracy and prediction coverage, explicitly treating abstentions. Repeated CV predictions from the same individual are correlated and cannot be treated as independent participants. A fixed-model/bootstrap analysis conditional on saved predictions cannot replace external validation.

A statistical association does not establish a diagnostic or geographic barcode. Laboratory and independent-population validation must be planned separately.
