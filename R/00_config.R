# PI4K Southwest Nigeria — Analysis configuration. Freeze before analysing state associations.
options(stringsAsFactors = FALSE)
ROOT <- Sys.getenv('PI4K_ROOT', unset='/mnt/c/PI4K_Nigeria_Study')
stopifnot(dir.exists(ROOT))
OUT <- file.path(ROOT, 'results', 'Southwest5', 'R_Analysis_v1')
for (x in c('tables','figures','rds','audit','model','assay')) dir.create(file.path(OUT,x), recursive=TRUE, showWarnings=FALSE)
STATES <- c('Oyo','Ogun','Lagos','Osun','Ondo')
FILE_ALL <- file.path(ROOT, 'variants/Southwest5/Sanger_SNP_complete/PI4K_Southwest5_282_Sanger_SNP_PASS_derived.vcf.gz')
FILE_2020 <- file.path(ROOT, 'variants/Southwest5/Sanger_SNP_complete/Step21_2020/PI4K_2020_119_PASS_sites.vcf.gz')
META <- file.path(ROOT, 'results/Southwest5/Southwest5_QC_pass_samples.tsv')
FWS <- file.path(ROOT, 'metadata/Pf8_fws.tsv')
SAMPLE_QC <- file.path(ROOT, 'results/Southwest5/Step18_QC/PI4K_QC_by_sample.tsv')
GFF <- file.path(ROOT, 'reference/PlasmoDB-55_Pfalciparum3D7.gff.gz')
FASTA <- file.path(ROOT, 'reference/PlasmoDB-54-Pfalciparum3D7-Genome.fasta')
SAMPLE_LIST_DIR <- file.path(ROOT, 'results/Southwest5/sample_lists')
COHORT_VCF_DIR <- file.path(ROOT, 'variants/Southwest5/Sanger_SNP_complete/Step21_2020')
REGION_CHROM <- 'Pf3D7_05_v3'; REGION_START <- 409403L; REGION_END <- 415758L
EXPECTED_ALL <- 282L; EXPECTED_2020 <- 119L; EXPECTED_PASS_SITES <- 1824L
SETTINGS <- list(seed=20261109L, call_rate=0.90, state_call_rate=0.80,
  minor_allele_count=3L, fws_low_diversity=0.95, sample_missingness_limit=0.10,
  bootstrap_reps=400L, fst_perm_reps=499L, global_perm_reps=999L,
  outer_folds=5L, outer_repeats=5L, inner_folds=3L,
  panel_sizes=c(1L,3L,5L,8L,12L,16L), barcode_perm_reps=99L,
  ld_prune_r2=0.90, abstain_margin=0.00)
set.seed(SETTINGS$seed)
