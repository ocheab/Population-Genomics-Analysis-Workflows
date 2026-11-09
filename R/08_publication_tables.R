source('R/01_helpers.R')
if(!requireNamespace('flextable',quietly=TRUE)||!requireNamespace('officer',quietly=TRUE)) {
  message('Publication Word tables skipped; install.packages(c("flextable","officer")) and rerun R/08_publication_tables.R')
} else {
  library(flextable);library(officer)
  specs<-list(
    list(title='Table 1. PI4K study population and quality metrics',file='tables/Table_1_study_population_and_QC.tsv',
      note='Same-study 2020 subset is the main geographic comparison; Fws and missingness are descriptive.'),
    list(title='Table 2. PI4K polymorphism and site filtering audit',file='tables/Table_2_site_attrition.tsv',
      note='Informativeness is defined on unambiguous haploid consensus calls, not source two-slot counts.'),
    list(title='Table 3. PI4K within-state diversity',file='tables/Table_3_within_PI4K_allelic_diversity.tsv',
      note='Site heterozygosity is estimated over assayed informative positions; not an unbiased genome-wide nucleotide diversity estimate.'),
    list(title='Table 4. Pairwise Hudson-type PI4K FST',file='tables/Table_4_primary2020_pairwise_Fst.tsv',
      note='Signed estimates with isolate resampling and pairwise label-permutation tests; a single linked locus limits generalisation.'),
    list(title='Table 5. Exploratory nested-CV PI4K barcode evaluation',file='model/Table_5_barcode_internal_validation_summary.tsv',
      note='Internal discovery only; panel is NOT externally or experimentally validated. Preliminary permutation counts must be increased for confirmatory claims.')
  )
  doc<-officer::read_docx()
  doc<-officer::body_add_par(doc,'PI4K Southwestern Nigeria — Manuscript-ready statistical tables',style='heading 1')
  doc<-officer::body_add_par(doc,'Analysis outputs generated in R from locally verified Pf8 SNP-only calls. All results require investigator scientific review.',style='Normal')
  for(z in specs) {
    p<-file.path(OUT,z$file)
    if(!file.exists(p)){message('Missing table:',p);next}
    d<-read_tsv(p)
    # Preserve raw precision in TSVs; format numerical values only for presentation.
    for(nm in names(d)) if(is.numeric(d[[nm]])&&!all(is.na(d[[nm]]))) {
      if(any(is.finite(d[[nm]]))) d[[nm]]<-formatC(d[[nm]],format='fg',digits=4)
    }
    doc<-officer::body_add_par(doc,z$title,style='heading 2')
    f<-flextable::flextable(d)
    f<-flextable::theme_booktabs(f)
    f<-flextable::bold(f,part='header')
    f<-flextable::fontsize(f,size=8,part='all')
    f<-flextable::autofit(f)
    doc<-flextable::body_add_flextable(doc,value=f)
    doc<-officer::body_add_par(doc,z$note,style='Normal')
  }
  path<-file.path(OUT,'tables','PI4K_manuscript_tables_R_generated.docx')
  print(doc,target=path)
  cat('Publication table document exported:',path,'\n')
}
