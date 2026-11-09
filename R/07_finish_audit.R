source('R/01_helpers.R')
cat('WP8: Versioning, output audit, independent validation protocol template\n')
files<-list.files(OUT,recursive=TRUE,full.names=TRUE)
files<-files[file.info(files)$isdir==FALSE]
manifest<-data.frame(relative_path=substring(files,nchar(OUT)+2L),bytes=file.info(files)$size,
  sha256=vapply(files,function(f)strsplit(system2('sha256sum',f,stdout=TRUE),' +')[[1]][1],character(1)))
write_tsv(manifest,'output_checksums.tsv','audit')
writeLines(capture.output(sessionInfo()),file.path(OUT,'audit','R_sessionInfo.txt'))
external<-data.frame(sample_id=character(),state_of_collection=character(),
  location_of_acquisition_evidence=character(),year=integer(),study_site=character(),
  age_group=character(),travel_history=character(),assay_batch=character(),
  genotype_call_qc_status=character(),is_independent_subject=logical())
write_tsv(external,'EXTERNAL_VALIDATION_DATA_TEMPLATE.tsv','audit')
cat('Versioned outputs audited. External validation remains prospective and NOT COMPLETE.\n')
