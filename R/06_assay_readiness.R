source('R/01_helpers.R')
cat('WP7: Candidate flanks and gene-structure overlap. NO primer validation inferred.\n')
candidate_path<-file.path(OUT,'assay','Candidate_assay_SNPs_NOT_VALIDATED.tsv')
if(!file.exists(candidate_path)) {
  cat('No ranked candidate SNPs produced; WP7 marker feasibility deferred.\n')
} else {
CAND<-read_tsv(candidate_path)
# Parse the public reference genome locally. Verify REF against 3D7 reference for every selected site.
fa<-readLines(required(FASTA),warn=FALSE)
h<-grep('^>',fa)
target<-which(grepl(paste0('^>',REGION_CHROM,'([[:space:]]|$)'),fa[h]))
assert(length(target)==1L,'Cannot uniquely identify Pf3D7 chromosome-5 FASTA header')
start<-h[target]+1L; end<-if(target<length(h))h[target+1L]-1L else length(fa)
refseq<-paste0(fa[start:end],collapse='')
assert(nchar(refseq)>=REGION_END+250,'Chromosome sequence shorter than expected')
# Only gene-structure overlap, NOT codon consequence or functional annotation.
gff<-readLines(gzfile(required(GFF),'rt'),warn=FALSE)
gff<-gff[!startsWith(gff,'#')]
gff<-gff[grepl('PF3D7_0509800',gff,fixed=TRUE) & grepl(REGION_CHROM,gff,fixed=TRUE)]
gff_df<-if(length(gff)) {
 z<-read.delim(text=gff,header=FALSE,sep='\t',quote='',comment.char='',stringsAsFactors=FALSE)
 if(ncol(z)>=9)z else NULL
} else NULL
out<-list()
for(i in seq_len(nrow(CAND))) {
 pos<-as.integer(CAND$pos[i]); expected<-toupper(CAND$ref[i]); observed<-toupper(substr(refseq,pos,pos+nchar(expected)-1L))
 assert(identical(expected,observed),paste('3D7 REF mismatch at',pos,'expected',expected,'observed',observed))
 left<-max(1L,pos-250L);right<-min(nchar(refseq),pos+250L)
 struct<-if(!is.null(gff_df)) {
   feature<-gff_df[gff_df[[4]]<=pos & gff_df[[5]]>=pos,3]
   paste(unique(feature),collapse=';')
 } else ''
 out[[i]]<-data.frame(pos=pos,ref=CAND$ref[i],alt=CAND$alt[i],
    flank_start=left,flank_end=right,
    reference_501bp=substr(refseq,left,right),
    gene_annotation_feature_overlap=struct,
    assay_status='NOT_SCREENED_FOR_PRIMERS_OFFTARGETS_INDELS_OR_CLINICAL_VALIDATION')
}
write_tsv(do.call(rbind,out),'Candidate_SNP_3D7_flanks_and_GFF_overlaps.tsv','assay')
cat('WP7 preliminary reference flanks generated. No validated primer design or functional consequence assigned.\n')

}
