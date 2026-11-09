# Core utilities; base R for calculations. ggplot2 is used for figures only.
source('R/00_config.R')
assert <- function(ok, msg) if (!isTRUE(ok)) stop(msg, call.=FALSE)
required <- function(path) {assert(file.exists(path), paste('Missing expected file:', path)); path}
write_tsv <- function(x, name, subfolder='tables') {
  p <- file.path(OUT,subfolder,name)
  write.table(x,p,sep='\t',row.names=FALSE,quote=TRUE,na='')
  p
}
read_tsv <- function(p) read.delim(required(p),check.names=FALSE,fileEncoding='UTF-8-BOM',quote='"',stringsAsFactors=FALSE)
read_gt_vcf <- function(path) {
  path <- required(path)
  lines <- readLines(gzfile(path,'rt'), warn=FALSE)
  header <- grep('^#CHROM\t',lines,value=TRUE)
  assert(length(header)==1L, 'VCF #CHROM header missing or repeated')
  col <- strsplit(sub('^#','',header), '\t',fixed=FALSE)[[1]]
  records <- lines[!startsWith(lines,'#')]
  assert(length(records)>0, paste('VCF has no variant records:',path))
  dat <- read.delim(text=records, header=FALSE, sep='\t', quote='',
                    comment.char='',check.names=FALSE,stringsAsFactors=FALSE)
  assert(ncol(dat)==length(col), 'VCF record/header column mismatch')
  names(dat) <- col
  sid <- names(dat)[-(1:9)]
  assert(length(sid)==length(unique(sid)), 'Duplicate VCF sample IDs')
  assert(identical(unique(dat$CHROM),REGION_CHROM),'Unexpected VCF chromosome')
  assert(all(dat$POS>=REGION_START & dat$POS<=REGION_END),'VCF positions out of PI4K bounds')
  assert(!is.unsorted(as.integer(dat$POS),strictly=TRUE),'VCF positions not strictly increasing')
  n <- nrow(dat); ns <- length(sid)
  gt <- matrix(NA_integer_,n,ns,dimnames=list(paste0(dat$CHROM,':',dat$POS),sid))
  ambiguous <- matrix(FALSE,n,ns,dimnames=dimnames(gt))
  partial <- matrix(FALSE,n,ns,dimnames=dimnames(gt))
  for (j in seq_len(ns)) {
    fields <- sub(':.*$', '', dat[[j+9L]])
    spl <- strsplit(fields,'[/|]')
    for (i in seq_len(n)) {
      alleles <- spl[[i]]
      if (length(alleles)==1L && grepl('^[0-9]+$',alleles[1])) {
        gt[i,j] <- as.integer(alleles[1]); next
      }
      if (length(alleles)!=2L) {partial[i,j] <- TRUE; next}
      v <- suppressWarnings(as.integer(alleles))
      if (all(!is.na(v)) && v[1]==v[2]) gt[i,j] <- v[1]
      else if (all(!is.na(v))) ambiguous[i,j] <- TRUE
      else if (sum(is.na(v))==1L) partial[i,j] <- TRUE
    }
  }
  meta <- data.frame(chrom=as.character(dat$CHROM), pos=as.integer(dat$POS),
     ref=as.character(dat$REF),alt=as.character(dat$ALT),
     filter=as.character(dat$FILTER), stringsAsFactors=FALSE)
  # REF plus all ALTs, so a 0 or 1 allele code always points to a specific base.
  allele_n <- 1L+lengths(strsplit(meta$alt,',',fixed=TRUE))
  allele_n[meta$alt=='.'] <- 1L
  bad <- which(vapply(seq_len(n),function(i) any(gt[i,,drop=TRUE]>=allele_n[i],na.rm=TRUE),logical(1)))
  assert(!length(bad),paste('Genotype allele code exceeds ALT enumeration:', paste(head(bad),collapse=',')))
  assert(all(meta$filter=='PASS'),'This analysis requires Pf8 PASS-only VCF')
  list(G=gt,ambiguous=ambiguous,partial=partial,sites=meta,ids=sid,path=path)
}
state_from_lists <- function(ids,cohort=c('all','2020')) {
  cohort <- match.arg(cohort)
  mapping <- setNames(rep(NA_character_,length(ids)),ids)
  for (s in STATES) {
    p <- if(cohort=='all') file.path(SAMPLE_LIST_DIR,paste0(s,'_sample_IDs.txt'))
         else file.path(COHORT_VCF_DIR,paste0('PI4K_2020_',s,'_PASS_sites.vcf.gz'))
    state_ids <- if(cohort=='all') trimws(readLines(required(p),warn=FALSE)) else {
      x<-readLines(gzfile(required(p),'rt'),warn=FALSE)
      h<-grep('^#CHROM\t',x,value=TRUE)
      assert(length(h)==1L,paste('Missing VCF header for state',s))
      strsplit(sub('^#','',h),'\t')[[1]][-(1:9)]
    }
    assert(all(state_ids %in% ids),paste('Unexpected',s,'sample ID in',cohort,'mapping'))
    assert(all(is.na(mapping[state_ids])),paste('State ID overlap:',s))
    mapping[state_ids] <- s
  }
  assert(!anyNA(mapping),'Unmapped sample IDs')
  unname(mapping[ids])
}
extract_col <- function(df,possible,label,required_col=TRUE) {
  hit<-intersect(possible,names(df))
  if (!length(hit)) {if(required_col)stop(paste('Missing column',label,'options:',paste(possible,collapse=',')));return(rep(NA,nrow(df)))}
  df[[hit[1]]]
}
load_sample_metadata <- function(ids,state) {
  m<-read_tsv(META)
  names(m)<-sub('^\\ufeff','',names(m))
  sample<-as.character(extract_col(m,c('Sample','sample_id'),'Sample'))
  ii<-match(ids,sample)
  assert(!anyNA(ii),'Missing IDs in source QC metadata')
  recorded_state <- as.character(extract_col(m,c('Admin level 1'),'state')[ii])
  assert(identical(recorded_state,as.character(state)), 'Sample list state disagrees with original Pf8 metadata')
  result<-data.frame(sample_id=ids,state=state,
    study=as.character(extract_col(m,c('Study'),'Study')[ii]),
    year=suppressWarnings(as.integer(as.numeric(extract_col(m,c('Year'),'Year')[ii]))),
    case_id=as.character(extract_col(m,c('All samples same case'),'case')[ii]),
    genome_callable=suppressWarnings(as.numeric(extract_col(m,c('% callable'),'callability')[ii])),
    stringsAsFactors=FALSE)
  f<-read_tsv(FWS)
  result$Fws<-suppressWarnings(as.numeric(extract_col(f,c('Fws'),'Fws')[match(ids,as.character(extract_col(f,c('Sample'),'Sample')))]))
  assert(!anyNA(result$Fws),'Missing Fws values')
  result
}
wilson_ci <- function(k,n,z=1.959963984540054) {
  if (n==0) return(c(NA_real_,NA_real_))
  p<-k/n; d<-1+z*z/n
  center<-(p+z*z/(2*n))/d
  half<-z*sqrt(p*(1-p)/n+z*z/(4*n*n))/d
  pmax(0,c(center-half,center+half))
}
site_stats <- function(G, state, cohort='primary2020') {
  nsite<-nrow(G); pnames<-rownames(G)
  cr<-rowMeans(!is.na(G))
  states_call<-sapply(STATES, function(s) rowMeans(!is.na(G[,state==s,drop=FALSE])))
  allele_count<-vapply(seq_len(nsite), function(i) length(unique(na.omit(G[i,]))),integer(1))
  minor_count<-vapply(seq_len(nsite), function(i) {
    f<-table(G[i,]); if (length(f)!=2L) return(NA_integer_)
    as.integer(min(f))
  },integer(1))
  keep<-cr>=SETTINGS$call_rate & rowSums(states_call>=SETTINGS$state_call_rate)==length(STATES) &
        allele_count==2L & !is.na(minor_count) & minor_count>=SETTINGS$minor_allele_count
  data.frame(id=pnames, callable=cr, states_call, allele_count=allele_count,
             minor_count=minor_count, informative=keep, cohort=cohort,check.names=FALSE)
}
assay_allele_count <- function(G, state, positions, cohort) {
  # Frequencies are per independent sample consensus, never double the encoded slots.
  out<-vector('list',length(STATES))
  for (ss in seq_along(STATES)) {
    s<-STATES[ss]; g<-G[,state==s,drop=FALSE]
    n_obs<-rowSums(!is.na(g)); n_alt<-rowSums(g==1L,na.rm=TRUE)
    out[[ss]]<-data.frame(cohort=cohort,state=s,pos=positions,n_observed=n_obs,
                         n_alt_code1=n_alt,alt_code1_frequency=ifelse(n_obs>0,n_alt/n_obs,NA_real_))
  }
  do.call(rbind,out)
}
export_plot <- function(p,name,w=8,h=5) {
  requireNamespace('ggplot2',quietly=TRUE) || stop('Install ggplot2: Rscript -e "install.packages(\"ggplot2\")"')
  ggplot2::ggsave(file.path(OUT,'figures',paste0(name,'.pdf')), plot=p, width=w,height=h,units='in',device='pdf',useDingbats=FALSE)
  ggplot2::ggsave(file.path(OUT,'figures',paste0(name,'.tiff')), plot=p,width=w,height=h,units='in',dpi=600,compression='lzw')
}
clean_theme <- function() {
  ggplot2::theme_minimal(base_size=11) + ggplot2::theme(panel.grid.minor=ggplot2::element_blank(),
    plot.title=ggplot2::element_text(face='bold',size=12),
    axis.title=ggplot2::element_text(face='bold'),legend.position='bottom')
}
