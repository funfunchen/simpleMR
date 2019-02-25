#' Standardize the input for two sampleMR analysis
#'
#' 
#'
#' @param rsid the SNPs ID
#' @param expob the BETAs of sample one (exposure)
#' @param expose the SEs of sample one (exposure)
#' @param outb the BETAs of sample two (outcome)
#' @param outse the SEs of sample one (outcome)
#' @param expo.id the id labeled to sample one
#' @param out.id the id labled to sample two
#'
#' @return A data frame ready for the mr function
#'
#' @export
format_mr <- function(rsid, expob, expose, outb, outse, expo.id="expo", out.id="outcome"){
  if(missing(rsid)) {
    rsid<-seq(from=1,to=length(expob),by=1)
    warning("Missing SNP IDs; Generating placeholders")
  } 
  
  Mr.Data<-tibble(rsid, expob, expose, outb, outse) %>% filter(expose > 0 & outse > 0) %>% 
    mutate(expo.id=expo.id, out.id=out.id) %>% as.data.frame()
  
  #Rename variables based on MRBase conventional names
  names(Mr.Data)<-c("SNP", "beta.exposure", "se.exposure", "beta.outcome", "se.outcome", "expo.id", "out.id")
  
  return(Mr.Data)
}


#' Screen the top hit SNPs
#'
#' Get the top hit SNPs in given \code{window.size} regions.
#'
#' @param infile the GWAS/meta-analysis file contains association summary statistic
#' @param window.size the size of region for each Top SNP, could be \code{10^6} or \code{2*10^6}
#' @param pvalue the pvaule you want to use as threshold, use the quotation mark("")
#'
#' @return A data frame with top hit SNPs
#'
#' @export
screen_top_snp <- function(infile, window.size, pvalue){
  if (nrow(infile) == 0) stop("no records in the file");
  new.data <- infile[0,];
  for (i in unique(infile$CHROM)){
    file.chrom <- filter(infile, CHROM == i) %>% arrange_(pvalue); #  use the standard evaluation versions of the dplyr functions (just append '_' to the function names, ie. group_by_ & summarise_)
    new.data <- bind_rows(new.data, file.chrom[1,], .id = NULL) %>% filter(!is.na(CHROM));  ##find the smallest p value, pass to new.data
    sign.exist=TRUE;  ## get rid of the snps within the "window" of the snp we choose
    while(sign.exist){
      min.pos <- file.chrom[1,]$POS;
      file.chrom <- filter(file.chrom, !(POS>(min.pos-window.size/2) & POS<(min.pos+window.size/2))) %>% arrange_(pvalue);
      new.data <- bind_rows(new.data, file.chrom[1,], .id = NULL);
      if(length(file.chrom$POS)==0) {
        sign.exist=FALSE;
      }
    }
  }
  new.data <- filter(new.data, !is.na(CHROM)) %>% arrange(CHROM,POS);
}
