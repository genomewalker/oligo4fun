#' Create a sequence alignment with only the 1st and 2nd codon
#' @param aln a sequence alignment 
#' @export
remove_3rd_codon <- function(aln, ...){
  aln_1st_2nd<-aln[, -seq(0, ncol(aln), by=3)]
}