#' Create a sequence alignment with only the 3rd codon
#' @param aln a sequence alignment 
#' @export
get_3rd_codon <- function(aln, ...){
  is_inFrame(aln)
  aln_3rd<-aln[, seq(0, ncol(aln), by=3)]
}