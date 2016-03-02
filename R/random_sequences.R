#' Pick random sequences from the alignment
#' @param aln a sequence alignment 
#' @param nseqs number of sequences to pick
#' @param seed choose a seed for the random number generator. This will allow you to repeat the analysis.
#' @return An object of class \dQuote{oligodiag} is a list containing at least the following components:  
#' \describe{
#'  \item{seed}{the seed used for picking the random sequences}
#'  \item{aln}{a matrix of the random selected sequences stored in binary format}
#' }
#' 
#' @examples rand_seqs <- random_sequences(aln)
#' @export
random_sequences <-  function(aln = aln, nseqs = 1000, seed = 0, verbose = FALSE, ...){
  is_inFrame(aln)
  aln_num_seqs <- dim(aln)[1]
  
  if (nseqs > aln_num_seqs){
    stop("Reduce the number of random sequences to pick. Max: ", aln_num_seqs, call.=FALSE)
  }
  
  if (seed == 0){
    seed<-sample(1000000:9999999,1)
    if (verbose){
    cat("Note: A random seed has been selected\n")
    }
  }
  
  if (verbose){
  cat(paste("Seed used:"),seed,"\n")
  }
  set.seed(seed)
  rand_seqs_ids <- sample(1:aln_num_seqs, nseqs)
  rand_seqs <- aln[rand_seqs_ids,]
  #cat("Number of random sequences extracted:", dim(rand_seqs)[1], "\n")
  results <- list(
    seed = seed,
    aln = rand_seqs
  )
  class(results)<-"oligodiag"
  return(results)
}