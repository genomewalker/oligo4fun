#' Draw a saturation plot of the oligotyping input file (n) times
#' @description Perform multiple random subsamplings of your alignemnt. 
#' @param aln a matrix containing the DNA sequences; this must be of class "DNAbin" 
#' @param rep number of random subsamplings
#' @param nseqs a value for the number of sequences to pick
#' @param model a character string specifying the evolutionary model to be used by \code{\link[ape]{dist.dna}}
#' @param all a logical indicating whether to use all codon positions; defaults to FALSE so only the third codon position is used.
#' @param parallel a logical indicating whether to do the random subsampling on a sequential way or use the \code{\href{https://github.com/tudo-r/BatchJobs}{BatchJobs}}
#' framework to distribute the random subsamplings on different cores or in computer cluster
#' @param reg_id Name of registry. Displayed e.g. in mails or in cluster queue
#' @param reg_dir Path where files regarding the registry / jobs should be saved
#' @param conf_file Location of the configuration file to load
#' @details You can calculate multiple saturation plots and its associated statistics. The different
#' random subsamplings can be easily distributed over multiple cores or in a computer cluster.
#' 
#' @return An object of class \dQuote{oligodiag} is a list containing at least the following components:  
#' \describe{
#'  \item{plot}{a ggplot object containing the saturation plots of each random sampling}
#'  \item{seed}{the seeds used for picking the random sequences}
#'  \item{aln}{a matrix of the random selected sequences stored in binary format for each random subsampling}
#'  \item{combined_stats}{mean and standard deviation of transitions and transversions for each random subsampling}
#'  \item{saturation}{whether your alignment presents saturation for each random subsampling}
#'  \item{raw}{raw results for each random subsampling}
#' }
#' 
#' 
#' @examples saturation_plots <- plot_saturation_n(aln, nseqs = 1000, rep = 100, parallel = F)
#' @examples saturation_plots <- plot_saturation_n(aln, nseqs = 1000, rep = 100, parallel = T, reg_id = "test_id", reg_dir = "test-dir", conf_file = ".BatchJobs.R")
#' @export
plot_saturation_n<-function(aln = aln, rep = 100, nseqs = 1000, model = "K80", parallel = FALSE,
                            all = FALSE, reg_id = NULL, reg_dir = NULL, conf_file = NULL, 
                            job_res = list(), ...){
  results <- list()
  if (parallel){
    if ((is.null(conf_file)) || (is.null(reg_dir)) || (is.null(reg_id))){
      stop("Please add the configuration file and registry values requiered for BatchJobs" , call. = FALSE)
    }
    
    function_args <- list()
    fun <- plot_saturation
    function_args <- list(aln = aln, nseqs = nseqs, model = model, all = all, verbose = FALSE, seed = 0, rsamp = TRUE)
    iterations <- 1:rep
    
    batch_function <- function(X) {
      tmp <- do.call(fun, function_args)
    }
    
    BatchJobs::loadConfig(conf_file) 
    reg <- BatchJobs::makeRegistry(id=reg_id, file.dir=reg_dir)
    id  <- BatchJobs::batchMap(reg, batch_function, iterations)
    plot_submission <- BatchJobs::submitJobs(reg, resources=job_res)
    plot_run <- BatchJobs::waitForJobs(reg, id)
    plot_runs <- reduceResultsList(reg)
    removeRegistry(reg, ask = "no")
    
    if (!plot_run) {
      stop('Error in batch jobs', call. = FALSE)
    }
    
  }else{
    plot_runs <- plyr::llply(1:rep, plot_saturation, ..., aln = aln, nseqs = nseqs, model = model, all = all, verbose = FALSE,
                             seed = 0, rsamp = TRUE, .parallel = F, .progress = plyr::progress_text(width = 80))
  }
  
  results$combined_stats <- dplyr::rbind_all(lapply(plot_runs, function(x) dplyr::rbind_list(x[["stats"]])))
  results$raw <- plot_runs
  results$plot <- lapply(plot_runs, function(x) x[["plot"]])
  results$aln <- lapply(plot_runs, function(x) x[["aln"]])
  results$aln <- lapply(results$aln, function(x) {class(x)<-"DNAbin"; x})
  results$seed <- lapply(plot_runs, function(x) x[["seed"]])
  results$saturation <- lapply(plot_runs, function(x) x[["saturation"]])
  results$aln_no_3rd <- remove_3rd_codon(aln)
  class(results)<-"oligodiag"
  return(results)
}