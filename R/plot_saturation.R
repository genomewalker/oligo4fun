#' Draw a saturation plot of the oligotyping input file
#' @param aln a matrix containing the DNA sequences; this must be of class "DNAbin" 
#' @param nseqs a value for the number of sequences to pick
#' @param seed a value for the random number generator. This will allow you to repeat the analysis
#' @param model a character string specifying the evolutionary model to be used by \code{\link[ape]{dist.dna}}
#' @param all a logical indicating whether to use all codon positions; defaults to FALSE so only the third codon position is used.
#' @param verbose a logical indicating whether to show in screen the progress; defaults to TRUE.
#' @details Usually oligotyping datasets are very large and perform a full analysis would be computationally really expensive. We
#' recommend to use up to 5000 random sequences in a desktop computer. You can calculate multiple saturation plots and its associated 
#' statistics using \code{\link[oligo4fun]{plot_saturation_n}} to get a better overview of your alignment.
#' 
#' @return An object of class \dQuote{oligodiag} is a list containing at least the following components:  
#' \describe{
#'  \item{plot}{a ggplot object containing the saturation plot}
#'  \item{seed}{the seed used for picking the random sequences}
#'  \item{aln}{a matrix of the random selected sequences stored in binary format}
#'  \item{stats}{mean and standard deviation of transitions and transversions}
#'  \item{all_codons}{logical if all codon positions have been used}
#'  \item{saturation}{whether your alignemnt possible presents saturation}
#' }
#' 
#' 
#' @examples saturation_plot <- plot_saturation(aln, nseqs = 1000, all = FALSE)
#' @export
plot_saturation<-function(aln = aln, nseqs = 1000, model = "K80", all = FALSE, verbose = TRUE, seed = 0, rsamp = FALSE, alns = FALSE, ...){
  aln_num_seqs <- dim(aln)[1]
  if (aln_num_seqs <= 1000) {
    rand_seqs <- aln
    seed <- NULL
    if (verbose){
      cat("Used all sequences\n\n")
    }
  }else{
    rand_seqs <- random_sequences(aln = aln, nseqs = nseqs, seed = seed, verbose = verbose, ...)
    if (verbose){
      cat(paste("Randomly selected", nseqs, "sequences\n\n"))
    }
  }
  
  if (!(all)) {
    rand_seqs$aln <- get_3rd_codon(rand_seqs$aln)
    if (verbose){
      cat("Sequence alignment cointaining only the 3rd codon position\n")
    }
  }
  
  if (verbose){
    cat("Calculating the number of pairwise transitions and transversions\n")
  }
  subs <- spider::titv(rand_seqs$aln)
  
  if (verbose){
    cat(paste("Computing matrix of pairwise distances using model",model,"\n"))
  }
  aln_dist <- ape::dist.dna(rand_seqs$aln, model)
  #Transversions
  tv <- t(subs)
  tv <- tv[lower.tri(tv)]
  #Transitions
  ti <- subs[lower.tri(subs)]
  ti_df <- data.frame(x=as.vector(aln_dist),y=ti, m="s")
  tv_df <- data.frame(x=as.vector(aln_dist),y=tv, m="v")
  ti_tv <- as.data.frame(rbind(ti_df,tv_df))
  
  if (verbose){
    cat("Calculating statistics\n")
  }
  
  ti_tv_stats <- list(
    mean_ti=mean(ti),
    sd_ti=sd(ti),
    mean_tv=mean(tv),
    sd_tv=sd(tv)
  )
  
  if (verbose){
    cat("Drawing saturation plot\n")
  }
  require(ggplot2)
  p <- ggplot2::ggplot(aes(x,y,color=m),data=unique(ti_tv)) +
    geom_hline(yintercept = mean(ti), size=0.5, color="#ff800e",) +
    geom_hline(yintercept = c(mean(ti)+(sd(ti)),mean(ti)-((sd(ti)))) , size=0.5, color="#ff800e",linetype=2) +
    geom_hline(yintercept = mean(tv) ,size=0.5, color="#0b86cf") +
    geom_hline(yintercept = c(mean(tv)+(sd(tv)),mean(tv)-((sd(tv)))) , size=0.5, color="#0b86cf",linetype=2)+
    geom_point(shape=15, alpha=0.8, size=3) + theme_bw() +
    scale_color_manual(values=c("#ff800e","#0b86cf"), labels=c("Transitions", "Transversions"), name="") +
    theme(legend.position="top",legend.key = element_blank(),
          legend.text = element_text( size = 11),
          axis.text=element_text(size=10),axis.title=element_text(size=12)) +
    xlab(paste("Genetic distance (", model, ")", sep = "")) + ylab("Transitions and Transversions") +
    guides(colour=guide_legend(override.aes=list(size=4)))
  
  saturation <- mean(ti) > mean(tv)
  
  if (saturation){
    aln_no_3rd <- remove_3rd_codon(aln)
    if (verbose){
      cat("Results suggests that there is saturation\n")
      cat("Removing 3rd codon position to the alignment\n")
      cat("Note: Substitution saturation occurs when the frequency of transitions (mean:", round(mean(ti), 2),  ") overtakes the frequency of transversions (mean:", round(mean(tv), 2),  ")\n", sep = "")
    }
  }else{
    if (verbose){
      cat("Results suggests that there is no saturation\n")
      cat("Note: Substitution saturation occurs when the frequency of transitions (", round(mean(ti), 2),  ") overtakes the frequency of transversions(", round(mean(ti), 2),  ")\n", sep = "")
    }
    aln_no_3rd <- NULL
  }
  
  if (rsamp){
    results<-list(
      plot = p,
      if (alns){
        aln = rand_seqs$aln,
      }
      seed = rand_seqs$seed,
      stats = ti_tv_stats,
      all_codons = all,
      ti_tv = dplyr::tbl_df(ti_tv),
      saturation = ifelse(saturation, "Results suggests that there is saturation.", "Results suggests that there is no saturation.")    )
  }else{
    results<-list(
      plot = p,
      if (alns){
        aln = rand_seqs$aln,
      }
      seed = rand_seqs$seed,
      stats = ti_tv_stats,
      all_codons = all,
      ti_tv = dplyr::tbl_df(ti_tv),
      saturation = ifelse(saturation, "Results suggests that there is saturation.", "Results suggests that there is no saturation."),
      aln_no_3rd = aln_no_3rd
    )
  }
  if (saturation & !(rsamp)){
    class(results$aln_no_3rd) <- "DNAbin"
  }
  if(alns){
    class(results$aln)<-"DNAbin"
  }
  class(results)<-"oligdiag"
  return(results)
}
