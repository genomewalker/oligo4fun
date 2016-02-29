#' Plot the distribution of entropy values
#' @param data a dataframe containing the results of running \href{https://github.com/meren/oligotyping}{entropy-analysis} 
#' command from oligotyping package. It could be a raw file obtained after the entropy-analysis or a processed raw file by
#' \code{\link{process_entropy_file}}
#' @param process a logical indicating whether to process the entropy data; defaults to TRUE so 
#' @return a ggplot2 object with  
#'
#' @references Eren, A. M., Maignien, L., Sul, W. J., Murphy, L. G., Grim, S. L., Morrison, H. G., 
#' and Sogin, M. L. (2013). Oligotyping: Differentiating between closely related microbial taxa 
#' using 16S rRNA gene data. Methods Ecol Evol. 4(12)
#'
#' @references Eren, A. M., Morrison, H. G., Lescault, P. J., Reveillaud, J., Vineis, J. H., & Sogin, M. L. (2015). 
#' Minimum entropy decomposition: Unsupervised oligotyping for sensitive partitioning of high-throughput marker gene 
#' sequences. The ISME Journal, 9(4), 968-979
#' 
#' @examples entropy_plot <- plot_entropy_density(entropy)
#' @export
plot_entropy_density <- function(data, process = TRUE){
  process <- .check_entropy_df(data)
  if (process){
    data <- process_entropy_file(data)
    cat("Note: Entropy dataframe has been processed\n")
  }
  colors<-c("#97d400","#00dadb","#db00ec")
  p<-ggplot2::ggplot(data=data, aes(x=entropy)) +
    geom_density(aes(group=codon, fill=codon), alpha=.3) +
    theme_bw() +
    xlab("Entropy") +
    ylab("Density") +
    scale_fill_manual(guide = guide_legend(title = "Codon"), values = colors)+ 
    theme(legend.position="top",legend.key = element_blank())
}

#' Check the entropy dataframe for consistency
#' @param data a dataframe containing the results of running \href{https://github.com/meren/oligotyping}{entropy-analysis} 
#' command from oligotyping package. It could be a raw file obtained after the entropy-analysis or a processed raw file by
#' \code{\link{process_entropy_file}} 
#' @return a logical indicating whether to process the entropy dataframe
#' 
#' @examples check_entropy_df(entropy)
#' @export
.check_entropy_df <- function(data){
  num_cols <- dim(data)[2]
  col_names <- names(data)
  if ((num_cols < 2) || (num_cols >3) || is.null(num_cols)){
    stop("The entropy dataframe doesn't have the correct format. Check the manual for instructions.", call. = FALSE)
  }
  if (( "position" %in% col_names) & 
        ( "entropy" %in% col_names) & 
        ( "codon" %in% col_names) &
        (num_cols == 3)){
    warning("Seems that the entropy dataframe has been already processed. Skipping process step.")
    process <- FALSE
  }else{
    process <- TRUE
  }
}