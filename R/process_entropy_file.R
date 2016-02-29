#' Pre-process entropy-analysis file
#' @description Pre-process the entropy-analysis file containing the results after running 
#' \href{https://github.com/meren/oligotyping}{entropy-analysis} command from oligotyping package
#' @param file a dataframe containing the results after running 
#' \href{https://github.com/meren/oligotyping}{entropy-analysis} command from oligotyping package 
#' 
#' @details Transforms the generated entropy file for plotting. 
#' 
#' @return A dataframe with the nucleotide positions, the entropy values and the codon position  
#'
#' @references Eren, A. M., Morrison, H. G., Lescault, P. J., Reveillaud, J., Vineis, J. H., & Sogin, M. L. (2015). 
#' Minimum entropy decomposition: Unsupervised oligotyping for sensitive partitioning of high-throughput marker gene 
#' sequences. The ISME Journal, 9(4), 968-979
#' 
#' @examples entropy_file <- process_entropy_file(entropy)
#' @export
process_entropy_file <- function(file = file){
  file <- file %>% 
    dplyr::rename(position = V1, entropy = V2) %>% 
    dplyr::arrange(position) %>% 
    dplyr::mutate(position = position + 1, codon = rep(c("1st", "2nd", "3rd"), dim(file)[1]/3)) 
}

