
#' \emph{Rhodopirellula} sp. alignemnts used for Zure et al., 2016
#'
#' A dataset containing the aligned reads for the study described in 
#'  Zure et al., 2016. The variables are as follows:
#'
#' @references Zure, M., Fernandez-Guerra, A., Munn C. and Harder, J. (2016). High resolution distribution of 
#'  closely related \emph{Rhodopirellula} species in European coastal sediments. Submitted.
#'  
#' @format A DNAbin object with 361890 sequences of 483 nucleotides each
"aln"

#' Entropy values for each position of the alignment used for Zure et al., 2016
#'
#' A dataset containing the entropy values for each 483 nucleotides contained in the 
#'  alignemnt. This is a raw dataframe containing the results of running \href{https://github.com/meren/oligotyping}{entropy-analysis} 
#' command from oligotyping package.
#'
#' @references Eren, A. M., Morrison, H. G., Lescault, P. J., Reveillaud, J., Vineis, J. H., & Sogin, M. L. (2015). 
#' Minimum entropy decomposition: Unsupervised oligotyping for sensitive partitioning of high-throughput marker gene 
#' sequences. The ISME Journal, 9(4), 968-979
#' 
#' @references Zure, M., Fernandez-Guerra, A., Munn C. and Harder, J. (2016). High resolution distribution of 
#'  closely related \emph{Rhodopirellula} species in European coastal sediments. Submitted.
#'
#' 
#' @format A data frame with 483 rows and 2 variables:
#' \itemize{
#'   \item V1: nucleotide positions of the alignment
#'   \item V2: entropy value
#' }
"entropy"