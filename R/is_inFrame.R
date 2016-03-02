#' @export
is_inFrame <- function(X){
  if ((class(X) != "DNAbin" && class(X) != "data.frame")){
    stop("Wrong class. Should be a DNAbin object for the alignment or a data.frame for the entropy", call. = FALSE)
  }
  
  if (class(X) == "data.frame" & dim(X)[1] %% 3 != 0){
    stop("Looks like your sequences is not in the correct frame. Please predict the ORFs on your sequence and provide the correct nuclotide file", call. = FALSE)
  }
  
  if (class(X) == "DNAbin" & dim(X)[2] %% 3 != 0){
    stop("Looks like your sequences is not in the correct frame. Please predict the ORFs on your sequence and provide the correct nuclotide file", call. = FALSE)
  }
}