#' Write alignments from each saturation test random sampling to disk
#' @param data a \dQuote{oligodiag} object 
#' @param dir path where the alignment should be saved
#' @param no3rd a logical wether to print the alignment with all codon position or the 
#' @examples write_saturation_alignments(data)
#' @export
write_saturation_alignments<-function(data, dir = NULL, no3rd = FALSE, ...){
  if (!(is(data) == "oligodiag")){
    stop("The object has incorrect class. Please run plot_saturation or plot_saturation_n to get the correct object class", call. = FALSE)
  }
  
  #   if (no3rd & (length(data$seed) > 1)) {
  #     stop("Exporting the alignment without the 3rd position is only possible with the results from plot_saturation", call. = FALSE)
  #   }
  if (is.null(dir)){
    dir <- getwd()
  }
  dir <- sub('/$', '', dir) 
  date_analysis <-  format(Sys.time(), "%m%d%Y%H%M")
  if (length(data$seed) > 1){
    if (no3rd){
      file_name <- paste(dir,"/saturation_test_alignment_no3rd_", date_analysis, ".fasta",sep = ""); 
      write.dna(data$aln_no_3rd, file = file_name, format = "fasta")     
      cat("Alignmet(s) written to ", file_name, "\n")
    }else{
      write_alns <- lapply(1:length(data$seed), function(x){
        file_name <- paste(dir,"/saturation_test_alignment_", x, "_seed", data$seed[[x]],"_", date_analysis, ".fasta",sep = ""); 
        write.dna(data$aln[[x]], file = file_name, format = "fasta")
      })
    }
  }else{
    if (no3rd){
      if (is.null(data$aln_no_3rd)) {
        stop("Most probably your alignment is not saturated. If you still want to get an alignment without the third codon, please use remove_3rd_codon function",call. = FALSE)
      }
      file_name <- paste(dir,"/saturation_test_alignment_no3rd_", date_analysis, ".fasta",sep = ""); 
      write.dna(data$aln_no_3rd, file = file_name, format = "fasta")
      cat("Alignmet(s) written to ", file_name, "\n")
    }else{
      file_name <- paste(dir,"/saturation_test_alignment_seed", data$seed,"_", date_analysis, ".fasta",sep = ""); 
      write.dna(data$aln, file = file_name, format = "fasta")
      cat("Alignmet(s) written to ", file_name, "\n")
    }
  }
  if (!(no3rd) & (length(data$seed) > 1)) {
    cat("Alignmet(s) written to ", dir, "\n")
  }
}