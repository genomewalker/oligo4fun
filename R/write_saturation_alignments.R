#' Write alignments from each saturation test random sampling to disk
#' @param data a \dQuote{oligodiag} object 
#' @param dir path where the alignment should be saved
#' @param no3rd a logical wether to print the original alignment with all codon position or without the 3rd position
#' @examples write_saturation_alignments(data)
#' @export
write_saturation_alignments<-function(data, dir = NULL, no3rd = FALSE, ...){
  if (!(is(data) == "oligodiag")){
    stop("The object has incorrect class. Please run estimate_saturation or estimate_saturation_n to get the correct object class", call. = FALSE)
  }
  
  if (is.null(dir)){
    dir <- getwd()
  }
  
  dir <- sub('/$', '', dir) 
  date_analysis <-  format(Sys.time(), "%m%d%Y%H%M")
  
  # if it's not 3rd we need to extract the sequences again using the seeds 
  
  if (length(data$seed) > 1){
    if (no3rd){
      if (is.null(data$aln_no_3rd)) {
        stop("Most probably your alignment is not saturated. If you still want to get an alignment without the third codon, please use remove_3rd_codon function",call. = FALSE)
      }
      cat("Saving original alignment without the 3rd codon position\n")
      file_name <- paste(dir,"/saturation_test_alignment_no3rd_", date_analysis, ".fasta",sep = ""); 
      write.dna(data$aln_no_3rd, file = file_name, format = "fasta")     
      cat("Alignmet(s) written to", file_name, "\n")
    }else{
      cat("Saving random alignments used for the saturation estimation\n")
      write_alns <- plyr::llply(1:length(data$seed), function(x){
        rand_aln <- random_sequences(aln = data$aln, seed = data$seed[[x]], nseqs = data$nseqs, verbose = FALSE)
        file_name <- paste(dir,"/saturation_test_alignment_seed", data$seed[[x]],"_", x, '_', date_analysis, ".fasta",sep = ""); 
        write.dna(rand_aln$aln, file = file_name, format = "fasta")
      }, .parallel = F, .progress = plyr::progress_text(width = 80))
    }
  }else{
    if (no3rd){
      if (is.null(data$aln_no_3rd)) {
        stop("Most probably your alignment is not saturated. If you still want to get an alignment without the third codon, please use remove_3rd_codon function",call. = FALSE)
      }
      cat("Saving original alignment without the 3rd codon position\n")
      file_name <- paste(dir,"/saturation_test_alignment_no3rd_", date_analysis, ".fasta",sep = ""); 
      write.dna(data$aln_no_3rd, file = file_name, format = "fasta")
      cat("Alignmet(s) written to", file_name, "\n")
    }else{
      cat("Saving random alignment used for the saturation estimation\n")
      rand_aln <- random_sequences(aln = data$aln, seed = data$seed)
      file_name <- paste(dir,"/saturation_test_alignment_seed", data$seed,"_", date_analysis, ".fasta",sep = ""); 
      write.dna(rand_aln$aln, file = file_name, format = "fasta")
      cat("Alignmet written to", file_name, "\n")
    }
  }
  if (!(no3rd) & (length(data$seed) > 1)) {
    cat("Alignmet(s) written to ", dir, "\n")
  }
}