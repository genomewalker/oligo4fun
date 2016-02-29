#' Write saturation plots from each saturation test random sampling to disk
#' @param data a \dQuote{oligodiag} object 
#' @param dir path where the plot should be saved
#' @param format image format for the saved plot. Should be \dQuote{png}, \dQuote{jpg} or \dQuote{pdf}
#' 
#' 
#' @examples write_saturation_plots(data)
#' @export
write_saturation_plots<-function(data, dir = NULL, format = "png", ...){
  if (!(is(data) == "oligodiag")){
    stop("The object has incorrect class. Please run plot_saturation or plot_saturation_n to get the correct object class", call. = FALSE)
  }
  
  if (is.null(dir)){
    dir <- getwd()
  }
  
  path <- sub('/$', '',dir)
  date <-  format(Sys.time(), "%m%d%Y%H%M")
  
  if (length(data$seed) > 1){
    write_plots <- lapply(1:length(data$seed), function(x){
      file_name <- paste(dir, "/saturation_test_alignment_plot_", x, "_seed", data$seed[[x]],"_", date, ".", format, sep = "")    
      .set_plot_format(format, file_name)
      print(data$plot[[x]])
      dev.off()
    })
  }else{
    file_name <- paste(dir, "/saturation_test_alignment_plot_seed", data$seed, "_", date, ".", format, sep = "")    
    .set_plot_format(format, file_name)
    print(data$plot)
    dev.off() 
  }
  cat("Plot(s) written to ", dir, "\n")
}