#' Write a density plot of the transitions and transversions mean values from each saturation test random sampling to disk
#' @param data a \dQuote{oligodiag} object 
#' @param dir path where the plot should be saved
#' @param format image format for the saved plot. Should be \dQuote{png}, \dQuote{jpg} or \dQuote{pdf}
#' @param show a logical whether to show the plot
#' 
#' @examples write_saturation_plot_density(data)
#' @export
write_saturation_plot_density<-function(data, dir = NULL, format = "png", show = TRUE, ...){
  if (!(is(data) == "oligodiag")){
    stop("The object has incorrect class. Please run estimation_saturation or estimation_saturation_n to get the correct object class", call. = FALSE)
  }
  
  if (is.null(dir)){
    dir <- getwd()
  }
  
  dir <- sub('/$', '', dir)
  date <-  format(Sys.time(), "%m%d%Y%H%M")
  
  file_name <- paste(dir, "/saturation_test_alignment_density_", date, ".", format, sep = "")    
  
  if (dim(data$combined_stats)[1] < 3){
    stop("At least 3 random searches needed to plot. Increase 'rep' to a value > 2", call. = FALSE)
  }
  
  combined_stats <- data$combined_stats  %>% 
    dplyr::mutate(run=row.names(data$combined_stats)) %>% 
    dplyr::select(run, mean_ti, mean_tv) %>% 
    tidyr::gather(mean, value, c(mean_ti, mean_tv))
  
  .set_plot_format(format, file_name)
  
  p<-ggplot2::ggplot(data=combined_stats, aes(x=value)) +
    geom_density(aes(group=mean, fill=mean), alpha=.3) +
    scale_fill_manual(values=c("#ff800e","#0b86cf"), labels=c("Transitions", "Transversions"), name="") +
    theme_bw() +
    xlab("Average value") +
    ylab("Density") +
    theme(legend.position="top",legend.key = element_blank())
  print(p)
  dev.off()
  if (show){
  print(p)
  }
  cat("Plots written to ", file_name, "\n")
}