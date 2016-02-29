#' @export
.set_plot_format <- function(format, file_name){
  if((format == "pdf")){
    pdf(file_name, height = 8.27, width = 11.69)
  }else if ((format == "png")){
    png(file_name, height = 8.27, width = 11.69, units = "in", res = 300)
  }else if ((format == "jpg")){
    jpeg(file_name, height = 8.27, width = 11.69, units = "in", res = 300)
  }else{
    stop("Format not valid. Please choose 'png', 'pdf' or 'jpg'\n", call. = FALSE)
  }
}