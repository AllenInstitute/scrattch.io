#' Write scrattch data to feather files
build_feather <- function(anno = NULL, data = NULL, desc = NULL, feather_dir = NULL) {

  library(dplyr)
  library(feather)

  if(is.null(anno)) {
    stop("anno table required.")
  } else if(is.null(data)) {
    stop("data table required.")
  } else if(is.null(desc)) {
    stop("desc table required.")
  } else {

    # Check for old data format
    if(names(data)[1] == "gene") {

      # if old format is used, transpose the table
      gn <- data$gene
      data <- t(data[,-1])
      colnames(data) <- gn
      data <- as.data.frame(data,stringsAsFactors=F)
      data <- cbind(data,sample_id = rownames(data),stringsAsFactors=F)

    }

    if(!dir.exists(feather_dir)) {
      dir.create(feather_dir)
    }

    datafile <- paste0(feather_dir,"/data.feather")
    cat("Writing data table to",datafile,"\n")
    feather::write_feather(data, datafile)

    annofile <- paste0(feather_dir,"/anno.feather")
    cat("Writing anno table to",annofile,"\n")
    feather::write_feather(anno, annofile)

    descfile <- paste0(feather_dir,"/desc.feather")
    cat("Writing desc table to",descfile,"\n")
    feather::write_feather(desc, descfile)

  }
}
