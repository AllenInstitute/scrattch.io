read_seurat_dgCMatrix <- function(seurat_obj,
                                  which_data = c("raw.data","data","scale.data")) {

  if(which_data == "raw.data") {
    if(is.null(suerat_obj@raw_data)) {
      stop("This seurat object doesn't have anything in the raw.data slot.")
    }

    return(seurat_obj@raw_data)

  } else if (which_data == "data") {
    if(is.null(suerat_obj@data)) {
      stop("This seurat object doesn't have anything in the data slot.")
    }

    return(seurat_obj@data)

  } else if (which_data == "scale.data") {
    if(is.null(suerat_obj@raw_data)) {
      stop("This seurat object doesn't have anything in the scale.data slot.")
    }

    return(seurat_obj@scale.data)

  }

}

read_seurat_anno <- function(seurat_obj) {

  if(is.null(suerat_obj@meta)) {
    stop("This seurat object doesn't have anything in the meta slot.")
  }

  meta <- seurat_obj@meta
  meta <- tibble::rownames_to_column(meta, "sample_name")
  meta <- dplyr::select(meta, c("sample_name", dplyr::everything()))

  return(meta)
}
