#' Write a dgCMatrix to an h5 file similar to cellRanger format
#'
#' @param mat a dgCMatrix to write.
#' @param cols_are Whether columns are "gene_names" or "sample_names". If "gene_names", mat will be transposed to match 10X conventions.
#' @param h5_target The target .h5 file for output.
#' @param ref_name Reference name for storing the data
#' @param gene_ids If available, ENSEMBL IDs for each gene
#'
#' @return an .h5 file with these mappings from dgCMatrix -> .h5:
#' \itemize{
#'   \item colnames(mat) -> /ref_name/barcodes (after transposition if cols_are == "gene_names")
#'   \item rownames(mat) -> /ref_name/gene_names (after transposition if cols_are == "gene_names")
#'   \item mat@x -> /ref_name/data
#'   \item mat@i -> /ref_name/indices
#'   \item mat@p -> /ref_name/indptr
#'   \item gene_ids -> /ref_name/gene
#' }
#'
write_dgCMatrix_h5 <- function(mat,
                               cols_are = "gene_names",
                               h5_target,
                               ref_name = "mm10-1.2.0_premrna",
                               gene_ids = NULL) {

  #library(Matrix)

  if(grepl("gene",cols_are)) {
    mat <- Matrix::t(mat)
  }

  # Create target file
  rhdf5::h5createFile(h5_target)
  # Create data group
  rhdf5::h5createGroup(h5_target,
                       ref_name)

  # Store sample ids (barcodes) and gene names
  rhdf5::h5write(colnames(mat),
                 h5_target,
                 paste0("/",ref_name,"/barcodes"))
  rhdf5::h5write(rownames(mat),
                 h5_target,
                 paste0("/",ref_name,"/gene_names"))

  if(is.null(gene_ids)) {
    gene_ids <- rownames(mat)
  }

  rhdf5::h5write(gene_ids,
                 h5_target,
                 paste0("/",ref_name,"/gene"))

  # Store dimensions as shape
  rhdf5::h5write(dim(mat),
                 h5_target,
                 paste0("/",ref_name,"/shape"))

  # Store values from mat@x as data
  rhdf5::h5createDataset(h5_target,
                         paste0("/",ref_name,"/data"),
                         dims = length(mat@x),
                         storage.mode = "integer",
                         chunk = 1000,
                         level = 4)
  rhdf5::h5write(mat@x,
                 h5_target,
                 paste0("/",ref_name,"/data"))

  # Store row indices from mat@i as indices
  rhdf5::h5createDataset(h5_target,
                         paste0("/",ref_name,"/indices"),
                         dims = length(mat@i),
                         storage.mode = "integer",
                         chunk = 1000,
                         level = 4)
  rhdf5::h5write(mat@i,
                 h5_target,
                 paste0("/",ref_name,"/indices"))

  # Store column pointers from mat@p as indptr
  rhdf5::h5write(mat@p,
                 h5_target,
                 paste0("/",ref_name,"/indptr"))

}
