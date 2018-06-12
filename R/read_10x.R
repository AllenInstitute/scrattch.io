#' Read a whole sparse matrix directly from a 10X-style .h5 file
#'
#' @param h5 .h5 file to read.
#' @param target The sparse matrix to read within the .h5 file.
#'
read_10x_dgCMatrix <- function(h5,
                               target) {
  library(rhdf5)
  library(purrr)
  library(dplyr)
  library(Matrix)

  H5close()

  root <- H5Fopen(h5)

  i_path <- paste0(target,"/indices")
  p_path <- paste0(target,"/indptr")
  x_path <- paste0(target,"/data")
  dims_path <- paste0(target,"/shape")

  print("Reading indices")
  i <- as.vector(h5read(root, i_path))
  print("Reading pointers")
  p <- h5read(root, p_path)
  print("Reading values")
  x <- as.vector(h5read(root, x_path))
  print("Reading dimensions")
  dims <- h5read(root, dims_path)

  H5Fclose(root)
  H5close()

  print("Assembling dgCMatrix")
  sparseMatrix(i = i,
               p = p,
               x = x,
               index1 = FALSE,
               dims = dims)

}
