#' Read a whole sparse matrix directly from a 10X-style .h5 file
#'
#' @param h5 .h5 file to read.
#' @param target The sparse matrix to read within the .h5 file.
#'
read_10x_dgCMatrix <- function(h5,
                               target) {
  #library(Matrix)

  root <- rhdf5::H5Fopen(h5)

  i_path <- paste0(target,"/indices")
  p_path <- paste0(target,"/indptr")
  x_path <- paste0(target,"/data")
  dims_path <- paste0(target,"/shape")

  print("Reading indices")
  i <- read_tome_vector(root, i_path)
  print("Reading pointers")
  p <- read_tome_vector(root, p_path)
  print("Reading values")
  x <- read_tome_vector(root, x_path)
  print("Reading dimensions")
  dims <- read_tome_vector(root, dims_path)

  H5Fclose(root)

  print("Assembling dgCMatrix")
  Matrix::sparseMatrix(i = i,
                       p = p,
                       x = x,
                       index1 = FALSE,
                       dims = dims)

}
