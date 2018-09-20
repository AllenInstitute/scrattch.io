#' Read a whole sparse matrix directly from a .h5ad file
#'
#' @param h5ad .h5ad file to read.
#' @param target The sparse matrix to read within the .h5ad file. Default = "/raw.X"
#'
read_h5ad_dgCMatrix <- function(h5ad,
                               target = "/raw.X") {
  library(rhdf5)
  library(purrr)
  library(dplyr)
  library(Matrix)

  H5close()

  root <- H5Fopen(h5ad)

  i_path <- paste0(target,"/indices")
  p_path <- paste0(target,"/indptr")
  x_path <- paste0(target,"/data")

  print("Reading indices")
  i <- as.vector(h5read(root, i_path))
  print("Reading pointers")
  p <- h5read(root, p_path)
  print("Reading values")
  x <- as.vector(h5read(root, x_path))
  print("Reading observations")
  o <- as.vector(h5read(root, "/obs")$index)
  print("Reading variables")
  v <- as.vector(h5read(root, "/var")$index)

  print("Reading dimensions")
  dims <- c(length(v), length(o))

  H5Fclose(root)
  H5close()

  print("Assembling dgCMatrix")
  m <- sparseMatrix(i = i,
               p = p,
               x = x,
               index1 = FALSE,
               dims = dims)

  rownames(m) <- v
  colnames(m) <- o

  return(m)

}
