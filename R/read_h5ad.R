#' Read a whole sparse matrix directly from a .h5ad file
#'
#' @param h5ad .h5ad file to read.
#' @param target The sparse matrix to read within the .h5ad file. Default = "/raw.X"
#'
read_h5ad_dgCMatrix <- function(h5ad,
                                target = "/raw.X") {
  library(rhdf5)
  library(Matrix)

  root <- rhdf5::H5Fopen(h5ad)

  i_path <- paste0(target,"/indices")
  p_path <- paste0(target,"/indptr")
  x_path <- paste0(target,"/data")

  print("Reading indices")
  i <- read_tome_vector(root, i_path)
  print("Reading pointers")
  p <- read_tome_vector(root, p_path)
  print("Reading values")
  x <- read_tome_vector(root, x_path)
  print("Reading observations")
  o <- as.vector(rhdf5::h5read(root, "/obs")$index)
  print("Reading variables")
  v <- as.vector(rhdf5::h5read(root, "/var")$index)

  print("Reading dimensions")
  dims <- c(length(v), length(o))

  H5Fclose(root)

  print("Assembling dgCMatrix")
  m <- Matrix::sparseMatrix(i = i,
                            p = p,
                            x = x,
                            index1 = FALSE,
                            dims = dims)

  rownames(m) <- v
  colnames(m) <- o

  return(m)

}
