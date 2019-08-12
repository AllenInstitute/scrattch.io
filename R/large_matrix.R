#' Convert a large matrix object to a dgCMatrix
#'
#' May be useful in cases where calling Matrix() directly gives a "negative length vectors are not allowed" error.
#'
#' This works through the matrix in column-based chunks becaues this is much faster for conversion to column-indexed dgCMatrices.
#'
#' @param mat The matrix object to convert
#' @param chunk_size The number of columns to read as a chunk. For ~30k genes, a chunk of 5000 (the default) is ~1GB in memory.
#'
#' @return A dgCMatrix object in the same orientation as the matrix object.
#'
large_matrix_to_dgCMatrix <- function(mat,
                                      chunk_size = 5000) {
  #library(Matrix)

  # Row names
  row_names <- rownames(mat)

  # Sample names are in col_attrs/CellID
  col_names <- colnames(mat)

  n_row <- length(row_names)
  n_col <- length(col_names)

  n_chunks <- floor(n_col/chunk_size)

  # Initial chunk
  chunk <- 1
  chunk_start <- 1 + chunk_size * (chunk - 1)
  chunk_end <- chunk_size * chunk
  print(paste0("Adding rows ",chunk_start," to ", chunk_end))
  chunk_mat <- mat[,chunk_start:chunk_end]
  all_sparse <- Matrix::Matrix(chunk_mat, sparse = T)

  # chunkation over chunks
  for(chunk in 2:n_chunks) {
    chunk_start <- 1 + chunk_size * (chunk - 1)
    chunk_end <- chunk_size * chunk
    print(paste0("Adding rows ",chunk_start," to ", chunk_end))
    chunk_mat <- mat[,chunk_start:chunk_end]
    chunk_sparse <- Matrix::Matrix(chunk_mat, sparse = T)
    all_sparse <- cbind(all_sparse, chunk_mat)
  }

  # Remaining samples
  chunk_mat <- mat[,(n_chunks*chunk_size + 1):n_col]
  chunk_sparse <- Matrix::Matrix(chunk_mat, sparse = T)
  all_sparse <- cbind(all_sparse, chunk_mat)

  rownames(all_sparse) <- row_names
  colnames(all_sparse) <- col_names

  return(all_sparse)

}

#' Convert a large dgCMatrix object to a matrix
#'
#' May be useful in cases where calling as.matrix() directly gives a "Cholmod error 'problem too large'" error.
#'
#' This works through the matrix in column-based chunks becaues this is much faster for conversion from column-indexed dgCMatrices.
#'
#' It should still work on other types of Matrix-package classes, but may be slower.
#'
#' @param mat The dgCMatrix object to convert
#' @param chunk_size The number of columns to read as a chunk. For ~30k genes, a chunk of 5000 (the default) is ~1GB in memory, and fits within the CHOLMOD limits.
#'
#' @return A standard matrix object in the same orientation as the dgCMatrix object.
#'
large_dgCMatrix_to_matrix <- function (mat,
                                       chunk_size = 5000) {

  row_names <- rownames(mat)
  col_names <- colnames(mat)
  n_row <- nrow(mat)
  n_col <- ncol(mat)

  print("Allocating initial dense matrix")
  all_dense <- matrix(0, ncol = n_col, nrow = n_row)

  n_chunks <- floor(n_col/chunk_size)
  chunk <- 1
  chunk_start <- 1 + chunk_size * (chunk - 1)
  chunk_end <- chunk_size * chunk
  print(paste0("Adding cols ", chunk_start, " to ", chunk_end))

  chunk_mat <- as(mat[, chunk_start:chunk_end], "matrix")
  all_dense[, chunk_start:chunk_end] <- chunk_mat

  for (chunk in 2:n_chunks) {
    chunk_start <- 1 + chunk_size * (chunk - 1)
    chunk_end <- chunk_size * chunk
    print(paste0("Adding cols ", chunk_start, " to ", chunk_end))
    chunk_mat <- as(mat[,chunk_start:chunk_end], "matrix")
    all_dense[,chunk_start:chunk_end] <- chunk_mat
  }
  chunk_start <- n_chunks * chunk_size + 1
  chunk_end <- n_col
  print(paste0("Adding cols ", chunk_start, " to ", chunk_end))
  chunk_mat <- as(mat[, (n_chunks * chunk_size + 1):n_col], "matrix")
  all_dense[,chunk_start:chunk_end] <- chunk_mat

  if(!is.null(row_names)) {
    rownames(all_dense) <- row_names
  }
  if(!is.null(col_names)) {
    colnames(all_dense) <- col_names
  }

  return(all_dense)

}
