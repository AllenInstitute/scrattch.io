#' Convert a large matrix object to a dgCMatrix
#'
#' May be useful in cases where calling Matrix() directly gives a "negative length vectors are not allowed" error.
#'
#' This works through the matrix in column-based chunks becaues this is much faster for conversion to column-indexed dgCMatrices.
#'
#' @param mat The matrix file to convert
#' @param chunk_size The number of columns to read as a chunk. For ~30k genes, a chunk of 5000 (the default) is ~1GB in memory.
#'
#' @return A dgCMatrix object in the same orientation as the matrix object.
#'
large_matrix_to_dgCMatrix <- function(mat,
                                      chunk_size = 5000) {
  library(Matrix)

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
  all_sparse <- Matrix(chunk_mat, sparse = T)

  # chunkation over chunks
  for(chunk in 2:n_chunks) {
    chunk_start <- 1 + chunk_size * (chunk - 1)
    chunk_end <- chunk_size * chunk
    print(paste0("Adding rows ",chunk_start," to ", chunk_end))
    chunk_mat <- mat[,chunk_start:chunk_end]
    chunk_sparse <- Matrix(chunk_mat, sparse = T)
    all_sparse <- cbind(all_sparse, chunk_mat)
  }

  # Remaining samples
  chunk_mat <- mat[,(n_chunks*chunk_size + 1):n_col]
  chunk_sparse <- Matrix(chunk_mat, sparse = T)
  all_sparse <- cbind(all_sparse, chunk_mat)

  rownames(all_sparse) <- row_names
  colnames(all_sparse) <- col_names

  return(all_sparse)

}
