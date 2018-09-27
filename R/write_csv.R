#' Write a dgCMatrix to CSV in reasonable chunks
#'
#' This function makes it so that we don't have to explode a dgCMatrix to a full matrix
#' in order to write the matrix to a CSV file with a named first column derived from rownames(mat).
#'
#' @param mat the dgCMatrix to write
#' @param filename the target CSV file
#' @param col1_name the name of the first column to write
#' @param chunk_size the number of rows to retrieve per chunk. Lower = less RAM usage.
#'
write_dgCMatrix_csv <- function(mat,
                                filename,
                                col1_name = "gene",
                                chunk_size = 1000) {

  #library(Matrix)
  #library(data.table)

  # Transpose so retrieval of "rows" is much faster
  mat <- Matrix::t(mat)

  # Row names
  row_names <- colnames(mat)

  # gene names are now columns
  col_names <- rownames(mat)

  n_row <- length(row_names)
  n_col <- length(col_names)

  n_chunks <- floor(n_row/chunk_size)

  # Initial chunk
  chunk <- 1
  chunk_start <- 1 + chunk_size * (chunk - 1)
  chunk_end <- chunk_size * chunk
  print(paste0("Writing rows ",chunk_start," to ", chunk_end))
  chunk_mat <- t(as.matrix(mat[,chunk_start:chunk_end]))
  chunk_df <- cbind(data.frame(col1 = row_names[chunk_start:chunk_end]),as.data.frame(chunk_mat))
  names(chunk_df)[1] <- col1_name
  data.table::fwrite(chunk_df, file = filename, append = F)

  # chunkation over chunks
  for(chunk in 2:n_chunks) {
    chunk_start <- 1 + chunk_size * (chunk - 1)
    chunk_end <- chunk_size * chunk
    print(paste0("Writing rows ",chunk_start," to ", chunk_end))
    chunk_mat <- t(as.matrix(mat[,chunk_start:chunk_end]))
    chunk_df <- cbind(data.frame(col1 = row_names[chunk_start:chunk_end]),as.data.frame(chunk_mat))
    data.table::fwrite(chunk_df, file = filename, append = T)
  }

  # Remaining samples
  chunk_start <- (n_chunks*chunk_size + 1)
  chunk_end <- n_row
  print(paste0("Writing rows ",chunk_start," to ", chunk_end))
  chunk_mat <- t(as.matrix(mat[,chunk_start:chunk_end]))
  chunk_df <- cbind(data.frame(col1 = row_names[chunk_start:chunk_end]),as.data.frame(chunk_mat))
  data.table::fwrite(chunk_df, file = filename, append = T)

}
