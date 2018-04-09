#' Reads a Loom matrix as a dgCMatrix
#'
#' @param loom_file The loom file to read
#' @param chunk_size The number of rows to read as a chunk. For ~30k genes, a chunk of 5000 is ~1GB in memory.
#'
#' @return A dgCMatrix object with genes as columns and samples as rows.
#'
read_loom_dgCMatrix <- function(loom_file,
                                chunk_size = 5000) {
  library(rhdf5)
  library(Matrix)

  # Gene names are in /row_attrs/Gene
  gene_names <- as.vector(unlist(h5read(loom_file, "/row_attrs/Gene")))

  # Sample names are in col_attrs/CellID
  sample_names <- as.vector(unlist(h5read(loom_file, "/col_attrs/CellID")))

  n_samples <- length(sample_names)
  n_genes <- length(gene_names)

  n_chunks <- floor(n_samples/chunk_size)

  # Initial chunk
  chunk <- 1
  chunk_start <- 1 + chunk_size * (chunk - 1)
  chunk_end <- chunk_size * chunk
  chunk_mat <- h5read(loom, "/matrix", index = list(chunk_start:chunk_end, 1:n_genes))
  all_sparse <- Matrix(chunk_mat, sparse = T)

  # chunkation over chunks
  for(chunk in 2:n_chunks) {
    chunk_start <- 1 + chunk_size * (chunk - 1)
    chunk_end <- chunk_size * chunk
    print(paste0("Reading samples ",chunk_start," to ", chunk_end))
    chunk_mat <- h5read(loom, "/matrix", index = list(chunk_start:chunk_end, 1:n_genes))
    chunk_sparse <- Matrix(chunk_mat, sparse = T)
    all_sparse <- rbind(all_sparse, chunk_mat)
  }

  # Remaining samples
  chunk_mat <- h5read(loom, "/matrix", index = list((n_chunks*chunk_size + 1):n_samples, 1:n_genes))
  chunk_sparse <- Matrix(chunk_mat, sparse = T)
  all_sparse <- rbind(all_sparse, chunk_mat)

  rownames(all_sparse) <- sample_names
  colnames(all_sparse) <- gene_names

  return(all_sparse)

}
