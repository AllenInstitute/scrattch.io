#' Read a Loom matrix as a dgCMatrix
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
  print(paste0("Reading samples ",chunk_start," to ", chunk_end))
  chunk_mat <- h5read(loom_file, "/matrix", index = list(chunk_start:chunk_end, 1:n_genes))
  all_sparse <- Matrix(chunk_mat, sparse = T)

  # chunkation over chunks
  for(chunk in 2:n_chunks) {
    chunk_start <- 1 + chunk_size * (chunk - 1)
    chunk_end <- chunk_size * chunk
    print(paste0("Reading samples ",chunk_start," to ", chunk_end))
    chunk_mat <- h5read(loom_file, "/matrix", index = list(chunk_start:chunk_end, 1:n_genes))
    chunk_sparse <- Matrix(chunk_mat, sparse = T)
    all_sparse <- rbind(all_sparse, chunk_mat)
  }

  # Remaining samples
  chunk_mat <- h5read(loom_file, "/matrix", index = list((n_chunks*chunk_size + 1):n_samples, 1:n_genes))
  chunk_sparse <- Matrix(chunk_mat, sparse = T)
  all_sparse <- rbind(all_sparse, chunk_mat)

  rownames(all_sparse) <- sample_names
  colnames(all_sparse) <- gene_names

  return(all_sparse)

}

#' Read Loom sample annotations
#'
#' @param loom_file The loom file to read
#'
#' @return A data.frame with annotations as columns and samples as rows. The Loom CellID becomes the first column, sample_name.
#'
read_loom_anno <- function(loom_file) {
  library(rhdf5)

  # annotations are stored in /col_attrs (Column attributes)
  anno <- read_tome_data.frame(loom_file, "/col_attrs", stored_as = "vectors")

  # projections are also stored here, and start with X_
  projection_columns <- names(anno)[grepl("^X_",names(anno))]

  anno <- anno %>%
    select(-one_of(projection_columns)) %>%
    rename_("sample_name" = "CellID") %>%
    select(sample_name, everything())

  return(anno)

}

#' Read Loom projections
#'
#' @param loom_file The loom file to read
#'
#' @return A data.frame with projection values as columns and samples as rows. The Loom CellID becomes the first column, sample_name.
#' Loom doesn't use a prefix to indicate paired coordinates, so you'll have to figure these out on your own.
#'
read_loom_projections <- function(loom_file) {
  library(rhdf5)

  # projtations are stored in /col_attrs (Column attributes)
  proj <- read_tome_data.frame(loom_file, "/col_attrs", stored_as = "vectors")

  # projections are also stored here, and start with X_
  projection_columns <- names(proj)[grepl("^X_",names(proj))]

  proj <- proj %>%
    select(one_of(c("CellID",projection_columns))) %>%
    rename_("sample_name" = "CellID") %>%
    select(sample_name, everything())

  names(proj) <- sub("^X_","",names(proj))

  return(proj)
}
