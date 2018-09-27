#' Read a Loom matrix as a dgCMatrix
#'
#' @param loom_file The loom file to read
#' @param chunk_size The number of rows to read as a chunk. For ~30k genes, a chunk of 5000 (the default) is ~1GB in memory.
#'
#' @return A dgCMatrix object with genes as columns and samples as rows.
#'
read_loom_dgCMatrix <- function(loom_file,
                                chunk_size = 5000,
                                row_names = "Gene",
                                col_names = "CellID") {
  library(rhdf5)
  library(Matrix)

  # Gene names are in /row_attrs/Gene
  gene_names <- read_tome_vector(loom_file, paste0("/row_attrs/", row_names))

  # Sample names are in col_attrs/CellID
  sample_names <- read_tome_vector(loom_file, paste0("/col_attrs/", col_names))

  n_samples <- length(sample_names)
  n_genes <- length(gene_names)

  if(n_samples < chunk_size) {
    chunk_size <- n_samples
    n_chunks <- 1
  } else {
    n_chunks <- floor(n_samples/chunk_size)
  }

  # Initial chunk
  chunk <- 1
  chunk_start <- 1 + chunk_size * (chunk - 1)
  chunk_end <- chunk_size * chunk
  cat(paste0("Reading samples ",chunk_start," to ", chunk_end))
  chunk_mat <- rhdf5::h5read(loom_file, "/matrix", index = list(chunk_start:chunk_end, 1:n_genes))
  all_sparse <- Matrix::Matrix(chunk_mat, sparse = T)

  if(n_chunks > 1) {

    # chunkation over chunks
    for(chunk in 2:n_chunks) {
      chunk_start <- 1 + chunk_size * (chunk - 1)
      chunk_end <- chunk_size * chunk
      cat(paste0("Reading samples ",chunk_start," to ", chunk_end))
      chunk_mat <- rhdf5::h5read(loom_file, "/matrix", index = list(chunk_start:chunk_end, 1:n_genes))
      chunk_sparse <- Matrix::Matrix(chunk_mat, sparse = T)
      all_sparse <- rbind(all_sparse, chunk_mat)
    }

    # Remaining samples
    chunk_mat <- rhdf5::h5read(loom_file, "/matrix", index = list((n_chunks*chunk_size + 1):n_samples, 1:n_genes))
    chunk_sparse <- Matrix::Matrix(chunk_mat, sparse = T)
    all_sparse <- rbind(all_sparse, chunk_mat)
  }

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
read_loom_anno <- function(loom_file,
                           sample_col = "CellID") {
  library(rhdf5)

  # annotations are stored in /col_attrs (Column attributes)
  anno <- read_tome_data.frame(loom_file, "/col_attrs", stored_as = "vectors")

  # projections are also stored here, and start with X_
  projection_columns <- names(anno)[grepl("^X_",names(anno))]

  anno <- anno %>%
    dplyr::select(-one_of(projection_columns)) %>%
    dplyr::rename_("sample_name" = sample_col) %>%
    dplyr::select(sample_name, dplyr::everything())

  return(anno)

}

#' Read Loom projections
#'
#' @param loom_file The loom file to read
#'
#' @return A data.frame with projection values as columns and samples as rows. The Loom CellID becomes the first column, sample_name.
#' Loom doesn't use a prefix to indicate paired coordinates, so you'll have to figure these out on your own.
#'
read_loom_projections <- function(loom_file,
                                  sample_col = "CellID") {
  library(rhdf5)

  # projtations are stored in /col_attrs (Column attributes)
  proj <- read_tome_data.frame(loom_file, "/col_attrs", stored_as = "vectors")

  # projections are also stored here, and start with X_
  projection_columns <- names(proj)[grepl("^X_",names(proj))]

  proj <- proj %>%
    dplyr::select(one_of(c(sample_col, projection_columns))) %>%
    dplyr::rename_("sample_name" = sample_col) %>%
    dplyr::select(sample_name, dplyr::everything())

  names(proj) <- sub("^X_","",names(proj))

  return(proj)
}
