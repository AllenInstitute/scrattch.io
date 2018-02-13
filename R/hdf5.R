#' Save a dgCMatrix to HDF5 format
#'
#' @param mat The matrix to be written
#' @param out_file The HDF5 file to write to
#' @param cols_are Specifies whether columns in the matrix are sample_ids or genes
#' @param overwrite Whether or not to overwrite an existing out_file
#' @param compression_level The data compression level for large HDF5 matrix objects
#'
save_sparse_matrix_h5 <- function(mat,
                                  out_file,
                                  cols_are = "sample_id",
                                  overwrite = F,
                                  compression_level = 0) {

  library(rhdf5)

  # Modified from cellrangerRkit
  H5close()
  if (file.exists(out_file) & overwrite) {
    # Delete old version and make a new file
    unlink(out_file)
    h5createFile(out_file)

  } else if (file.exists(out_file) & !overwrite) {
    # Stop if overwrite = F
    stop(paste0(out_file," exists. Set overwrite = TRUE to overwrite."))
    return()

  } else if(!file.exists(out_file)) {
    # If the file doesn't exist, make a new file
    h5createFile(out_file)

  }

  if(grepl("sample",cols_are) | grepl("cell",cols_are)) {
    print("Columns are samples; Rows are genes. Matrix will be transposed.")
    tmat <- mat
    mat <- t(mat)
  } else {
    print("Columns are genes; Rows are samples. Matrix will be used as-is.")
    tmat <- t(mat)
  }

  root <- H5Fopen(out_file, flags = "H5F_ACC_RDWR")
  H5Fclose(root)

  h5createGroup(out_file, "data")
  h5createGroup(out_file, "t_data")

  suppressWarnings({
    # Rows = Samples, Columns = Genes (Fast sample retrieval)
    print("Writing data/x.")
    # data values
    h5createDataset(out_file,
                    dataset = "data/x",
                    dims = length(mat@x),
                    chunk = 1000,
                    level = compression_level)
    h5write(mat@x,
            out_file,
            "data/x")

    # data indices
    print("Writing data/i.")
    h5createDataset(out_file,
                    dataset = "data/i",
                    dims = length(mat@x),
                    chunk = 1000, level = compression_level)
    h5write(mat@i,
            out_file,
            "data/i")

    # data index pointers
    print("Writing data/p.")
    h5write(mat@p,
            out_file,
            "data/p")

    h5write(c(nrow(mat), ncol(mat)),
            out_file,
            "data/dim")

    # Rows = Genes, Columns = Samples (Fast gene retrieval)

    # t_data values
    print("Writing t_data/x.")
    h5createDataset(out_file,
                    dataset = "t_data/x",
                    dims = length(tmat@x),
                    chunk = 1000,
                    level = compression_level)
    h5write(tmat@x,
            out_file,
            "t_data/x")

    # t_data indices
    print("Writing t_data/i.")
    h5createDataset(out_file,
                    dataset = "t_data/i",
                    dims = length(tmat@i),
                    chunk = 1000,
                    level = compression_level)
    h5write(tmat@i,
            out_file,
            "t_data/i")

    # t_data index pointers
    print("Writing t_data/p.")
    h5write(tmat@p,
            out_file,
            "t_data/p")
    h5write(c(nrow(tmat), ncol(tmat)),
            out_file,
            "t_data/dim")

    # genes and sample_ids
    print("Writing gene_names.")
    h5write(colnames(mat),
            out_file,
            "gene_names")

    print("Writing sample_names.")
    h5write(rownames(mat),
            out_file,
            "sample_names")

    print("Calculating total Counts per sample")
    total_counts <- unname(unlist(apply(tmat, 2, sum)))
    print("Writing total_counts.")
    h5write(total_counts,
            out_file,
            "total_counts")

  })

  H5close()

}


#' Read Gene Expression Data from an HDF5 file
#'
#' @param hdf5_file HDF5 file to read
#' @param genes A vector of gene names to read
#'
#' @return A data.frame with sample_id as the first column and each subsequent column
#' containing gene expression values and named for the genes.
#'
read_gene_data_hdf5 <- function(hdf5_file,
                                genes,
                                type = c("counts","cpm")) {
  library(rhdf5)
  library(purrr)
  library(dplyr)

  root <- H5Fopen(hdf5_file)
  gene_index <- match(genes, h5read(root,"/genes"))

  gene_starts <- h5read(root, "/data/indptr")[gene_index]
  # Python-like indexing; -1
  gene_ends <- h5read(root, "/data/indptr")[(gene_index + 1)] - 1

  gene_values <- map2(gene_starts, gene_ends, function(start, end) {
    if(end > start) {
      values <- h5read(root, "/data/values", index = list(start:end))
    } else {
      values <- NA
    }
    values
  })

  gene_sample_indexes <- map2(gene_starts, gene_ends, function(start, end) {
    if(end > start) {
      values <- h5read(root, "/data/indices", index = list(start:end))
    } else {
      values <- NA
    }
    values
  })

  gene_samples <- map(gene_sample_indexes, function(gene_indexes) {
    if(length(gene_indexes) == 1) {
      if(!is.na(gene_indexes)) {
        sample_ids <- NA
      } else {
        # R-like indexing; Previous Python-like indexing +1
        sample_ids <- h5read(root, "/sample_ids")[gene_indexes + 1]
      }
    } else {
      sample_ids <- h5read(root, "/sample_ids")[gene_indexes + 1]
    }
    sample_ids
  })

  gene_dfs <- pmap(list(x = gene_samples,
                        y = gene_values,
                        z = genes),
                   function(x,
                            y,
                            z) {
                     if(!is.na(gene_samples)[1]) {
                       df <- data.frame(sample_id = x,
                                        values = y)
                       names(df)[2] <- z
                     } else {
                       df <- NA
                     }
                     df
                   })

  out_df <- data.frame(sample_id = h5read(root, "/sample_ids"))

  H5Fclose(root)

  suppressWarnings({
    out_df <- reduce(c(list(out_df), gene_dfs),
                     left_join,
                     by = "sample_id")
  })

  out_df[is.na(out_df)] <- 0

  out_df
}

#' Unlink an HDF5 object
#'
#' @param h5_file The HDF5 file to edit.
#' @param h5_object The full name of the object within the HDF5 file to unlink.
#'
h5_delete <- function(h5_file,
                      h5_object) {
  f <- h5::h5file(h5_file)
  h5::h5unlink(f, h5_object)
  h5::h5close(f)
}
