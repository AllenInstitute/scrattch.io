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

#' Save a both exon and intron counts to an HDF5 file (dttch)
#'
#' @param exon_mat The exon matrix to store in dgCMatrix format
#' @param intron_mat The intron matrix to store in dgCMatrix format
#' @param out_file The HDF5 file to write to
#' @param cols_are Specifies whether columns in the matrix are sample_ids or genes
#' @param overwrite Whether or not to overwrite an existing out_file
#' @param compression_level The data compression level for large HDF5 matrix objects. default = 4.
#'
write_dttch_data <- function(exon_mat = NULL,
                             intron_mat = NULL,
                             out_file = "counts.dttch",
                             cols_are = "sample_id",
                             overwrite = F,
                             compression_level = 4) {

  library(Matrix)
  library(rhdf5)

  if(is.null(exon_mat) & is.null(intron_mat)) {
    stop("Provide at least one of exon_mat or intron_mat.")
  }

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
    print("Columns are samples; Rows are genes. Matrices will be transposed.")
    if(!is.null(exon_mat)) {
      t_exon_mat <- exon_mat
      exon_mat <- Matrix::t(exon_mat)
    }

    if(!is.null(intron_mat)) {
      t_intron_mat <- intron_mat
      intron_mat <- Matrix::t(intron_mat)
    }

  } else {
    print("Columns are genes; Rows are samples. Matrices will be used as-is.")
    if(!is.null(exon_mat)) {
      t_exon_mat <- Matrix::t(exon_mat)
    }
    if(!is.null(intron_mat)) {
      t_intron_mat <- Matrix::t(intron_mat)
    }
  }

  root <- H5Fopen(out_file, flags = "H5F_ACC_RDWR")
  H5Fclose(root)

  ## Exon Data
  if(!is.null(exon_mat)) {
    h5createGroup(out_file, "exon")
    h5createGroup(out_file, "t_exon")
    suppressWarnings({
      # Rows = Samples, Columns = Genes (Fast gene retrieval)
      print("Writing exon/x.")
      # data values
      h5createDataset(out_file,
                      dataset = "exon/x",
                      dims = length(exon_mat@x),
                      chunk = 1000,
                      level = compression_level)
      h5write(exon_mat@x,
              out_file,
              "exon/x")

      # data indices
      print("Writing exon/i.")
      h5createDataset(out_file,
                      dataset = "exon/i",
                      dims = length(exon_mat@x),
                      chunk = 1000, level = compression_level)
      h5write(exon_mat@i,
              out_file,
              "exon/i")

      # data index pointers
      print("Writing exon/p.")
      h5write(exon_mat@p,
              out_file,
              "exon/p")

      h5write(c(nrow(exon_mat), ncol(exon_mat)),
              out_file,
              "exon/dim")

      # Rows = Genes, Columns = Samples (Fast sample retrieval)

      # t_data values
      print("Writing t_exon/x.")
      h5createDataset(out_file,
                      dataset = "t_exon/x",
                      dims = length(t_exon_mat@x),
                      chunk = 1000,
                      level = compression_level)
      h5write(t_exon_mat@x,
              out_file,
              "t_exon/x")

      # t_data indices
      print("Writing t_exon/i.")
      h5createDataset(out_file,
                      dataset = "t_exon/i",
                      dims = length(t_exon_mat@i),
                      chunk = 1000,
                      level = compression_level)
      h5write(t_exon_mat@i,
              out_file,
              "t_exon/i")

      # t_data index pointers
      print("Writing t_exon/p.")
      h5write(t_exon_mat@p,
              out_file,
              "t_exon/p")
      h5write(c(nrow(t_exon_mat), ncol(t_exon_mat)),
              out_file,
              "t_exon/dim")
    })
  }

  if(!is.null(intron_mat)) {
    h5createGroup(out_file, "intron")
    h5createGroup(out_file, "t_intron")

    suppressWarnings({

      ## Intron data
      # Rows = Samples, Columns = Genes (Fast gene retrieval)
      print("Writing intron/x.")
      # data values
      h5createDataset(out_file,
                      dataset = "intron/x",
                      dims = length(intron_mat@x),
                      chunk = 1000,
                      level = compression_level)
      h5write(intron_mat@x,
              out_file,
              "intron/x")

      # data indices
      print("Writing intron/i.")
      h5createDataset(out_file,
                      dataset = "intron/i",
                      dims = length(intron_mat@x),
                      chunk = 1000, level = compression_level)
      h5write(intron_mat@i,
              out_file,
              "intron/i")

      # data index pointers
      print("Writing intron/p.")
      h5write(intron_mat@p,
              out_file,
              "intron/p")

      h5write(c(nrow(intron_mat), ncol(intron_mat)),
              out_file,
              "intron/dim")

      # Rows = Genes, Columns = Samples (Fast sample retrieval)

      # t_data values
      print("Writing t_intron/x.")
      h5createDataset(out_file,
                      dataset = "t_intron/x",
                      dims = length(t_intron_mat@x),
                      chunk = 1000,
                      level = compression_level)
      h5write(t_intron_mat@x,
              out_file,
              "t_intron/x")

      # t_data indices
      print("Writing t_intron/i.")
      h5createDataset(out_file,
                      dataset = "t_intron/i",
                      dims = length(t_intron_mat@i),
                      chunk = 1000,
                      level = compression_level)
      h5write(t_intron_mat@i,
              out_file,
              "t_intron/i")

      # t_data index pointers
      print("Writing t_intron/p.")
      h5write(t_intron_mat@p,
              out_file,
              "t_intron/p")
      h5write(c(nrow(t_intron_mat), ncol(t_intron_mat)),
              out_file,
              "t_intron/dim")
    })
  }


    # genes and sample_ids
    print("Writing gene_names.")
    if(!is.null(exon_mat)) {
      h5write(colnames(exon_mat),
              out_file,
              "gene_names")
    } else {
      h5write(colnames(intron_mat),
              out_file,
              "gene_names")
    }


    print("Writing sample_names.")
    if(!is.null(exon_mat)) {
      h5write(rownames(exon_mat),
              out_file,
              "sample_names")
    } else {
      h5write(rownames(intron_mat),
              out_file,
              "sample_names")
    }

    if(!is.null(exon_mat)) {
      print("Calculating total exon counts per sample")
      total_exon_counts <- unname(unlist(apply(t_exon_mat, 2, sum)))
      print("Writing total_exon_counts.")
      h5write(total_exon_counts,
              out_file,
              "total_exon_counts")
    }

    if(!is.null(intron_mat)) {
      print("Calculating total intron counts per sample")
      total_intron_counts <- unname(unlist(apply(t_intron_mat, 2, sum)))
      print("Writing total_intron_counts.")
      h5write(total_intron_counts,
              out_file,
              "total_intron_counts")
    }

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

#' Read Gene Expression Data from a dttch file
#'
#' @param dttch_file dttch file to read
#' @param genes A vector of gene names to read. If NULL, will read all genes.
#' @param regions The gene regions to use. Can be "exon", "intron", or "both". Default = "exon".
#' @param values The type of values to return. Can be "counts" or "cpm". Default = "counts".
#' @param transform Transformation to apply to values. Can be "none", "log", "log2", "log10". Each of these will add 1 to values before transformation. Default = "none".
#'
#' @return A data.frame with sample_name as the first column and each subsequent column
#' containing gene expression values and named for the genes.
#'
read_dttch_gene_data <- function(dttch_file,
                                 genes,
                                 regions = "exon",
                                 type = "counts",
                                 transform = "none") {
  library(rhdf5)
  library(purrr)
  library(dplyr)

  root <- H5Fopen(dttch_file)
  gene_names <- h5read(root,"/gene_names")
  sample_names <- h5read(root,"/sample_names")

  gene_index <- match(genes, gene_names)

  ## Exon values
  if(regions == "exon" | regions == "both") {

    exon_starts <- h5read(root, "/exon/p")[gene_index] + 1
    exon_ends <- h5read(root, "/exon/p")[(gene_index + 1)]

    exon_values <- map2(exon_starts, exon_ends, function(start, end) {
      if(end > start) {
        values <- h5read(root, "/exon/x", index = list(start:end))
      } else {
        values <- NA
      }
      values
    })

    exon_sample_indexes <- map2(exon_starts, exon_ends, function(start, end) {
      if(end > start) {
        values <- h5read(root, "/exon/i", index = list(start:end))
      } else {
        values <- NA
      }
      values
    })

    exon_sample_names <- map(exon_sample_indexes, function(exon_indexes) {
      if(length(exon_indexes) == 1) {
        if(!is.na(exon_indexes)) {
          sample_names <- NA
        } else {
          # R-like indexing; Previous Python-like indexing +1
          sample_names <- sample_names[exon_indexes + 1]
        }
      } else {
        sample_names <- sample_names[exon_indexes + 1]
      }
      sample_names
    })

    exon_df <- data.frame(sample_name = sample_names)

    exon_dfs <- pmap(list(x = exon_sample_names,
                          y = exon_values,
                          z = genes),
                     function(x,
                              y,
                              z) {
                       if(!is.na(exon_sample_names)[1]) {
                         df <- data.frame(sample_name = x,
                                          values = y)
                         names(df)[2] <- z
                       } else {
                         df <- NA
                       }
                       df
                     })



    suppressWarnings({
      exon_df <- reduce(c(list(exon_df), exon_dfs),
                        left_join,
                        by = "sample_name")
    })

    if(type == "cpm") {
      total_exon_counts <- h5read(root, "/total_exon_counts")
      exon_df[,genes] <- exon_df[,genes]/(total_exon_counts/1e6)
    }

    exon_df[is.na(exon_df)] <- 0

  }

  ## Intron values
  if(regions == "intron" | regions == "both") {

    intron_starts <- h5read(root, "/intron/p")[gene_index] + 1
    intron_ends <- h5read(root, "/intron/p")[(gene_index + 1)]

    intron_values <- map2(intron_starts, intron_ends, function(start, end) {
      if(end > start) {
        values <- h5read(root, "/intron/x", index = list(start:end))
      } else {
        values <- NA
      }
      values
    })

    intron_sample_indexes <- map2(intron_starts, intron_ends, function(start, end) {
      if(end > start) {
        values <- h5read(root, "/intron/i", index = list(start:end))
      } else {
        values <- NA
      }
      values
    })

    intron_sample_names <- map(intron_sample_indexes, function(intron_indexes) {
      if(length(intron_indexes) == 1) {
        if(!is.na(intron_indexes)) {
          sample_names <- NA
        } else {
          # R-like indexing; Previous Python-like indexing +1
          sample_names <- sample_names[intron_indexes + 1]
        }
      } else {
        sample_names <- sample_names[intron_indexes + 1]
      }
      sample_names
    })

    intron_df <- data.frame(sample_name = sample_names)

    intron_dfs <- pmap(list(x = intron_sample_names,
                            y = intron_values,
                            z = genes),
                       function(x,
                                y,
                                z) {
                         if(!is.na(intron_sample_names)[1]) {
                           df <- data.frame(sample_name = x,
                                            values = y)
                           names(df)[2] <- z
                         } else {
                           df <- NA
                         }
                         df
                       })



    suppressWarnings({
      intron_df <- reduce(c(list(intron_df), intron_dfs),
                          left_join,
                          by = "sample_name")
    })

    if(type == "cpm") {
      total_intron_counts <- h5read(root, "/total_intron_counts")
      intron_df[,genes] <- intron_df[,genes]/(total_intron_counts/1e6)
    }

    intron_df[is.na(intron_df)] <- 0

  }

  H5Fclose(root)

  if(regions == "exon") {
    out_df <- exon_df
  } else if(regions == "intron") {
    out_df <- intron_df
  } else if(regions == "both") {
    out_df <- exon_df
    out_df[,genes] <- out_df[,genes] + intron_df[,genes]
    out_df
  }

  if(transform == "log") {
    out_df[,genes] <- log(out_df[,genes] + 1)
  } else if(transform == "log2") {
    out_df[,genes] <- log2(out_df[,genes] + 1)
  } else if(transform == "log10") {
    out_df[,genes] <- log10(out_df[,genes] + 1)
  }

  out_df

}

#' Read Sample Expression Data from a dttch file
#'
#' @param dttch_file dttch file to read
#' @param samples A vector of sample names to read. If NULL, will read all samples.
#' @param regions The gene regions to use. Can be "exon", "intron", or "both". Default = "exon".
#' @param values The type of values to return. Can be "counts" or "cpm". Default = "counts".
#' @param transform Transformation to apply to values. Can be "none", "log", "log2", "log10". Each of these will add 1 to values before transformation. Default = "none".
#'
#' @return A data.frame with gene_name as the first column and each subsequent column
#' containing gene expression values and named for the samples.
#'
read_dttch_sample_data <- function(dttch_file,
                                   samples,
                                   regions = "exon",
                                   type = "counts",
                                   transform = "none") {
  library(rhdf5)
  library(purrr)
  library(dplyr)

  root <- H5Fopen(dttch_file)
  sample_names <- h5read(root,"/sample_names")
  gene_names <- h5read(root,"/gene_names")

  sample_index <- match(samples, sample_names)

  ## Exon values
  if(regions == "exon" | regions == "both") {

    exon_starts <- h5read(root, "/t_exon/p")[sample_index] + 1
    exon_ends <- h5read(root, "/t_exon/p")[(sample_index + 1)]

    exon_values <- map2(exon_starts, exon_ends, function(start, end) {
      if(end > start) {
        values <- h5read(root, "/t_exon/x", index = list(start:end))
      } else {
        values <- NA
      }
      values
    })

    exon_gene_indexes <- map2(exon_starts, exon_ends, function(start, end) {
      if(end > start) {
        values <- h5read(root, "/t_exon/i", index = list(start:end))
      } else {
        values <- NA
      }
      values
    })

    exon_gene_names <- map(exon_gene_indexes, function(exon_indexes) {
      if(length(exon_indexes) == 1) {
        if(!is.na(exon_indexes)) {
          gene_ids <- NA
        } else {
          # R-like indexing; Previous Python-like indexing +1
          gene_ids <- gene_names[exon_indexes + 1]
        }
      } else {
        gene_ids <- gene_names[exon_indexes + 1]
      }
      gene_ids
    })

    exon_df <- data.frame(gene_id = gene_names)

    exon_dfs <- pmap(list(x = exon_gene_names,
                          y = exon_values,
                          z = samples),
                     function(x,
                              y,
                              z) {
                       if(!is.na(exon_gene_names)[1]) {
                         df <- data.frame(gene_id = x,
                                          values = y)
                         names(df)[2] <- z
                       } else {
                         df <- NA
                       }
                       df
                     })



    suppressWarnings({
      exon_df <- reduce(c(list(exon_df), exon_dfs),
                        left_join,
                        by = "gene_id")
    })

    if(type == "cpm") {
      total_exon_counts <- h5read(root, "/total_exon_counts")[sample_index]
      exon_df[,samples] <- sweep(exon_df[,samples], 2, (total_exon_counts/1e6), "/")
    }

    exon_df[is.na(exon_df)] <- 0

  }

  ## Intron values
  if(regions == "intron" | regions == "both") {

    intron_starts <- h5read(root, "/t_intron/p")[sample_index] + 1
    intron_ends <- h5read(root, "/t_intron/p")[(sample_index + 1)]

    intron_values <- map2(intron_starts, intron_ends, function(start, end) {
      if(end > start) {
        values <- h5read(root, "/t_intron/x", index = list(start:end))
      } else {
        values <- NA
      }
      values
    })

    intron_gene_indexes <- map2(intron_starts, intron_ends, function(start, end) {
      if(end > start) {
        values <- h5read(root, "/t_intron/i", index = list(start:end))
      } else {
        values <- NA
      }
      values
    })

    intron_gene_names <- map(intron_gene_indexes, function(intron_indexes) {
      if(length(intron_indexes) == 1) {
        if(!is.na(intron_indexes)) {
          gene_ids <- NA
        } else {
          # R-like indexing; Previous Python-like indexing +1
          gene_ids <- gene_names[intron_indexes + 1]
        }
      } else {
        gene_ids <- gene_names[intron_indexes + 1]
      }
      gene_ids
    })

    intron_df <- data.frame(gene_id = gene_names)

    intron_dfs <- pmap(list(x = intron_gene_names,
                            y = intron_values,
                            z = samples),
                       function(x,
                                y,
                                z) {
                         if(!is.na(intron_gene_names)[1]) {
                           df <- data.frame(gene_id = x,
                                            values = y)
                           names(df)[2] <- z
                         } else {
                           df <- NA
                         }
                         df
                       })



    suppressWarnings({
      intron_df <- reduce(c(list(intron_df), intron_dfs),
                          left_join,
                          by = "gene_id")
    })

    if(type == "cpm") {
      total_intron_counts <- h5read(root, "/total_intron_counts")[sample_index]
      intron_df[,samples] <- sweep(intron_df[,samples], 2, (total_intron_counts/1e6), "/")
    }

    intron_df[is.na(intron_df)] <- 0

  }

  H5Fclose(root)

  if(regions == "exon") {
    out_df <- exon_df
  } else if(regions == "intron") {
    out_df <- intron_df
  } else if(regions == "both") {
    out_df <- exon_df
    out_df[,samples] <- out_df[,samples] + intron_df[,samples]
    out_df
  }

  if(transform == "log") {
    out_df[,samples] <- log(out_df[,samples] + 1)
  } else if(transform == "log2") {
    out_df[,samples] <- log2(out_df[,samples] + 1)
  } else if(transform == "log10") {
    out_df[,samples] <- log10(out_df[,samples] + 1)
  }

  out_df

}
