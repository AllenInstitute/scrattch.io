#' Read a data.frame from a tome file
#'
#' @param tome tome file to read
#' @param df_name character, the name of the data.frame object in the tome file structure
#' @param stored_as character, the storage mode of the data.frame used in write_tome_data.frame. Default = "data.frame".
#' If "data.frame", will read the df_name as a compound object. If "vectors" will read separate column vectors and build a data.frame.
#' @param columns character vector, specific columns to read (only if the object was stored as vectors).
#' @param get_all logical, whether or not to append all other columns after the specified columns (only if the object was stored as vectors).
#'
read_tome_data.frame <- function(tome,
                                 df_name,
                                 stored_as = "data.frame",
                                 columns = NULL,
                                 match_type = "exact",
                                 get_all = FALSE) {

  library(rhdf5)
  library(purrr)

  H5close()

  ls <- h5ls(tome)

  if(stored_as == "vectors") {
    all_columns <- ls$name[ls$group == df_name]

    # Filter cols if columns are provided.
    if(!is.null(columns)) {
      if(match_type == "grep") {
        if(length(columns) > 1) {
          column_pattern <- paste(columns, collapse = "|")
        } else {
          column_pattern <- columns
        }
        selected_columns <- all_columns[grepl(column_pattern, all_columns)]
      } else if(match_type == "exact") {
        selected_columns <- all_columns[all_columns %in% columns]
      }


      # If get_all, get the selected columns first, then all of the others
      if(get_all) {
        selected_columns <- c(selected_columns, setdiff(all_columns, selected_columns))
      }

    } else {
      selected_columns <- all_columns
    }

    df <- map(selected_columns,
              function(x) {
                h5read(tome, paste0(df_name,"/",x))
              }
    )
    names(df) <- selected_columns
    df <- as.data.frame(df)
  } else if(stored_as == "data.frame") {
    df <- h5read(tome, df_name)
  }

  H5close()

  df

}

#' Read a serialized object from a tome file
#'
#' @param tome tome file to read
#' @param target character, the name of the serialized object in the tome file structure
#'
read_tome_serialized <- function(tome,
                                 target) {

  H5close()

  serial_obj <- h5read(tome,
                       target)

  obj <- unserialize(charToRaw(serial_obj))

  H5close()

  obj
}


#' Read tome gene count data as a jagged list
#'
#' @param tome the tome file to read.
#' @param genes a character vector of genes to read from the tome.
#' @param regions Which regions to retrieve. Can be "exon", "intron", or "both".
#'
read_tome_genes_jagged <- function(tome,
                                   genes,
                                   regions = "exon") {
  library(rhdf5)
  library(purrr)
  library(dplyr)
  library(Matrix)

  H5close()

  root <- H5Fopen(tome)
  gene_names <- h5read(root,"/gene_names")
  sample_names <- h5read(root,"/sample_names")

  gene_index <- match(genes, gene_names)

  ## Exon values
  if(regions == "exon" | regions == "both") {

    exon_starts <- h5read(root, "data/exon/p")[gene_index] + 1
    exon_ends <- h5read(root, "data/exon/p")[(gene_index + 1)]

    exon_values <- map2(exon_starts, exon_ends, function(start, end) {
      if(end > start) {
        values <- h5read(root, "data/exon/x", index = list(start:end))
      } else {
        values <- NA
      }
      values
    })

    exon_sample_indexes <- map2(exon_starts, exon_ends, function(start, end) {
      if(end > start) {
        values <- h5read(root, "data/exon/i", index = list(start:end))
      } else {
        values <- NA
      }
      values
    })

    exon <- list(x = exon_values,
                 i = exon_sample_indexes,
                 p = c(0, cumsum(map_int(exon_values, length))))


  }

  ## Intron values
  if(regions == "intron" | regions == "both") {

    intron_starts <- h5read(root, "data/intron/p")[gene_index] + 1
    intron_ends <- h5read(root, "data/intron/p")[(gene_index + 1)]

    intron_values <- map2(intron_starts, intron_ends, function(start, end) {
      if(end > start) {
        values <- h5read(root, "data/intron/x", index = list(start:end))
      } else {
        values <- NA
      }
      values
    })

    intron_sample_indexes <- map2(intron_starts, intron_ends, function(start, end) {
      if(end > start) {
        values <- h5read(root, "data/intron/i", index = list(start:end))
      } else {
        values <- NA
      }
      values
    })

    intron <- list(x = intron_values,
                   i = intron_sample_indexes,
                   p = c(0, cumsum(map_int(intron_values, length))))

  }

  H5Fclose(root)

  if(regions == "exon") {
    out <- exon
  } else if(regions == "intron") {
    out <- intron
  } else if(regions == "both") {
    out <- exon
    xi <- pmap(list(exon_x = exon$x,
                    exon_i = exon$i,
                    intron_x = intron$x,
                    intron_i = intron$i),
               function(exon_x, exon_i,
                        intron_x, intron_i) {
                 both_i <- union(exon_i, intron_i)
                 both_i <- sort(both_i)
                 both_x <- numeric(length(both_i))
                 both_x[match(exon_i, both_i)] <- both_x[match(exon_i, both_i)] + exon_x
                 both_x[match(intron_i, both_i)] <- both_x[match(intron_i, both_i)] + intron_x

                 list(x = both_x,
                      i = both_i)
               })
    for(j in 1:length(xi)) {
      out$x[[j]] <- xi[[j]]$x
      out$i[[j]] <- xi[[j]]$i
    }
    out$p <- c(0, cumsum(map_int(out$x, length)))
  }

  out$gene_names <- genes
  out$sample_names <- read_tome_sample_names(tome)
  out$dims <- read_tome_data_dims(tome)
  out$dims[2] <- length(genes)

  H5close()

  out

}

#' Read tome sample count data as a jagged list
#'
#' @param tome the tome file to read.
#' @param samples a character vector of genes to read from the tome.
#' @param regions Which regions to retrieve. Can be "exon", "intron", or "both".
#'
read_tome_samples_jagged <- function(tome,
                                     samples,
                                     regions = "exon") {
  library(rhdf5)
  library(purrr)
  library(dplyr)
  library(Matrix)

  H5close()

  root <- H5Fopen(tome)
  gene_names <- h5read(root,"/gene_names")
  sample_names <- h5read(root,"/sample_names")

  sample_index <- match(samples, sample_names)

  ## Exon values
  if(regions == "exon" | regions == "both") {

    exon_starts <- h5read(root, "data/t_exon/p")[sample_index] + 1
    exon_ends <- h5read(root, "data/t_exon/p")[(sample_index + 1)]

    exon_values <- map2(exon_starts, exon_ends, function(start, end) {
      if(end > start) {
        values <- h5read(root, "data/t_exon/x", index = list(start:end))
      } else {
        values <- NA
      }
      values
    })

    exon_sample_indexes <- map2(exon_starts, exon_ends, function(start, end) {
      if(end > start) {
        values <- h5read(root, "data/t_exon/i", index = list(start:end))
      } else {
        values <- NA
      }
      values
    })

    exon <- list(x = exon_values,
                 i = exon_sample_indexes,
                 p = c(0, cumsum(map_int(exon_values, length))))


  }

  ## Intron values
  if(regions == "intron" | regions == "both") {

    intron_starts <- h5read(root, "data/t_intron/p")[sample_index] + 1
    intron_ends <- h5read(root, "data/t_intron/p")[(sample_index + 1)]

    intron_values <- map2(intron_starts, intron_ends, function(start, end) {
      if(end > start) {
        values <- h5read(root, "data/t_intron/x", index = list(start:end))
      } else {
        values <- NA
      }
      values
    })

    intron_sample_indexes <- map2(intron_starts, intron_ends, function(start, end) {
      if(end > start) {
        values <- h5read(root, "data/t_intron/i", index = list(start:end))
      } else {
        values <- NA
      }
      values
    })

    intron <- list(x = intron_values,
                   i = intron_sample_indexes,
                   p = c(0, cumsum(map_int(intron_values, length))))

  }

  H5Fclose(root)

  if(regions == "exon") {
    out <- exon
  } else if(regions == "intron") {
    out <- intron
  } else if(regions == "both") {
    out <- exon
    xi <- pmap(list(exon_x = exon$x,
                    exon_i = exon$i,
                    intron_x = intron$x,
                    intron_i = intron$i),
               function(exon_x, exon_i,
                        intron_x, intron_i) {
                 both_i <- union(exon_i, intron_i)
                 both_i <- sort(both_i)
                 both_x <- numeric(length(both_i))
                 both_x[match(exon_i, both_i)] <- both_x[match(exon_i, both_i)] + exon_x
                 both_x[match(intron_i, both_i)] <- both_x[match(intron_i, both_i)] + intron_x

                 list(x = both_x,
                      i = both_i)
               })
    for(j in 1:length(xi)) {
      out$x[[j]] <- xi[[j]]$x
      out$i[[j]] <- xi[[j]]$i
    }
    out$p <- c(0, cumsum(map_int(out$x, length)))
  }

  out$gene_names <- read_tome_gene_names(tome)
  out$sample_names <- samples
  out$dims <- read_tome_data_dims(tome, transpose = TRUE)
  out$dims[2] <- length(samples)

  H5close()

  out

}

#' Convert a jagged list of gene counts to a matrix
#'
#' @param jagged The jagged list object to convert.
#' @param rows Character, either "sample_names" or "gene_names".
#' @param cols Character, either "sample_names" or "gene_names".
#'
jagged_to_matrix <- function(jagged,
                             rows = c("sample_names","gene_names"),
                             cols = c("gene_names", "sample_names")) {

  out <- matrix(0,
                nrow = jagged$dims[1],
                ncol = jagged$dims[2])
  rownames(out) <- jagged[[rows]]
  colnames(out) <- jagged[[cols]]

  for(j in 1:length(jagged$x)) {
    x <- jagged$x[[j]]
    i <- jagged$i[[j]] + 1
    col <- rep(j, length(x))
    out[i, col] <- x
  }

  out
}

#' Convert a jagged list of gene counts to a data.frame
#'
#' @param jagged The jagged list object to convert.
#' @param rows Character, either "sample_names" or "gene_names".
#' @param cols Character, either "sample_names" or "gene_names".
#'
jagged_to_data.frame <- function(jagged,
                                 rows = c("sample_names","gene_names"),
                                 cols = c("gene_names", "sample_names")) {
  if(rows == "sample_names") {
    out <- data.frame(sample_name = jagged$sample_names,
                      matrix(0,
                             nrow = jagged$dims[1],
                             ncol = jagged$dims[2]))
    names(out)[-1] <- jagged$gene_names

  } else {
    out <- data.frame(gene_name = jagged$gene_names,
                      matrix(0,
                             nrow = jagged$dims[1],
                             ncol = jagged$dims[2]))
    names(out)[-1] <- jagged$sample_names

  }

  for(j in 1:length(jagged$x)) {
    x <- jagged$x[[j]]
    i <- jagged$i[[j]] + 1
    out[i, j + 1] <- x
  }

  out
}

#' Convert a jagged list of gene counts to a sparse, dgCMatrix
#'
#' @param jagged The jagged list object to convert.
#' @param rows Character, either "sample_names" or "gene_names".
#' @param cols Character, either "sample_names" or "gene_names".
#'
jagged_to_dgCMatrix <- function(jagged,
                                rows = c("sample_names","gene_names"),
                                cols = c("gene_names", "sample_names")) {

  sparseMatrix(i = unlist(jagged$i),
               p = jagged$p,
               x = unlist(jagged$x),
               index1 = FALSE,
               dims = jagged$dims,
               dimnames = list(jagged[[rows]],
                               jagged[[cols]]))

}
