#' Read a vector from a tome (or other HDF5 object)
#'
#' The h5read function sticks close to the actual object stored in an HDF5 file - that is,
#' it will return a 1d array instead of a true vector object. This function is a simple
#' wrapper around h5read that will apply unlist and as.vector to the output of h5read so that
#' the result is a standard R vector.
#'
#' @param tome A tome file (other HDF5 files will work as well).
#' @param name The name of the object in the HDF5 hierarchy.
#'
#' @return A vector of the same type as the object in the HDF5 file.
#'
read_tome_vector <- function(tome,
                             name) {
  as.vector(unlist(rhdf5::h5read(tome, name)))

}

#' Read a data.frame from a tome file
#'
#' @param tome tome file to read
#' @param df_name character, the name of the data.frame object in the tome file structure
#' @param stored_as character, the storage mode of the data.frame used in write_tome_data.frame. Default = "data.frame".
#' If "data.frame", will read the df_name as a compound object. If "vectors" will read separate column vectors and build a data.frame.
#' @param columns character vector, specific columns to read (only if the object was stored as vectors).
#' @param match_type Whether to match column names exactly ("exact") or using grep ("grep"). Default is "exact".
#' @param get_all logical, whether or not to append all other columns after the specified columns (only if the object was stored as vectors).
#'
read_tome_data.frame <- function(tome,
                                 df_name,
                                 stored_as = "data.frame",
                                 columns = NULL,
                                 match_type = "exact",
                                 get_all = FALSE) {

  ls <- rhdf5::h5ls(tome)

  if(stored_as == "vectors") {
    all_columns <- ls$name[ls$group == df_name]
  } else if (stored_as == "data.frame") {
    df <- rhdf5::h5read(tome, df_name)
    all_columns <- names(df)
  }
  # Filter cols if columns are provided.

  if(is.null(columns)) {
    selected_columns <- all_columns
  } else {
    if(match_type == "grep") {

      if(length(columns) > 1) {
        selected_columns <- purrr::map(columns,
                                       function(x) {
                                         all_columns[grepl(x, all_columns)]
                                       })

        unmatched_columns <- columns[map_int(selected_columns, length) == 0]

        if(length(unmatched_columns) > 0) {
          warning(paste("Warning: no match found for columns:",
                        paste(unmatched_columns, collapse = ", ")))
        }

        selected_columns <- unique(unlist(selected_columns))

      } else {
        selected_columns <- all_columns[grepl(columns, all_columns)]
      }

      # Stop if no matches found
      if(length(selected_columns) == 0) {
        stop("Error: No columns match the columns argument.")
      }

    } else if(match_type == "exact") {
      selected_columns <- all_columns[all_columns %in% columns]

      # Stop if no matches found
      if(length(selected_columns) == 0) {
        stop("Error: No columns match the columns argument.")
      }

      # Warn if some columns aren't matched
      unmatched_columns <- setdiff(columns, all_columns)
      if(length(unmatched_columns) > 0) {
        warning(paste("Warning: no match found for columns:",
                      paste(unmatched_columns, collapse = ", ")))
      }

    }

    # If get_all, get the selected columns first, then all of the others
    if(get_all) {
      selected_columns <- c(selected_columns, setdiff(all_columns, selected_columns))
    }

  }

  if(stored_as == "vectors") {

    df <- purrr::map_dfc(selected_columns,
                         function(x) {
                           read_tome_vector(tome, paste0(df_name,"/",x))
                         }
    )
    names(df) <- selected_columns

  } else if(stored_as == "data.frame") {
    df <- df[,selected_columns]
  }

  df

}

#' Read a serialized object from a tome file
#'
#' @param tome tome file to read
#' @param target character, the name of the serialized object in the tome file structure
#'
read_tome_serialized <- function(tome,
                                 target) {

  serial_obj <- rhdf5::h5read(tome,
                              target)

  obj <- unserialize(charToRaw(serial_obj))

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

  gene_names <- read_tome_vector(tome,"/gene_names")

  if(is.null(genes)) {
    genes <- gene_names
  } else {

    # Check for genes not found in the dataset.
    missing_genes <- setdiff(genes, gene_names)

    if(length(missing_genes) > 0) {
      stop("Some genes not found in tome: ", paste(missing_genes, collapse = ", "))
    }

  }

  sample_names <- read_tome_vector(tome,"/sample_names")

  gene_index <- match(genes, gene_names)

  ## Exon values
  if(regions == "exon" | regions == "both") {

    exon_starts <- read_tome_vector(tome, "data/exon/p")[gene_index] + 1
    exon_ends <- read_tome_vector(tome, "data/exon/p")[(gene_index + 1)]

    index_list <- unlist(map2(exon_starts, exon_ends, function(x,y) { x:y }))
    split_list <- unlist(map(1:length(exon_starts), function(x) { rep(x, length(exon_starts[x]:exon_ends[x])) }))

    exon_values <- rhdf5::h5read(tome, "data/exon/x", index = list(index_list))
    exon_values <- split(exon_values, split_list)
    names(exon_values) <- NULL

    exon_sample_indexes <- rhdf5::h5read(tome, "data/exon/i", index = list(index_list))
    exon_sample_indexes <- split(exon_sample_indexes, split_list)
    names(exon_sample_indexes) <- NULL

    exon_pointers <- c(0, cumsum(map_int(exon_values, length)))
    names(exon_pointers) <- NULL


    exon <- list(x = exon_values,
                 i = exon_sample_indexes,
                 p = exon_pointers)


  }

  ## Intron values
  if(regions == "intron" | regions == "both") {

    intron_starts <- read_tome_vector(tome, "data/intron/p")[gene_index] + 1
    intron_ends <- read_tome_vector(tome, "data/intron/p")[(gene_index + 1)]

    index_list <- unlist(map2(intron_starts, intron_ends, function(x,y) { x:y }))
    split_list <- unlist(map(1:length(intron_starts), function(x) { rep(x, length(intron_starts[x]:intron_ends[x])) }))

    intron_values <- rhdf5::h5read(tome, "data/intron/x", index = list(index_list))
    intron_values <- split(intron_values, split_list)
    names(intron_values) <- NULL

    intron_sample_indexes <- rhdf5::h5read(tome, "data/intron/i", index = list(index_list))
    intron_sample_indexes <- split(intron_sample_indexes, split_list)
    names(intron_sample_indexes) <- NULL

    intron_pointers <- c(0, cumsum(map_int(intron_values, length)))
    names(intron_pointers) <- NULL


    intron <- list(x = intron_values,
                 i = intron_sample_indexes,
                 p = intron_pointers)

  }

  if(regions == "exon") {
    out <- exon
  } else if(regions == "intron") {
    out <- intron
  } else if(regions == "both") {
    out <- exon
    xi <- purrr::pmap(list(exon_x = exon$x,
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
    out$p <- c(0, cumsum(purrr::map_int(out$x, length)))
  }

  out$gene_names <- genes
  out$sample_names <- read_tome_sample_names(tome)
  out$dims <- read_tome_data_dims(tome)
  out$dims[2] <- length(genes)

  rhdf5::h5closeAll()

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

  gene_names <- read_tome_vector(tome,"/gene_names")
  sample_names <- read_tome_vector(tome,"/sample_names")

  sample_index <- match(samples, sample_names)

  ## Exon values
  if(regions == "exon" | regions == "both") {

    exon_starts <- read_tome_vector(tome, "data/t_exon/p")[sample_index] + 1
    exon_ends <- read_tome_vector(tome, "data/t_exon/p")[(sample_index + 1)]

    index_list <- unlist(map2(exon_starts, exon_ends, function(x,y) { x:y }))
    split_list <- unlist(map(1:length(exon_starts), function(x) { rep(x, length(exon_starts[x]:exon_ends[x])) }))

    exon_values <- rhdf5::h5read(tome, "data/t_exon/x", index = list(index_list))
    exon_values <- split(exon_values, split_list)
    names(exon_values) <- NULL

    exon_sample_indexes <- rhdf5::h5read(tome, "data/t_exon/i", index = list(index_list))
    exon_sample_indexes <- split(exon_sample_indexes, split_list)
    names(exon_sample_indexes) <- NULL

    exon_pointers <- c(0, cumsum(map_int(exon_values, length)))
    names(exon_pointers) <- NULL


    exon <- list(x = exon_values,
                 i = exon_sample_indexes,
                 p = exon_pointers)


  }

  ## Intron values
  if(regions == "intron" | regions == "both") {

    intron_starts <- read_tome_vector(tome, "data/t_intron/p")[sample_index] + 1
    intron_ends <- read_tome_vector(tome, "data/t_intron/p")[(sample_index + 1)]

    index_list <- unlist(map2(intron_starts, intron_ends, function(x,y) { x:y }))
    split_list <- unlist(map(1:length(intron_starts), function(x) { rep(x, length(intron_starts[x]:intron_ends[x])) }))

    intron_values <- rhdf5::h5read(tome, "data/t_intron/x", index = list(index_list))
    intron_values <- split(intron_values, split_list)
    names(intron_values) <- NULL

    intron_sample_indexes <- rhdf5::h5read(tome, "data/t_intron/i", index = list(index_list))
    intron_sample_indexes <- split(intron_sample_indexes, split_list)
    names(intron_sample_indexes) <- NULL

    intron_pointers <- c(0, cumsum(map_int(intron_values, length)))
    names(intron_pointers) <- NULL


    intron <- list(x = intron_values,
                 i = intron_sample_indexes,
                 p = intron_pointers)

  }

  if(regions == "exon") {
    out <- exon
  } else if(regions == "intron") {
    out <- intron
  } else if(regions == "both") {
    out <- exon
    xi <- purrr::pmap(list(exon_x = exon$x,
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
    out$p <- c(0, cumsum(purrr::map_int(out$x, length)))
  }

  out$gene_names <- read_tome_gene_names(tome)
  out$sample_names <- samples
  out$dims <- read_tome_data_dims(tome, transpose = TRUE)
  out$dims[2] <- length(samples)

  rhdf5::h5closeAll()

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

  Matrix::sparseMatrix(i = unlist(jagged$i),
                       p = jagged$p,
                       x = unlist(jagged$x),
                       index1 = FALSE,
                       dims = jagged$dims,
                       dimnames = list(jagged[[rows]],
                                       jagged[[cols]]))

}

#' Read a whole sparse matrix directly from a tome file
#'
#' @param tome Tome file to read.
#' @param target The sparse matrix to read within the tome file.
#'
read_tome_dgCMatrix <- function(tome,
                                target) {

  root <- rhdf5::H5Fopen(tome)

  i_path <- paste0(target,"/i")
  p_path <- paste0(target,"/p")
  x_path <- paste0(target,"/x")
  dims_path <- paste0(target,"/dims")

  i <- read_tome_vector(root, i_path)
  p <- read_tome_vector(root, p_path)
  x <- read_tome_vector(root, x_path)
  dims <- read_tome_vector(root, dims_path)

  H5Fclose(root)

  h5closeAll()

  Matrix::sparseMatrix(i = i,
                       p = p,
                       x = x,
                       index1 = FALSE,
                       dims = dims)

}

#' Check if an object in a tome file exists
#'
#' @param tome the tome file to check
#' @param name the name of the object in the tome file to check
#'
#' @return logical. TRUE if exists, FALSE if not.
#'
check_tome_existence <- function(tome,
                                 name) {

  ls <- rhdf5::h5ls(tome)
  h5closeAll()
  tome_names <- paste(ls$group, ls$name, sep = "/")
  tome_names <- sub("^//","/",tome_names)

  name %in% tome_names

}
