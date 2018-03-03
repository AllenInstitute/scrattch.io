
#' Read Gene Expression Data from a tome file
#'
#' @param tome tome file to read. Required.
#' @param genes A vector of gene names to read. Required.
#' @param regions The gene regions to use. Can be "exon", "intron", or "both". Default = "exon".
#' @param values The type of values to return. Can be "counts" or "cpm". Default = "counts".
#' @param transform Transformation to apply to values. Can be "none", "log", "log2", "log10". Log transforms will add 1 to values before transformation. Default = "none".
#' @param form The format of the output. Can be "data.frame", or "matrix". Default is "data.frame".
#'
#' @return A data.frame with sample_name as the first column and each subsequent column
#' containing gene expression values and named for the genes; Or a matrix with columns as genes and rows as samples.
#'
read_tome_gene_data <- function(tome,
                                genes,
                                regions = "exon",
                                type = "counts",
                                transform = "none",
                                format = "data.frame") {
  library(rhdf5)
  library(purrr)
  library(dplyr)
  library(Matrix)

  jagged <- read_tome_genes_jagged(tome,
                                     genes,
                                     regions)

  if(format == "data.frame") {
    out <- jagged_to_data.frame(jagged,
                                rows = "sample_names",
                                cols = "gene_names")
  } else if(format == "matrix") {
    out <- jagged_to_matrix(jagged,
                            rows = "sample_names",
                            cols = "gene_names")
  } else if(format == "dgCMatrix") {
    out <- jagged_to_dgCMatrix(jagged,
                               rows = "sample_names",
                               cols = "gene_names")
  }

  if(type == "cpm") {

    if(regions == "exon") {
      total_counts <- h5read(root, "/data/total_exon_counts")
    } else if(regions == "intron") {
      total_counts <- h5read(root, "/data/total_intron_counts")
    } else if(regions == "both") {
      total_counts <- h5read(root, "/data/total_exon_counts") + h5read(root, "/data/total_intron_counts")
    }
    out[,genes] <- out[,genes]/(total_counts/1e6)
  }

  if(transform == "log") {
    out[,genes] <- log(out[,genes] + 1)
  } else if(transform == "log2") {
    out[,genes] <- log2(out[,genes] + 1)
  } else if(transform == "log10") {
    out[,genes] <- log10(out[,genes] + 1)
  }

  out

}

#' Read Sample Expression Data from a tome file
#'
#' @param tome tome file to read. Required.
#' @param samples A vector of sample names to read. Required
#' @param regions The gene regions to use. Can be "exon", "intron", or "both". Default = "exon".
#' @param values The type of values to return. Can be "counts" or "cpm". Default = "counts".
#' @param transform Transformation to apply to values. Can be "none", "log", "log2", "log10". Log transforms will add 1 to values before transformation. Default = "none".
#' @param format The format of the output. Can be "data.frame", "matrix", or "dgCMatrix". Default is "data.frame".
#'
#' @return A data.frame with gene_name as the first column and each subsequent column
#' containing gene expression values and named for the samples; Or a matrix with columns as samples and rows as genes.
#'
read_tome_sample_data <- function(tome,
                                  samples,
                                  regions = "exon",
                                  type = "counts",
                                  transform = "none",
                                  format = "data.frame") {
  library(rhdf5)
  library(purrr)
  library(dplyr)

  jagged <- read_tome_samples_jagged(tome,
                                     samples,
                                     regions)

  if(format == "data.frame") {
    out <- jagged_to_data.frame(jagged,
                                rows = "gene_names",
                                cols = "sample_names")
  } else if(format == "matrix") {
    out <- jagged_to_matrix(jagged,
                            rows = "gene_names",
                            cols = "sample_names")
  } else if(format == "dgCMatrix") {
    out <- jagged_to_dgCMatrix(jagged,
                               rows = "gene_names",
                               cols = "sample_names")
  }

  if(type == "cpm") {
    root <- H5Fopen(tome)
    sample_names <- h5read(root,"/sample_names")
    sample_index <- match(samples, sample_names)

    if(regions == "exon") {
      total_counts <- h5read(root, "/data/total_exon_counts")[sample_index]
    } else if(regions == "intron") {
      total_counts <- h5read(root, "/data/total_intron_counts")[sample_index]
    } else if(regions == "both") {
      total_counts <- h5read(root, "/data/total_exon_counts")[sample_index] + h5read(root, "/data/total_intron_counts")[sample_index]
    }
    out[,samples] <- sweep(out[,samples], 2, (total_counts/1e6), "/")
  }

  if(transform == "log") {
    out[,samples] <- log(out[,samples] + 1)
  } else if(transform == "log2") {
    out[,samples] <- log2(out[,samples] + 1)
  } else if(transform == "log10") {
    out[,samples] <- log10(out[,samples] + 1)
  }

  out

}

#' Get all gene names in a tome file
#'
#' @param tome Tome file to read.
#'
read_tome_gene_names <- function(tome) {
  root <- H5Fopen(tome)
  gene_names <- h5read(root,"/gene_names")
  H5Fclose(root)
  gene_names
}

read_tome_sample_names <- function(tome) {
  root <- H5Fopen(tome)
  sample_names <- h5read(root,"/sample_names")
  H5Fclose(root)
  sample_names
}

read_tome_total_exon_counts <- function(tome) {
  root <- H5Fopen(tome)
  total_exon_counts <- h5read(root,"data/total_exon_counts")
  H5Fclose(root)
  total_exon_counts
}

read_tome_total_intron_counts <- function(tome) {
  root <- H5Fopen(tome)
  total_intron_counts <- h5read(root,"data/total_intron_counts")
  H5Fclose(root)
  total_intron_counts
}

read_tome_total_both_counts <- function(tome) {
  root <- H5Fopen(tome)
  total_exon_counts <- h5read(root,"data/total_exon_counts")
  total_intron_counts <- h5read(root,"data/total_intron_counts")
  H5Fclose(root)
  total_both_counts <- total_exon_counts + total_intron_counts
}

read_tome_data_dims <- function(tome,
                                transpose = F) {
  root <- H5Fopen(tome)

  if(transpose == F) {
    dims <- h5read(root,"data/exon/dims")
  } else if(transpose == T) {
    dims <- h5read(root,"data/t_exon/dims")
  }

  H5Fclose(root)

  dims
}

#' Read annotations table from a tome file
#'
#' @param tome The location of the tome file to read.
#' @param groups The groups to read - matches column names using grep. Can provide multiple with c(). If NULL, will get all columns. Default is NULL.
#'
read_tome_anno <- function(tome,
                           groups = NULL) {

  library(rhdf5)
  library(purrr)

  H5close()

  ls <- h5ls(tome)

  annos <- ls$name[ls$group == "/sample_meta/anno"]

  if(!is.null(groups)) {
    if(length(groups) > 1) {
      groups <- paste(groups, collapse = "|")
    }
    annos <- c("sample_name", annos[grepl(groups, annos)])
  }

  annos <- c("sample_name",annos[annos != "sample_name"])

  anno <- map(annos, function(x) h5read(tome,
                                        paste0("/sample_meta/anno/",x)))
  names(anno) <- annos
  anno <- as.data.frame(anno)
  anno

}


read_tome_genes_jagged <- function(tome,
                                  genes,
                                  regions = "exon") {
  library(rhdf5)
  library(purrr)
  library(dplyr)
  library(Matrix)

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

  out

}

read_tome_samples_jagged <- function(tome,
                                     samples,
                                     regions = "exon",
                                     type = "counts",
                                     transform = "none") {
  library(rhdf5)
  library(purrr)
  library(dplyr)
  library(Matrix)

  root <- H5Fopen(tome)
  gene_names <- h5read(root,"/gene_names")
  sample_names <- h5read(root,"/sample_names")

  gene_index <- match(genes, gene_names)

  ## Exon values
  if(regions == "exon" | regions == "both") {

    exon_starts <- h5read(root, "data/t_exon/p")[gene_index] + 1
    exon_ends <- h5read(root, "data/t_exon/p")[(gene_index + 1)]

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

    intron_starts <- h5read(root, "t_data/intron/p")[gene_index] + 1
    intron_ends <- h5read(root, "data/t_intron/p")[(gene_index + 1)]

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

  out$gene_names <- genes
  out$sample_names <- read_tome_sample_names(tome)
  out$dims <- read_tome_data_dims(tome, transpose = true)
  out$dims[2] <- length(samples)

  out

}

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
