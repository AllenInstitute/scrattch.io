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
    root <- H5Fopen(tome)

    if(regions == "exon") {
      total_counts <- c(h5read(root, "/data/total_exon_counts"))
    } else if(regions == "intron") {
      total_counts <- c(h5read(root, "/data/total_intron_counts"))
    } else if(regions == "both") {
      total_counts <- c(h5read(root, "/data/total_exon_counts") + h5read(root, "/data/total_intron_counts"))
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
      total_counts <- c(h5read(root, "/data/total_exon_counts")[sample_index])
    } else if(regions == "intron") {
      total_counts <- c(h5read(root, "/data/total_intron_counts")[sample_index])
    } else if(regions == "both") {
      total_counts <- c(h5read(root, "/data/total_exon_counts")[sample_index] + h5read(root, "/data/total_intron_counts")[sample_index])
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
#' @param tome tome file to read.
#'
read_tome_gene_names <- function(tome) {
  root <- H5Fopen(tome)
  gene_names <- h5read(root,"/gene_names")
  H5Fclose(root)
  gene_names
}

#' Get all sample names in a tome file
#'
#' @param tome tome file to read.
#'
read_tome_sample_names <- function(tome) {
  root <- H5Fopen(tome)
  sample_names <- h5read(root,"/sample_names")
  H5Fclose(root)
  sample_names
}

#' Get total per-gene counts from a tome-file.
#'
#' Useful for computing CPM or FPKM
#'
#' @param tome tome file to read.
#' @param region the gene regions to use. Can be "exon", "intron", or "both". Default = "exon"
#'
read_tome_total_counts <- function(tome,
                                   region = "exon") {
  root <- H5Fopen(tome)
  if(region == "exon") {

    total_counts <- h5read(root,"data/total_exon_counts")

  } else if(region == "intron") {

    total_counts <- h5read(root,"data/total_intron_counts")

  } else if(region == "both") {

    total_exon_counts <- h5read(root,"data/total_exon_counts")
    total_intron_counts <- h5read(root,"data/total_intron_counts")
    total_counts <- total_exon_counts + total_intron_counts

  }

  H5Fclose(root)
  total_counts
}

#' Get dims for data stored in a tome file.
#'
#' @param tome tome file to read.
#' @param transpose logical indicating which orientation to get dims. Default = FALSE.
#' If FALSE, gives dimensions with samples as rows and genes as columns.
#' If TRUE,  gives dimensions for genes as rows and samples as columns.
#'
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

  if(!is.null(groups)) {
    anno <- read_tome_data.frame(tome,
                                 "/sample_meta/anno",
                                 stored_as = "vectors",
                                 columns = c("sample_name", groups),
                                 get_all = FALSE)
  } else {
    anno <- read_tome_data.frame(tome,
                                 "/sample_meta/anno",
                                 stored_as = "vectors",
                                 columns = "sample_name",
                                 get_all = TRUE)
  }

  anno

}

#' Read annotation descriptions table from a tome file
#'
#' @param tome The location of the tome file to read.
#'
read_tome_anno_desc <- function(tome) {

  library(rhdf5)
  library(purrr)

  H5close()

  desc <- read_tome_data.frame(tome,
                               "/sample_meta/desc",
                               stored_as = "data.frame")

  desc

}

#' Read dendrogram descriptions table from a tome file
#'
#' @param tome The location of the tome file to read.
#'
read_tome_dend_desc <- function(tome) {

  library(rhdf5)
  library(purrr)

  H5close()

  desc <- read_tome_data.frame(tome,
                               "/dend/desc",
                               stored_as = "data.frame")

  desc

}

#' Read stats table from a tome file
#'
#' @param tome the location of the tome file to read.
#' @param stats_name the name of the stats table to read
#' @param columns selected columns to read. If NULL, reads all columns. Default = NULL.
#' @param get_all logical, whether or not to append all other columns after the specified columns. Default = FALSE.
#'
read_tome_stats <- function(tome,
                            stats_name,
                            columns = NULL,
                            get_all = FALSE) {
  library(rhdf5)
  library(purrr)

  stats_target <- paste0("stats/", stats_name)

  if(!is.null(columns)) {
    stats <- read_tome_data.frame(tome,
                                  stats_target,
                                  stored_as = "vectors",
                                  columns = c("sample_name", columns),
                                  get_all = get_all)
  } else {
    stats <- read_tome_data.frame(tome,
                                  stats_target,
                                  stored_as = "vectors",
                                  columns = "sample_name",
                                  get_all = get_all)
  }

  stats
}

#' Read stats descriptions table from a tome file
#'
#' @param tome The location of the tome file to read.
#'
read_tome_stats_desc <- function(tome) {

  library(rhdf5)
  library(purrr)

  H5close()

  desc <- read_tome_data.frame(tome,
                               "/stats/desc",
                               stored_as = "data.frame")

  desc

}
