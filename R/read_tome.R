#' Read Gene Expression Data from a tome file
#'
#' @param tome tome file to read. Required.
#' @param genes A vector of gene names to read. Required.
#' @param regions The gene regions to use. Can be "exon", "intron", or "both". Default = "exon".
#' @param values The type of values to return. Can be "counts" or "cpm". Default = "counts".
#' @param transform Transformation to apply to values. Can be "none", "log", "log2", "log10". Log transforms will add 1 to values before transformation. Default = "none".
#' @param format The format of the output. Can be "data.frame", "matrix", or "dgcMatrix" (sparse matrix). Default is "data.frame".
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
#' @param format The format of the output. Can be "data.frame", "matrix", or "dgcMatrix" (sparse matrix). Default is "data.frame".
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
  H5close()
  gene_names <- h5read(tome,"/gene_names")
  H5close()
  gene_names
}

#' Get all sample names in a tome file
#'
#' @param tome tome file to read.
#'
read_tome_sample_names <- function(tome) {
  H5close()
  sample_names <- h5read(tome,"/sample_names")
  H5close()
  sample_names
}

#' Get exon lengths from a tome file
#'
#' @param tome tome file to read.
#' @param genes specific genes to retrieve
#' @param return_as If "vector", will return only the length values as a vector.
#' If "data.frame", will return a data.frame with two columns: "gene_name" and "exon_length".
#'
read_tome_exon_lengths <- function(tome,
                                   genes = NULL,
                                   return_as = "vector") {
  H5close()
  exon_lengths <- h5read(tome,"/data/exon_lengths")

  if(return_as == "data.frame" | !is.null(genes)) {
    tome_genes <- read_tome_gene_names(tome)
    exon_lengths <- data.frame(gene_name = tome_genes,
                               exon_length = exon_lengths)
  }

  if(!is.null(genes)) {
    exon_lengths <- exon_lengths[gene_name %in% genes,]
  }

  H5close()

  if(return_as == "vector") {
    unlist(exon_lengths$exon_length)
  } else if(return_as == "data.frame") {
    exon_lengths
  }

}

#' Get intron lengths from a tome file
#'
#' @param tome tome file to read.
#' @param genes specific genes to retrieve
#' @param return_as If "vector", will return only the length values as a vector.
#' If "data.frame", will return a data.frame with two columns: "gene_name" and "intron_length".
#'
read_tome_intron_lengths <- function(tome,
                                   genes = NULL,
                                   return_as = "vector") {

  H5close()

  intron_lengths <- h5read(tome,"/data/intron_lengths")

  if(return_as == "data.frame" | !is.null(genes)) {
    tome_genes <- read_tome_gene_names(tome)
    intron_lengths <- data.frame(gene_name = tome_genes,
                               intron_length = intron_lengths)
  }

  if(!is.null(genes)) {
    intron_lengths <- intron_lengths[gene_name %in% genes,]
  }

  H5close()

  if(return_as == "vector") {
    unlist(intron_lengths$intron_length)
  } else if(return_as == "data.frame") {
    intron_lengths
  }

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
  H5close()
  if(region == "exon") {

    total_counts <- h5read(tome,"data/total_exon_counts")

  } else if(region == "intron") {

    total_counts <- h5read(tome,"data/total_intron_counts")

  } else if(region == "both") {

    total_exon_counts <- h5read(tome,"data/total_exon_counts")
    total_intron_counts <- h5read(tome,"data/total_intron_counts")
    total_counts <- total_exon_counts + total_intron_counts

  }

  H5close()
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
  H5close()

  if(transpose == F) {
    dims <- h5read(tome,"data/exon/dims")
  } else if(transpose == T) {
    dims <- h5read(tome,"data/t_exon/dims")
  }

  H5close()

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
                                 match_type = "grep",
                                 get_all = FALSE)
  } else {
    anno <- read_tome_data.frame(tome,
                                 "/sample_meta/anno",
                                 stored_as = "vectors",
                                 columns = "sample_name",
                                 get_all = TRUE)
  }

  anno <- anno %>%
    select(sample_name, everything())

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

#' Read gene metadata table from a tome file
#'
#' @param tome The location of the tome file to read.
#' @param columns Specific columns to read If NULL, will get all columns. Default is NULL.
#'
read_tome_gene_meta <- function(tome,
                                columns = NULL) {

  library(rhdf5)
  library(purrr)

  if(!is.null(columns)) {
    genes <- read_tome_data.frame(tome,
                                  "/gene_meta/genes",
                                  stored_as = "vectors",
                                  columns = c("gene_name", columns),
                                  get_all = FALSE)
  } else {
    genes <- read_tome_data.frame(tome,
                                  "/gene_meta/genes",
                                  stored_as = "vectors",
                                  columns = "gene_name",
                                  get_all = TRUE)
  }

  genes

}

#' Read gene metadata descriptions table from a tome file
#'
#' @param tome The location of the tome file to read.
#'
read_tome_gene_meta_desc <- function(tome) {

  library(rhdf5)
  library(purrr)

  H5close()

  desc <- read_tome_data.frame(tome,
                               "/gene_meta/desc",
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
                            stats_name = NULL,
                            columns = NULL,
                            get_all = FALSE) {
  library(rhdf5)
  library(purrr)

  H5close()

  ls <- h5ls(tome)
  stats_names <- ls$name[ls$group == "/stats"]
  stats_names <- stats_names[stats_names != "desc"]

  if(is.null(stats_name)) { stats_name <- ".namenotfound"}

  if(stats_name %in% stats_names) {

    stats_target <- paste0("/stats/", stats_name)

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
                                    columns = columns,
                                    get_all = TRUE)
    }
  } else {
    H5close()

    if(length(stats_names) > 0) {
      stats_message <- paste0("A stats table name (stats_name) is required. Available in this dataset are: ", paste(stats_names,collapse = ", "))
    } else {
      stats_message <- "No stats tables are found in this dataset."
    }

    stop(stats_message)

  }
  H5close()


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

#' Read projection coordinates table from a tome file
#'
#' @param tome the location of the tome file to read.
#' @param proj_name the projection to read
#'
read_tome_projection <- function(tome,
                                 proj_name = NULL) {

  library(rhdf5)
  H5close()

  ls <- h5ls(tome)
  proj_names <- ls$name[ls$group == "/projection"]
  proj_names <- proj_names[proj_names != "desc"]

  if(is.null(proj_name)) { proj_name <- ".namenotfound" }

  if(proj_name %in% proj_names) {
    proj_target <- paste0("projection/",proj_name)

    proj <- read_tome_data.frame(tome,
                                 proj_target,
                                 stored_as = "data.frame")
    H5close()

    proj
  } else {
    H5close()

    if(length(proj_names) > 0) {
      proj_message <- paste0("A projection name (proj_name) is required. Available in this dataset are: ", paste(proj_names,collapse = ", "))
    } else {
      proj_message <- "No projections are found in this dataset."
    }
    stop(proj_message)
  }

}

#' Read projection descriptions table from a tome file
#'
#' @param tome The location of the tome file to read.
#'
read_tome_projection_desc <- function(tome) {
  library(rhdf5)

  desc <- read_tome-data.frame(tome,
                               "/projection/desc",
                               stored_as = "data.frame")

  desc
}

#' Read dendrogram object from a tome file
#'
#' @param tome the location of the tome file to read.
#' @param dend_name the dendrogram to read.
#'
read_tome_dend <- function(tome,
                           dend_name = NULL) {
  library(rhdf5)
  H5close()

  ls <- h5ls(tome)
  dend_names <- ls$name[ls$group == "/dend"]
  dend_names <- dend_names[dend_names != "desc"]

  if(is.null(dend_name)) { dend_name <- ".namenotfound" }

  if(dend_name %in% dend_names) {
    dend_target <- paste0("dend/",dend_name)

    dend <- read_tome_serialized(tome,
                                 dend_target)
    H5close()

    dend
  } else {
    H5close()

    if(length(dend_names) > 0) {
      dend_message <- paste0("A dendrogram name (dend_name) is required. Available in this dataset are: ", paste(dend_names,collapse = ", "))
    } else {
      dend_message <- "No dendrograms are found in this dataset."
    }
    stop(dend_message)
  }

}


#' Read mapping table from a tome file
#'
#' @param tome the location of the tome file to read.
#' @param mapping_name the name of the mapping table to read
#' @param columns selected columns to read. If NULL, reads all columns. Default = NULL.
#' @param get_all logical, whether or not to append all other columns after the specified columns. Default = FALSE.
#'
read_tome_mapping <- function(tome,
                              mapping_name = NULL,
                              columns = NULL,
                              get_all = FALSE) {
  library(rhdf5)
  library(purrr)
  H5close()

  ls <- h5ls(tome)
  mapping_names <- ls$name[ls$group == "/mapping"]
  mapping_names <- mapping_names[mapping_names != "desc"]

  if(is.null(mapping_name)) { mapping_name <- ".namenotfound"}

  if(mapping_name %in% mapping_names) {

    mapping_target <- paste0("mapping/", mapping_name)

    if(!is.null(columns)) {
      mapping <- read_tome_data.frame(tome,
                                      mapping_target,
                                      stored_as = "vectors",
                                      columns = c("sample_name", columns),
                                      get_all = get_all)
    } else {
      mapping <- read_tome_data.frame(tome,
                                      mapping_target,
                                      stored_as = "vectors",
                                      columns = "sample_name",
                                      get_all = get_all)
    }
  } else {
    H5close()

    if(length(mapping_names) > 0) {
      mapping_message <- paste0("A mapping table name (mapping_name) is required. Available in this dataset are: ", paste(mapping_names,collapse = ", "))
    } else {
      mapping_message <- "No mapping tables are found in this dataset."
    }

    stop(mapping_message)

  }
  H5close()

  mapping
}

#' Read mapping descriptions table from a tome file
#'
#' @param tome The location of the tome file to read.
#'
read_tome_mapping_desc <- function(tome) {

  library(rhdf5)
  library(purrr)

  H5close()

  desc <- read_tome_data.frame(tome,
                               "/mapping/desc",
                               stored_as = "data.frame")

  desc

}
