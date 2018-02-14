
#' Read Gene Expression Data from a tome file
#'
#' @param tome_file tome file to read. Required.
#' @param genes A vector of gene names to read. Required.
#' @param regions The gene regions to use. Can be "exon", "intron", or "both". Default = "exon".
#' @param values The type of values to return. Can be "counts" or "cpm". Default = "counts".
#' @param transform Transformation to apply to values. Can be "none", "log", "log2", "log10". Log transforms will add 1 to values before transformation. Default = "none".
#' @param form The format of the output. Can be "data.frame", or "matrix". Default is "data.frame".
#'
#' @return A data.frame with sample_name as the first column and each subsequent column
#' containing gene expression values and named for the genes; Or a matrix with columns as genes and rows as samples.
#'
read_tome_gene_data <- function(tome_file,
                                 genes,
                                 regions = "exon",
                                 type = "counts",
                                 transform = "none",
                                 form = "data.frame") {
  library(rhdf5)
  library(purrr)
  library(dplyr)
  library(Matrix)

  root <- H5Fopen(tome_file)
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

    if(form == "data.frame") {
      exon <- data.frame(sample_name = sample_names)

      exons <- pmap(list(x = exon_sample_names,
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
        exon <- reduce(c(list(exon), exons),
                       left_join,
                       by = "sample_name")
      })

      if(type == "cpm") {
        total_exon_counts <- h5read(root, "/total_exon_counts")
        exon[,genes] <- exon[,genes]/(total_exon_counts/1e6)
      }

      exon[is.na(exon)] <- 0
    } else if(form == "matrix") {

      exon <- matrix(0, nrow = length(sample_names), ncol = length(genes))
      rownames(exon) <- sample_names
      colnames(exon) <- genes
      for(i in 1:length(genes)) {
        exon[exon_sample_names[[i]], i] <- exon_values[[i]]
      }

      if(type == "cpm") {
        total_exon_counts <- h5read(root, "/total_exon_counts")
        exon <- apply(exon, 2, function(x) x/(total_exon_counts/1e6))
        rownames(exon) <- sample_names
      }

    }



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

    if(form == "data.frame") {
      intron <- data.frame(sample_name = sample_names)

      introns <- pmap(list(x = intron_sample_names,
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
        intron <- reduce(c(list(intron), introns),
                         left_join,
                         by = "sample_name")
      })

      if(type == "cpm") {
        total_intron_counts <- h5read(root, "/total_intron_counts")
        intron[,genes] <- intron[,genes]/(total_intron_counts/1e6)
      }

      intron[is.na(intron)] <- 0

    } else if(form == "matrix") {

      intron <- matrix(0, nrow = length(sample_names), ncol = length(genes))
      rownames(intron) <- sample_names
      colnames(intron) <- genes
      for(i in 1:length(genes)) {
        intron[intron_sample_names[[i]], i] <- intron_values[[i]]
      }

      if(type == "cpm") {
        total_intron_counts <- h5read(root, "/total_intron_counts")
        intron <- apply(intron, 2, function(x) x/(total_intron_counts/1e6))
        rownames(intron) <- sample_names
      }

    }

  }

  H5Fclose(root)

  if(regions == "exon") {
    out <- exon
  } else if(regions == "intron") {
    out <- intron
  } else if(regions == "both") {
    out <- exon
    out[,genes] <- out[,genes] + intron[,genes]
    out
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
#' @param tome_file tome file to read. Required.
#' @param samples A vector of sample names to read. Required
#' @param regions The gene regions to use. Can be "exon", "intron", or "both". Default = "exon".
#' @param values The type of values to return. Can be "counts" or "cpm". Default = "counts".
#' @param transform Transformation to apply to values. Can be "none", "log", "log2", "log10". Log transforms will add 1 to values before transformation. Default = "none".
#' @param form The format of the output. Can be "data.frame", or "matrix". Default is "data.frame".
#'
#' @return A data.frame with gene_name as the first column and each subsequent column
#' containing gene expression values and named for the samples; Or a matrix with columns as samples and rows as genes.
#'
read_tome_sample_data <- function(tome_file,
                                   samples,
                                   regions = "exon",
                                   type = "counts",
                                   transform = "none",
                                   form = "data.frame") {
  library(rhdf5)
  library(purrr)
  library(dplyr)

  root <- H5Fopen(tome_file)
  sample_names <- h5read(root,"/sample_names")
  gene_names <- h5read(root,"/gene_names")

  if(is.null(samples)) {
    samples <- sample_names
  }

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

    if(form == "data.frame") {

      exon <- data.frame(gene_id = gene_names)

      exons <- pmap(list(x = exon_gene_names,
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
        exon <- reduce(c(list(exon), exons),
                       left_join,
                       by = "gene_id")
      })

      if(type == "cpm") {
        total_exon_counts <- h5read(root, "/total_exon_counts")[sample_index]
        exon[,samples] <- sweep(exon[,samples], 2, (total_exon_counts/1e6), "/")
      }

      exon[is.na(exon)] <- 0

    } else if(form == "matrix") {

      exon <- matrix(0, nrow = length(gene_names), ncol = length(samples))
      rownames(exon) <- gene_names
      colnames(exon) <- samples
      for(i in 1:length(genes)) {
        exon[exon_gene_names[[i]], i] <- exon_values[[i]]
      }

      if(type == "cpm") {
        total_exon_counts <- h5read(root, "/total_exon_counts")[sample_index]
        exon <- sweep(exon, 2, (total_exon_counts/1e6), "/")

      }

    }


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

    if(form == "data.frame") {

      intron <- data.frame(gene_id = gene_names)

      introns <- pmap(list(x = intron_gene_names,
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
        intron <- reduce(c(list(intron), introns),
                         left_join,
                         by = "gene_id")
      })

      if(type == "cpm") {
        total_intron_counts <- h5read(root, "/total_intron_counts")[sample_index]
        intron[,samples] <- sweep(intron[,samples], 2, (total_intron_counts/1e6), "/")
      }

      intron[is.na(intron)] <- 0

    } else if(form == "matrix") {

      intron <- matrix(0, nrow = length(gene_names), ncol = length(samples))
      rownames(intron) <- gene_names
      colnames(intron) <- samples
      for(i in 1:length(genes)) {
        intron[intron_gene_names[[i]], i] <- intron_values[[i]]
      }

      if(type == "cpm") {
        total_intron_counts <- h5read(root, "/total_intron_counts")[sample_index]
        intron <- sweep(intron, 2, (total_intron_counts/1e6), "/")

      }

    }

  }

  H5Fclose(root)

  if(regions == "exon") {
    out <- exon
  } else if(regions == "intron") {
    out <- intron
  } else if(regions == "both") {
    out <- exon
    out[,samples] <- out[,samples] + intron[,samples]
    out
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
#' @param tome file.
#'
read_tome_gene_names <- function(tome_file) {
  root <- H5Fopen(tome_file)
  gene_names <- h5read(root,"/gene_names")
  H5Fclose(root)
  gene_names
}

read_tome_sample_names <- function(tome_file) {
  root <- H5Fopen(tome_file)
  sample_names <- h5read(root,"/sample_names")
  H5Fclose(root)
  sample_names
}

read_tome_total_exon_counts <- function(tome_file) {
  root <- H5Fopen(tome_file)
  total_exon_counts <- h5read(root,"/total_exon_counts")
  H5Fclose(root)
  total_exon_counts
}

read_tome_total_intron_counts <- function(tome_file) {
  root <- H5Fopen(tome_file)
  total_intron_counts <- h5read(root,"/total_intron_counts")
  H5Fclose(root)
  total_intron_counts
}

read_tome_total_both_counts <- function(tome_file) {
  root <- H5Fopen(tome_file)
  total_exon_counts <- h5read(root,"/total_exon_counts")
  total_intron_counts <- h5read(root,"/total_intron_counts")
  H5Fclose(root)
  total_both_counts <- total_exon_counts + total_intron_counts
}

read_tome_data_dims <- function(tome_file,
                                 transpose = F) {
  root <- H5Fopen(tome_file)

  if(transpose == F) {
    dims <- h5read(root,"/exon/dims")
  } else if(transpose == T) {
    dims <- h5read(root,"/t_exon/dims")
  }

  H5Fclose(root)

  dims
}

#' Read annotations table from a tome file
#'
#' @param tome_file The location of the tome file to read.
#' @param groups The groups to read - matches column names using grep. Can provide multiple with c(). If NULL, will get all columns. Default is NULL.
#'
read_tome_anno <- function(tome_file,
                            groups = NULL) {

  library(rhdf5)
  library(purrr)

  H5close()

  ls <- h5ls(tome_file)

  annos <- ls$name[ls$group == "/sample_meta/anno"]

  if(!is.null(groups)) {
    if(length(groups) > 1) {
      groups <- paste(groups, collapse = "|")
    }
    annos <- c("sample_name", annos[grepl(groups, annos)])
  }

  annos <- c("sample_name",annos[annos != "sample_name"])

  anno <- map(annos, function(x) h5read(tome_file,
                                        paste0("/sample_meta/anno/",x)))
  names(anno) <- annos
  anno <- as.data.frame(anno)
  anno

}
