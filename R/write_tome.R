
#' Save a both exon and intron counts to an HDF5 file (tome)
#'
#' @param exon_mat The exon matrix to store in dgCMatrix format
#' @param intron_mat The intron matrix to store in dgCMatrix format
#' @param tome The HDF5 file to write to
#' @param cols_are Specifies whether columns in the matrix are sample_ids or genes
#' @param overwrite Whether or not to overwrite an existing tome
#' @param compression_level The data compression level for large HDF5 matrix objects. default = 4.
#'
write_tome_data <- function(exon_mat = NULL,
                            intron_mat = NULL,
                            tome = "counts.tome",
                            cols_are = "sample_name",
                            overwrite = F,
                            compression_level = 4) {

  library(Matrix)
  library(rhdf5)

  if(is.null(exon_mat) & is.null(intron_mat)) {
    stop("Provide at least one of exon_mat or intron_mat.")
  }

  H5close()
  if (file.exists(tome) & overwrite) {
    # Delete old version and make a new file
    unlink(tome)
    h5createFile(tome)

  } else if (file.exists(tome) & !overwrite) {
    # Stop if overwrite = F
    stop(paste0(tome," exists. Set overwrite = TRUE to overwrite."))
    return()

  } else if(!file.exists(tome)) {
    # If the file doesn't exist, make a new file
    h5createFile(tome)

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

  root <- H5Fopen(tome, flags = "H5F_ACC_RDWR")
  H5Fclose(root)

  h5createGroup(tome, "data")

  ## Exon Data
  if(!is.null(exon_mat)) {
    h5createGroup(tome, "data/exon")
    h5createGroup(tome, "data/t_exon")
    suppressWarnings({
      # Rows = Samples, Columns = Genes (Fast gene retrieval)
      print("Writing data/exon/x.")
      # data values
      h5createDataset(tome,
                      dataset = "data/exon/x",
                      dims = length(exon_mat@x),
                      chunk = 1000,
                      level = compression_level)
      h5write(exon_mat@x,
              tome,
              "data/exon/x")

      # data indices
      print("Writing data/exon/i.")
      h5createDataset(tome,
                      dataset = "data/exon/i",
                      dims = length(exon_mat@x),
                      chunk = 1000, level = compression_level)
      h5write(exon_mat@i,
              tome,
              "data/exon/i")

      # data index pointers
      print("Writing data/exon/p.")
      h5write(exon_mat@p,
              tome,
              "data/exon/p")

      h5write(c(nrow(exon_mat), ncol(exon_mat)),
              tome,
              "data/exon/dims")

      # Rows = Genes, Columns = Samples (Fast sample retrieval)

      # t_data values
      print("Writing data/t_exon/x.")
      h5createDataset(tome,
                      dataset = "data/t_exon/x",
                      dims = length(t_exon_mat@x),
                      chunk = 1000,
                      level = compression_level)
      h5write(t_exon_mat@x,
              tome,
              "data/t_exon/x")

      # t_data indices
      print("Writing data/t_exon/i.")
      h5createDataset(tome,
                      dataset = "data/t_exon/i",
                      dims = length(t_exon_mat@i),
                      chunk = 1000,
                      level = compression_level)
      h5write(t_exon_mat@i,
              tome,
              "data/t_exon/i")

      # t_data index pointers
      print("Writing data/t_exon/p.")
      h5write(t_exon_mat@p,
              tome,
              "data/t_exon/p")
      h5write(c(nrow(t_exon_mat), ncol(t_exon_mat)),
              tome,
              "data/t_exon/dims")
    })
  }

  if(!is.null(intron_mat)) {
    h5createGroup(tome, "data/intron")
    h5createGroup(tome, "data/t_intron")

    suppressWarnings({

      ## Intron data
      # Rows = Samples, Columns = Genes (Fast gene retrieval)
      print("Writing data/intron/x.")
      # data values
      h5createDataset(tome,
                      dataset = "data/intron/x",
                      dims = length(intron_mat@x),
                      chunk = 1000,
                      level = compression_level)
      h5write(intron_mat@x,
              tome,
              "data/intron/x")

      # data indices
      print("Writing data/intron/i.")
      h5createDataset(tome,
                      dataset = "data/intron/i",
                      dims = length(intron_mat@x),
                      chunk = 1000, level = compression_level)
      h5write(intron_mat@i,
              tome,
              "data/intron/i")

      # data index pointers
      print("Writing data/intron/p.")
      h5write(intron_mat@p,
              tome,
              "data/intron/p")

      h5write(c(nrow(intron_mat), ncol(intron_mat)),
              tome,
              "data/intron/dims")

      # Rows = Genes, Columns = Samples (Fast sample retrieval)

      # t_data values
      print("Writing data/t_intron/x.")
      h5createDataset(tome,
                      dataset = "data/t_intron/x",
                      dims = length(t_intron_mat@x),
                      chunk = 1000,
                      level = compression_level)
      h5write(t_intron_mat@x,
              tome,
              "data/t_intron/x")

      # t_data indices
      print("Writing data/t_intron/i.")
      h5createDataset(tome,
                      dataset = "data/t_intron/i",
                      dims = length(t_intron_mat@i),
                      chunk = 1000,
                      level = compression_level)
      h5write(t_intron_mat@i,
              tome,
              "data/t_intron/i")

      # t_data index pointers
      print("Writing data/t_intron/p.")
      h5write(t_intron_mat@p,
              tome,
              "data/t_intron/p")
      h5write(c(nrow(t_intron_mat), ncol(t_intron_mat)),
              tome,
              "data/t_intron/dims")
    })
  }


  # genes and sample_ids
  print("Writing gene_names.")
  if(!is.null(exon_mat)) {
    h5write(colnames(exon_mat),
            tome,
            "gene_names")
  } else {
    h5write(colnames(intron_mat),
            tome,
            "gene_names")
  }


  print("Writing sample_names.")
  if(!is.null(exon_mat)) {
    h5write(rownames(exon_mat),
            tome,
            "sample_names")
  } else {
    h5write(rownames(intron_mat),
            tome,
            "sample_names")
  }

  if(!is.null(exon_mat)) {
    print("Calculating total exon counts per sample")
    total_exon_counts <- unname(unlist(apply(t_exon_mat, 2, sum)))
    print("Writing data/total_exon_counts.")
    h5write(total_exon_counts,
            tome,
            "data/total_exon_counts")
  }

  if(!is.null(intron_mat)) {
    print("Calculating total intron counts per sample")
    total_intron_counts <- unname(unlist(apply(t_intron_mat, 2, sum)))
    print("Writing data/total_intron_counts.")
    h5write(total_intron_counts,
            tome,
            "data/total_intron_counts")
  }

  H5close()

}

#' Write an annotations table to a tome file.
#'
#' @param anno The annotations data.frame to write. The first column must be "sample_id" or "sample_name".
#' @param tome Path to the target tome file.
#'
write_tome_anno <- function(anno,
                            tome) {

  if(names(anno)[1] == "sample_id") {
    names(anno)[1] <- "sample_name"
  }

  write_tome_data.frame(df = anno,
                        tome = tome,
                        target = "/sample_meta/anno",
                        store_as = "vectors")

}

#' Write annotation desc table to a tome file.
#'
#' @param anno_desc The desc data.frame to write.
#' @param tome Path to the target tome file.
#'
write_tome_anno_desc <- function(anno_desc,
                                 tome) {

  write_tome_data.frame(df = anno_desc,
                        tome = tome,
                        target = "/sample_meta/desc",
                        store_as = "data.frame")

}

#' Write projection coordinates (e.g. tSNE or PCA) to a tome file.
#'
#' @param proj A data.frame with projection coordinates to write.
#' Requires columns: sample_name, x, y (, z optional)
#' @param proj_name The base name of the projection. Should match the projection description table
#' @param tome Path to the target tome file.
#'
write_tome_projection <- function(proj,
                                  proj_name = NULL,
                                  tome) {
  if(!is.null(proj_name)) {
    if(names(proj)[1] == "sample_id") {
      names(proj)[1] <- "sample_name"
    }

    write_tome_data.frame(df = proj,
                          tome = tome,
                          target = paste0("/projection/",proj_name),
                          store_as = "data.frame")
  } else {
    stop("A name for the projection (proj_name) is required.")
  }

}

#' Write a projection descriptions table to a tome file.
#'
#' @param proj_desc The desc data.frame to write.
#' @param tome Path to the target tome file.
#'
write_tome_projection_desc <- function(proj_desc,
                                       tome) {
  write_tome_data.frame(df = proj_desc,
                        tome = tome,
                        target = "/projection/desc",
                        store_as = "data.frame")

}

#' Write a stats table (e.g. median expression per cluster) to a tome file.
#'
#' @param stats The stats data.frame to write.
#' @param tome Path to the target tome file.
#'
write_tome_stats <- function(stats,
                             stats_name = NULL,
                             tome) {

  if(!is.null(stats_name)) {

    write_tome_data.frame(df = stats,
                          tome = tome,
                          target = paste0("/stats/",stats_name),
                          store_as = "vectors")
  } else {
    stop("A name for the stats table (stats_name) is required.")
  }

}

#' Write a stats descriptions table to a tome file.
#'
#' @param stats_desc The desc data.frame to write.
#' @param tome Path to the target tome file.
#'
write_tome_stats_desc <- function(stats_desc,
                                       tome) {

  write_tome_data.frame(df = stats_desc,
                        tome = tome,
                        target = "/stats/desc",
                        store_as = "data.frame")

}

#' Write a dendrogram object as serialized ASCII to a tome file.
#'
#' @param dend The desc data.frame to write.
#' @param dend_name The name of the dendrogram to store.
#' @param tome Path to the target tome file.
#'
write_tome_dend <- function(dend,
                            dend_name,
                            tome) {

  if(!is.null(dend_name)) {

    ls <- h5ls(tome) %>%
      mutate(full_name = ifelse(group == "/",
                                paste0(group, name),
                                paste(group, name, sep = "/")))

    if(!"/dend" %in% ls$full_name) {
      print("Creating group /dend")
      h5createGroup(tome,
                    "/dend")
    }

    dend_target <- paste0("/dend/", dend_name)

    if(dend_target %in% ls$full_name) {
      print(paste0("Removing existing ",dend_target))
      h5_delete(tome, dend_target)
    }

    print(paste0("Writing ",dend_target))

    serial_dend <- rawToChar(serialize(dend, NULL, ascii = T))

    h5write(serial_dend,
            tome,
            dend_target)

  } else {
    stop("A name for the dendrogram (dend_name) is required.")
  }

}

#' Write a dend descriptions table to a tome file.
#'
#' @param dend_desc The desc data.frame to write.
#' @param tome Path to the target tome file.
#'
write_tome_dend_desc <- function(dend_desc,
                                  tome) {

  write_tome_data.frame(df = dend_desc,
                        tome = tome,
                        target = "/dend/desc",
                        store_as = "data.frame")

}

write_tome_exon_lengths <- function(exon_lengths,
                                    tome) {



}

write_tome_intron_lengths <- function(intron_lengths,
                                      tome) {

}
