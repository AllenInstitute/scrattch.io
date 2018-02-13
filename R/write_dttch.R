
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
                             cols_are = "sample_name",
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
              "exon/dims")

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
              "t_exon/dims")
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
              "intron/dims")

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
              "t_intron/dims")
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

#' Write an annotations table to a dttch file.
#'
#' @param anno The annotations data.frame to write. The first column must be "sample_id" or "sample_name".
#' @param out_file The dttch file to write to
#' @param mode The write mode. "replace" will remove existing annotations and replace it with anno. "add" will only add columns not already found. Default = "add".
#'
write_dttch_anno <- function(anno,
                             out_file,
                             mode = "add") {

  library(rhdf5)
  library(purrr)

  H5close()

  if(names(anno)[1] == "sample_id") {
    names(anno)[1] <- "sample_name"
  }

  if(!file.exists(out_file)) {
    print(paste0(out_file," doesn't exist. Creating new file."))
    h5createFile(out_file)
    H5close()
  }

  ls <- h5ls(out_file) %>%
    mutate(full_name = paste0(group, name))

  # check for sample_meta
  if(!"/sample_meta" %in% ls$group) {
    print("Creating Group /sample_meta")
    h5createGroup(out_file, "/sample_meta")
  }

  # check for anno
  if(!"/sample_meta/anno" %in% ls$group) {
    print("Creating Group /sample_meta/anno")
    h5createGroup(out_file, "/sample_meta/anno")
  }

  existing_anno <- ls %>% filter(group == "/sample_meta/anno")
  anno_cols <- names(anno)

  # If overwrite, remove everything
  if(mode == "replace") {
    print("Removing existing /sample_meta/anno")

    existing_anno <- ls %>% filter(group == "/sample_meta/anno")
    print(nrow(existing_anno))
    if(nrow(existing_anno) > 0) {
      walk(existing_anno$full_name,
           function(x) {
             h5_delete(out_file, x)
           }
      )
    }

    print("Writing /sample_meta/anno")
    walk(anno_cols, function(x) {
      target <- paste0("/sample_meta/anno/",x)
      h5write(anno[[x]],
              out_file,
              target)
    })

  } else if(mode == "add") {

    if(nrow(existing_anno) > 0) {
      anno_cols <- anno_cols[!anno_cols %in% existing_anno$name]

      if(length(anno_cols) > 0) {
        print(paste0("Adding columns ", paste(anno_cols, collapse = ", "), " to /sample_meta/anno."))

        # Order to match existing annotations
        existing_sample_names <- h5read(out_file, "/sample_meta/anno/sample_name")

        anno <- anno %>%
          select(one_of("sample_name", anno_cols))

        anno <- anno[match(anno$sample_name, existing_sample_names),]

        # Write the new columns
        walk(anno_cols, function(x) {
          target <- paste0("/sample_meta/anno/",x)
          h5write(anno[[x]],
                  out_file,
                  target)
        })
      } else {
        print("No new columns to add. If you want to replace values, try mode = 'replace' to replace the whole anno table.")
      }
    } else {
      print("Writing /sample_meta/anno")
      anno_cols <- names(anno)
      walk(anno_cols, function(x) {
        target <- paste0("/sample_meta/anno/",x)
        h5write(anno[[x]],
                out_file,
                target)
      })
    }
  }

}
