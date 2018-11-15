#' Generate colors and ids for numeric annotations
#'
#' @param df data frame to annotate
#' @param col name of the numeric column to annotate
#' @param base base name for the annotation, which wil be used in the desc table. If not provided, will use col as base.
#' @param scale The scale to use for assigning colors. Options are "linear","log10","log2, and "zscore"
#' @param na_val The value to use to replace NAs. default = 0.
#' @param colorset A vector of colors to use for the color gradient. default = c("darkblue","white","red")
#'
#' @return A modified data frame: the annotated column will be renamed base_label, and base_id and base_color columns will be appended
#'
annotate_num <- function (df,
                          col = NULL, base = NULL,
                          scale = "log10", na_val = 0,
                          colorset = c("darkblue", "white", "red")) {

  #library(lazyeval)
  #library(dplyr)

  if(class(try(is.character(col), silent = T)) == "try-error") {
    col <- lazyeval::expr_text(col)
  } else if(class(col) == "NULL") {
    stop("Specify a column (col) to annotate.")
  }

  if(class(try(is.character(base), silent = T)) == "try-error") {
    base <- lazyeval::expr_text(base)
  } else if(class(base) == "NULL") {
    base <- col
  }

  if (!is.numeric(df[[col]])) {
    df[[col]] <- as.numeric(df[[col]])
  }

  df[[col]][is.na(df[[col]])] <- na_val

  x <- df[[col]]

  annotations <- data.frame(label = unique(x)) %>%
    dplyr::arrange(label) %>%
    dplyr::mutate(id = 1:dplyr::n())

  if (scale == "log10") {
    colors <- values_to_colors(log10(annotations$label + 1), colorset = colorset)
  } else if(scale == "log2") {
    colors <- values_to_colors(log2(annotations$label + 1), colorset = colorset)
  } else if(scale == "zscore") {
    colors <- values_to_colors(scale(annotations$label), colorset = colorset)
  } else if(scale == "linear") {
    colors <- values_to_colors(annotations$label, colorset = colorset)
  }
  annotations <- mutate(annotations, color = colors)
  names(annotations) <- paste0(base, c("_label", "_id", "_color"))

  names(df)[names(df) == col] <- paste0(base,"_label")
  df <- dplyr::left_join(df, annotations, by = paste0(base,"_label"))
  df
}

#' Generate colors and ids for categorical annotations
#'
#' @param df data frame to annotate
#' @param col name of the character column to annotate
#' @param base base name for the annotation, which wil be used in the desc
#'   table. If not provided, will use col as base.
#' @param sort_label a logical value to determine if the data in col should be
#'   arranged alphanumerically before ids are assigned. default = T.
#' @param na_val The value to use to replace NAs. default = "ZZ_Missing".
#' @param colorset The colorset to use for assigning category colors. Options
#'   are "rainbow","viridis","inferno","magma", and "terrain"
#' @param color_order The order in which colors should be assigned. Options are
#'   "sort" and "random". "sort" assigns colors in order; "random" will randomly
#'   assign colors.
#'
#' @return A modified data frame: the annotated column will be renamed
#'   base_label, and base_id and base_color columns will be appended
#'
annotate_cat <- function(df,
                         col = NULL, base = NULL,
                         sort_label = T, na_val = "ZZ_Missing",
                         colorset = "varibow", color_order = "sort") {

  #library(dplyr)
  #library(lazyeval)
  #library(viridisLite)

  if(class(try(is.character(col), silent = T)) == "try-error") {
    col <- lazyeval::expr_text(col)
  } else if(class(col) == "NULL") {
    stop("Specify a column (col) to annotate.")
  }

  if(class(try(is.character(base), silent = T)) == "try-error") {
    base <- lazyeval::expr_text(base)
  } else if(class(base) == "NULL") {
    base <- col
  }

  if(!is.character(df[[col]])) {
    df[[col]] <- as.character(df[[col]])
  }

  df[[col]][is.na(df[[col]])] <- na_val

  x <- df[[col]]

  annotations <- data.frame(label = unique(x), stringsAsFactors = F)

  if(sort_label) {
    annotations <- annotations %>% dplyr::arrange(label)
  }

  annotations <- annotations %>%
    dplyr::mutate(id = 1:n())

  if(colorset == "varibow") {
    colors <- varibow(nrow(annotations))
  } else if(colorset == "rainbow") {
    colors <- sub("FF$","",grDevices::rainbow(nrow(annotations)))
  } else if(colorset == "viridis") {
    colors <- sub("FF$","",viridisLite::viridis(nrow(annotations)))
  } else if(colorset == "magma") {
    colors <- sub("FF$","",viridisLite::magma(nrow(annotations)))
  } else if(colorset == "inferno") {
    colors <- sub("FF$","",viridisLite::inferno(nrow(annotations)))
  } else if(colorset == "plasma") {
    colors <- sub("FF$","",viridisLite::plasma(nrow(annotations)))
  } else if(colorset == "terrain") {
    colors <- sub("FF$","",grDevices::terrain.colors(nrow(annotations)))
  } else if(is.character(colorset)) {
    colors <- grDevices::colorRampPalette(colorset)(nrow(annotations))
  }

  if(color_order == "random") {

    colors <- sample(colors, length(colors))

  }

  annotations <- dplyr::mutate(annotations, color = colors)

  names(annotations) <- paste0(base, c("_label","_id","_color"))

  names(df)[names(df) == col] <- paste0(base,"_label")

  df <- dplyr::left_join(df, annotations, by = paste0(base, "_label"))

  df
}

#'Group annotation columns
#'
#'@param df the annotation dataframe to arrange
#'@param sample_col the column with unique sample ids. Default is "sample_name".
#'@param keep_order a logical value. If FALSE, will sort the annotations alphanumerically by base.
#'
#'
#'@return an annotation data frame with reordered columns
#'
group_annotations <- function(df, sample_col = "sample_name", keep_order = TRUE) {
  labels <- names(df)[grepl("_label",names(df))]
  if(!keep_order) {
    labels <- labels[order(labels)]
  }
  bases <- sub("_label","",labels)

  anno_cols <- c(paste0(rep(bases,each=3),c("_id","_label","_color")))
  extras <- setdiff(names(df),anno_cols)

  anno <- select(df,one_of(c(sample_col,anno_cols,extras)))

}

#' Convert a cl factor object to an annotation data.frame
#'
#' @param cl The cl factor object generated by clustering with scrattch.hicat
#' @param cl.df Optional: A cluster annotation data.frame.
#' @param cl_col The column to use to match cl.df to cl. Default is "cl", but if not found, will try to match using rownames(cl.df).
#' @param base The base of the cluster annotation. Default is "cluster", which will result in "cluster_id","cluster_label", and "cluster_color" columns.
#'
#' @return a data.frame with columns "sample_name", base"_id", base"_label", and base"_color"
#'
cl_to_anno <- function(cl,
                       cl.df = NULL,
                       cl_col = "cl",
                       base = "cluster") {

  anno_names <- c(paste0(base, "_id"), paste0(base,"_label"), paste0(base, "_color"))

  if(class(cl) != "factor") {
    cl_names <- names(cl)
    cl <- as.factor(cl)
    names(cl) <- cl_names
  }

  if(is.null(cl.df)) {
    color_cl <- cl
    levels(color_cl) <- varibow(length(levels(cl)))

    anno <- data.frame(sample_name = names(cl),
                       id = match(cl, levels(cl)),
                       label = as.character(cl),
                       color = as.character(color_cl),
                       stringsAsFactors = FALSE)

    names(anno) <- c("sample_name", anno_names)

  } else {

    # Check to see if there's a cl column in cl.df
    if(!cl_col %in% names(cl.df)) {

      # If not, check to see if levels are stored in the rownames of the df
      unique_rownames <- unique(rownames(cl.df))
      if(sum(unique_rownames %in% levels(cl)) == length(levels(cl))) {
        cl.df[[cl_col]] <- rownames(cl.df)
      } else {
        stop(paste0("No '",cl_col,"' column and rownames don't match levels(cl)."))
      }

    }

    anno_df <- cl.df[, names(cl.df) %in% c(cl_col, anno_names), drop = FALSE]

    anno <- data.frame(sample_name = names(cl),
                       cl = as.character(cl),
                       stringsAsFactors = FALSE)
    names(anno)[2] <- cl_col

    anno <- dplyr::left_join(anno, anno_df, by = cl_col)
    anno <- anno[,!names(anno) == cl_col]

    # Check and fix missing columns
    if(!paste0(base,"_id") %in% names(anno)) {
      anno[[paste0(base,"_id")]] <- match(cl, levels(cl))
    }

    if(!paste0(base,"_label") %in% names(anno)) {
      anno[[paste0(base,"_label")]] <- as.character(cl)
    }

    if(!paste0(base,"_color") %in% names(anno)) {
      color_cl <- cl
      levels(color_cl) <- varibow(length(levels(cl)))
      anno[[paste0(base,"_color")]] <- as.character(color_cl)
    }

  }

  group_annotations(anno)
}


#' Create a generic description file
#'
#' @param dat any data frame that you would like to create a description file for
#' @param names desired names of each element in the description file (default is the column names)
#' @param use_label_columns should only columns containing "_label" be included (default = FALSE)
#'
#' @return a data.frame with columns "base", "name", and "type" for writing to tome
#'
create_desc <- function(dat, name = colnames(dat), use_label_columns = FALSE) {

  if (use_label_columns) {
    dat <- dat[, grepl("_label", colnames(dat))]
    colnames(dat) <- gsub("_label", "", colnames(dat))
  }

  desc <- data.frame(base = colnames(dat), name = name, type = "cat")
  for (i in 1:dim(dat)[2]) if (is.element(class(dat[, i]), c("numeric", "integer"))) {
      desc[i, 3] <- "num"
    }
  desc
}


#' Automatically format an annotation file
#'
#' This takes an anno file as input at any stage and properly annotates it for compatability with
#'   shiny and other scrattch functions.  In particular, it ensures that columns have a label,
#'   an id, and a color, and that there are no factors.  It won't overwrite columns that have
#'   already been properly process.
#'
#' @param anno an existing annotation data frame
#' @param remove_factors should factors be converted to character (currently, should be kept as TRUE)
#'
#' @return an updated data frame that has been automatically annotated properly
#'
auto_annotate <- function(anno, remove_factors = TRUE) {

  ## Convert sample_id to sample_name, if needed
  anno_out <- anno
  colnames(anno_out) <- gsub("sample_id", "sample_name", colnames(anno_out))

  ## Determine which columns are not already formatted?
  cn <- colnames(anno_out)
  convertColumns <- cn[(!grepl("_label", cn)) & (!grepl("_id", cn)) & (!grepl("_color", cn))]
  convertColumns <- setdiff(convertColumns, "sample_name")

  ## Remove all factors
  if (remove_factors) {
    for (cc in convertColumns) if (is.factor(anno_out[, cc])) {
        anno_out[, cc] <- as.character(anno_out[, cc])
      }
  }

  ## Automatically annotate the columns
  for (cc in convertColumns) {
    if (is.numeric(anno_out[, cc])) {
      anno_out <- annotate_num(anno_out, cc)
    } else {
      anno_out <- annotate_cat(anno_out, cc)
    }
  }

  ## Reorganize the anno file (MIGHT NEED TO BE UPDATED)
  anno_out <- group_annotations(anno_out)
}
