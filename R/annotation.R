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
#'   are "varibow" (default), "rainbow","viridis","inferno","magma", and "terrain"
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
#' @param name desired names of each element in the description file (default is the column names)
#' @param use_label_columns should only columns containing "_label" be included (default = FALSE)
#' @param start_columns character vector of variables to include first in the list 
#'	(default = NULL, but "cluster" would be a common choice)
#'
#' @return a data.frame with columns "base", "name", and "type" for writing to tome
#'
create_desc <- function(dat, name = colnames(dat), use_label_columns = FALSE, start_columns = NULL) {
  dat <- as.data.frame(dat)
  if (use_label_columns) {
    dat <- dat[, grepl("_label", colnames(dat))]
    colnames(dat) <- gsub("_label", "", colnames(dat))
  }

  desc <- data.frame(base = colnames(dat), name = name, type = "cat")
  for (i in 1:dim(dat)[2]) if (is.element(class(dat[, i]), c("numeric", "integer"))) {
      desc[i, 3] <- "num"
  }
	
  ## Reorder colums as requested
  cn   <- c(intersect(start_columns,desc$base),setdiff(desc$base,start_columns))
  desc <- desc[match(cn,desc$base),]
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
#' @param scale_num should color scaling of numeric values be "predicted" (default and highly recommended;
#'   will return either "linear" or "log10" depending on scaling), "linear","log10","log2", or "zscore".
#' @param na_val_num The value to use to replace NAs for numeric columns. default = 0.
#' @param colorset_num A vector of colors to use for the color gradient.
#'   default = c("darkblue","white","red")
#' @param sort_label_cat a logical value to determine if the data in category columns
#'   should be arranged alphanumerically before ids are assigned. default = T.
#' @param na_val_cat The value to use to replace NAs in category and factor variables.
#'   default = "ZZ_Missing".
#' @param colorset_cat The colorset to use for assigning category and factor colors.
#'   Options are "varibow" (default), "rainbow","viridis","inferno","magma", and "terrain"
#' @param color_order_cat The order in which colors should be assigned for cat and
#'   factor variables. Options are "sort" and "random". "sort" (default) assigns colors
#'   in order; "random" will randomly assign colors.
#'
#' @return an updated data frame that has been automatically annotated properly
#'
auto_annotate <- function(anno, scale_num = "predicted", na_val_num = 0,
                          colorset_num = c("darkblue","white","red"),
                          sort_label_cat = TRUE, na_val_cat = "ZZ_Missing",
                          colorset_cat = "varibow", color_order_cat = "sort") {

  ## Convert sample_id to sample_name, if needed
  anno_out <- anno
  if(!is.element("sample_name",colnames(anno_out))){
    colnames(anno_out) <- gsub("sample_id", "sample_name", colnames(anno_out))
  }
  
  ## Determine which columns are not already formatted?
  cn <- colnames(anno_out)
  convertColumns <- cn[(!grepl("_label", cn)) & (!grepl("_id", cn)) & (!grepl("_color", cn))]
  convertColumns <- setdiff(convertColumns, "sample_name")
  convertColumns <- setdiff(convertColumns,gsub("_label","",cn[grepl("_label",cn)])) # Avoid column duplication error

  ## Automatically annotate the columns
  for (cc in convertColumns) {
    value <- anno_out[, cc]
		if (sum(!is.na(value))==0)
		    value = rep("N/A",length(value))  # Account for all NA values
    if (is.numeric(value)) {
	    if(length(table(value))==1)  
        value = jitter(value,0.000001)    # Avoid entirely constant values
	    val2 <- value[!is.na(value)]        # To avoid calculuating mean, max, min on NAs
      if (is.element(scale_num, c("linear", "log10", "log2", "zscore"))) {
        # If scale_num is pre-set for all numeric values...
        anno_out <- annotate_num(df = anno_out, col = cc, scale = scale_num,
                                 na_val = na_val_num, colorset = colorset_num)
      } else {
        # If scale_num is predicted...
        scalePred <- ifelse(min(val2) < 0, "linear", "log10") # Avoid NA values for log scale
        if ((max(val2 + 1) / min(val2 + 1)) < 100) {
          scalePred <- "linear"
        }  # Use log scale if large range
        if (mean((val2 - min(val2)) / diff(range(val2))) < 0.01) {
          scalePred <- "log10"
        }  # Use log scale if large skew towards high end
        anno_out <- annotate_num(df = anno_out, col = cc, scale = scalePred,
                                 na_val = na_val_num, colorset = colorset_num)
      }
    } else {
      if (is.factor(value)){
        anno_out <- annotate_factor(df = anno_out, col = cc, base = cc, na_val = na_val_cat,
                                    colorset = colorset_cat, color_order = color_order_cat)
      } else {
        anno_out <- annotate_cat(df = anno_out, col = cc, base = cc, na_val = na_val_cat,
                                 colorset = colorset_cat, color_order = color_order_cat,
                                 sort_label = sort_label_cat)
      }
    }
  }

  ## Reorganize the anno file (MIGHT NEED TO BE UPDATED)
  anno_out <- group_annotations(anno_out)

  anno_out
}


#' Generate colors and ids for categorical annotations that are factors
#'
#' @param df data frame to annotate
#' @param col name of the factor column to annotate
#' @param base base name for the annotation, which wil be used in the desc
#'   table. If not provided, will use col as base.
#' @param na_val The value to use to replace NAs. default = "ZZ_Missing".
#' @param colorset The colorset to use for assigning category colors. Options
#'   are "varibow" (default), "rainbow","viridis","inferno","magma", and "terrain"
#' @param color_order The order in which colors should be assigned. Options are
#'   "sort" and "random". "sort" assigns colors in order; "random" will randomly
#'   assign colors.
#'
#' @return A modified data frame: the annotated column will be renamed
#'   base_label, and base_id and base_color columns will be appended
#'
annotate_factor <- function(df,
                            col = NULL, base = NULL,
                            na_val = "ZZ_Missing",
                            colorset = "varibow", color_order = "sort") {

  # library(dplyr)
  # library(lazyeval)
  # library(viridisLite)

  if (class(try(is.character(col), silent = T)) == "try-error") {
    col <- lazyeval::expr_text(col)
  } else if (class(col) == "NULL") {
    stop("Specify a column (col) to annotate.")
  }

  if (class(try(is.character(base), silent = T)) == "try-error") {
    base <- lazyeval::expr_text(base)
  } else if (class(base) == "NULL") {
    base <- col
  }

  if (!is.factor(df[[col]])) {
    df[[col]] <- as.factor(df[[col]])
  }

  # Convert NA values and add NA level to the end
  if (sum(is.na(df[[col]])) > 0) {
    lev <- c(levels(df[[col]]), na_val)
    levels(df[[col]]) <- lev
    df[[col]][is.na(df[[col]])] <- na_val
  }

  x <- df[[col]]

  annotations <- data.frame(label = as.character(levels(x)), stringsAsFactors = F)

  annotations <- annotations %>%
    dplyr::mutate(id = 1:n())

  if (colorset == "varibow") {
    colors <- varibow(nrow(annotations))
  } else if (colorset == "rainbow") {
    colors <- sub("FF$", "", grDevices::rainbow(nrow(annotations)))
  } else if (colorset == "viridis") {
    colors <- sub("FF$", "", viridisLite::viridis(nrow(annotations)))
  } else if (colorset == "magma") {
    colors <- sub("FF$", "", viridisLite::magma(nrow(annotations)))
  } else if (colorset == "inferno") {
    colors <- sub("FF$", "", viridisLite::inferno(nrow(annotations)))
  } else if (colorset == "plasma") {
    colors <- sub("FF$", "", viridisLite::plasma(nrow(annotations)))
  } else if (colorset == "terrain") {
    colors <- sub("FF$", "", grDevices::terrain.colors(nrow(annotations)))
  } else if (is.character(colorset)) {
    colors <- grDevices::colorRampPalette(colorset)(nrow(annotations))
  }

  if (color_order == "random") {
    colors <- sample(colors, length(colors))
  }

  annotations <- dplyr::mutate(annotations, color = colors)

  names(annotations) <- paste0(base, c("_label", "_id", "_color"))

  names(df)[names(df) == col] <- paste0(base, "_label")

  df[[paste0(col,"_label")]] <- as.character(df[[paste0(col,"_label")]]) # convert the factor to a character in the anno

  df <- dplyr::left_join(df, annotations, by = paste0(base, "_label"))

  df
}
