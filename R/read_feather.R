#' Read data from a directory of feather files
#'
#' @param feather_dir Directory containing feather files.
#' @param genes Genes to read
#' @param group_by desc_bases to use for grouping samples
#' @param group_ids ID values to use for filtering
#'
get_feather_data <- function(feather_dir, genes, group_by, group_ids) {

  library(dplyr)
  library(feather)

  data_file <- paste0(feather_dir, "/data.feather")
  anno_file <- paste0(feather_dir, "/anno.feather")

  data <- feather::feather(data_file)

  # Read annotations and convert factors
  anno <- feather::read_feather(anno_file) %>%
    dplyr::mutate_if(is.factor, as.character)

  # If an _id column was a factor, it's now a character. Convert to numeric for sorting.
  id_cols <- names(anno)[grepl("_id$", names(anno)) & names(anno) != "sample_id"]
  anno[id_cols] <- lapply(anno[id_cols], as.numeric)

  # Check the provided genes against the column names in data_file
  data_names <- names(data)

  if(sum(genes %in% data_names) != length(genes)) {
    # Report if names don't match after ignorning case
    not_found <- genes[!toupper(genes) %in% toupper(data_names)]

    warning(paste(paste0(not_found, collapse = ", "), "not found in feather data!"))

    # Update genes to use names as formatted in data
    genes <- data_names[toupper(data_names) %in% toupper(genes)]
  }

  # Find column indexes for sample_id and the matched genes
  # This seems to be faster than calling data[,c("sample_id",genes)] directly
  data_cols <- which(data_names %in% c("sample_id", genes))

  # Read the data from the data feather file into memory
  gene_data <- data[,data_cols]

  # Change - to . in column names and genes
  colnames(gene_data) <- gsub("-",".",colnames(gene_data))
  genes <- gsub("-",".",genes)

  # rename the _id, _label, and _color for the group_by values for use in plotting
  all_anno <- anno %>%
    dplyr::rename_("plot_id" = paste0(group_by,"_id"),
            "plot_label" = paste0(group_by,"_label"),
            "plot_color" = paste0(group_by,"_color"))

  # use the group_ids to retain the order provided by the group_ids argument
  cluster_order <- data.frame(group_ids = group_ids) %>%
    dplyr::mutate(cluster_x = 1:n())

  # Filter and order the rows
  data <- dplyr::left_join(all_anno, gene_data, by = "sample_id") %>%
    dplyr::filter(plot_id %in% group_ids) %>%
    dplyr::left_join(cluster_order, by = c("plot_id" = "group_ids")) %>%
    dplyr::arrange(cluster_x) %>%
    dplyr::mutate(xpos = 1:n()) %>%
    dplyr::select(-plot_id) %>%
    dplyr::rename_("plot_id" = "cluster_x")

  return(data)
}

feather_to_list <- function(feather_dir = NULL, oldformat = F) {

  library(feather)

  if(is.null(feather_dir)) {
    stop("feather directory required.")
  }

  if(!dir.exists(feather_dir)) {
    stop("feather directory doesn't exist")
  }

  datafile <- paste0(feather_dir,"/data.feather")
  annofile <- paste0(feather_dir,"/anno.feather")
  descfile <- paste0(feather_dir,"/desc.feather")

  if(!file.exists(datafile)) {
    stop("data file not found.")
  } else if(!file.exists(annofile)) {
    stop("anno file not found")
  } else if(!file.exists(descfile)) {
    stop("desc file not found.")
  }

  cat("Reading ",descfile,"\n")
  desc <- feather::read_feather(descfile)

  cat("Reading ",annofile,"\n")
  anno <- feather::read_feather(annofile)

  cat("Reading ",datafile,"\n")
  main <- feather::read_feather(datafile)

  if(oldformat == T) {
    # Transforming for compatibility with get_list_data
    # This step will slow things down substantially.
    # In future, better to change get_list_data.
    cat("Transforming data table for compatibility\n")
    gn <- colnames(data)
    data <- t(data)
    data <- cbind(gene = gn, data)
    colnames(data) <- c("gene",anno$sample_id)
  }

  out_list <- list(anno = anno,
                   data = data,
                   desc = desc)

  return(out_list)

}
