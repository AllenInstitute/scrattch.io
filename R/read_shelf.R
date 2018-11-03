#' Read gene data across a set of tomes that share the same genes
#'
#' @param tomes a character vector of tome file locations
#' @param genes A vector of gene names to read. Required.
#' @param regions The gene regions to use. Can be "exon", "intron", or "both". Default = "exon".
#' @param units The type of values to return. Can be "counts" or "cpm". Default = "counts".
#' @param transform Transformation to apply to values. Can be "none", "log", "log2", "log10". Log transforms will add 1 to values before transformation. Default = "none".
#' @param format The format of the output. Can be "data.frame", "matrix", or "dgcMatrix" (sparse matrix). Default is "data.frame".
#'
#' @return A data.frame with sample_name as the first column and each subsequent column
#' containing gene expression values and named for the genes; Or a matrix with columns as genes and rows as samples.
#'
read_shelf_gene_data <- function(tomes,
                                 genes,
                                 regions = "exon",
                                 units = "counts",
                                 transform = "none",
                                 format = "data.frame") {

  available_genes_list <- purrr::map(tomes,
                                function(tome) {
                                  read_tome_gene_names(tome)
                                })

  common_genes <- Reduce(intersect, c(list(genes), available_genes_list))

  missing_genes <- setdiff(genes, common_genes)

  if(length(missing_genes) > 0) {
    cat("Some genes not found in all tomes: ",paste(missing_genes, collapse = ", "), "\n")
  }

  purrr::map_dfr(tomes,
                 function(tome) {
                   read_tome_gene_data(tome,
                                       common_genes,
                                       regions,
                                       units,
                                       transform,
                                       format)
                 })

}


#' Read sample data across a set of tomes that share the same genes
#'
#' @param tomes a character vector of tome file locations
#' @param samples A vector of sample names to read. Required.
#' @param regions The gene regions to use. Can be "exon", "intron", or "both". Default = "exon".
#' @param units The type of values to return. Can be "counts" or "cpm". Default = "counts".
#' @param transform Transformation to apply to values. Can be "none", "log", "log2", "log10". Log transforms will add 1 to values before transformation. Default = "none".
#' @param format The format of the output. Can be "data.frame", "matrix", or "dgcMatrix" (sparse matrix). Default is "data.frame".
#'
#' @return A data.frame with sample_name as the first column and each subsequent column
#' containing gene expression values and named for the genes; Or a matrix with columns as genes and rows as samples.
#'
read_shelf_sample_data <- function(tomes,
                                 samples,
                                 regions = "exon",
                                 units = "counts",
                                 transform = "none",
                                 format = "data.frame") {

  available_genes_list <- purrr::map(tomes,
                                     function(tome) {
                                       read_tome_gene_names(tome)
                                     })

  common_genes <- Reduce(intersect, available_genes_list)

  missing_genes <- setdiff(genes, common_genes)

  if(length(missing_genes) > 0) {
    cat(length(missing_genes), " genes not found in all tomes.\n")
  }

  purrr::map_dfc(tomes,
                 function(tome) {
                   tome_available_samples <- read_tome_sample_names(tome)
                   tome_samples <- intersect(samples, tome_available_samples)
                   data <- read_tome_sample_data(tome,
                                                 tome_samples,
                                                 regions,
                                                 units,
                                                 transform,
                                                 format)
                   if(format == "data.frame") {
                     data[data$gene_name %in% common_genes,]
                   } else {
                     data[common_genes, ]
                   }
                 })

}
