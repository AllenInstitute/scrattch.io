#' Cumulatively collapse along a vector
#'
#' @param x The character vector to collapse
#' @param collaps The character to use for collapsing
#'
#' @examples
#'
#' x <- c("","data","exon","i")
#' collapse_along(x)
#'
collapse_along <- function(x,
                           collapse = "/") {

  out <- vector(length = length(x))

   for(i in 1:length(x)) {
    out[i] <- paste(x[1:i], collapse = collapse)
  }

  out
}


#' Transpose a gene x sample data.frame without losing a gene_name or sample_name column
#'
#' @param df The data.frame to transpose
#' @param gene_col The column used for gene names. Default = "gene_name".
#' @param sample_col The column used for sample names. Default = "sample_name".
#'
flip_table <- function(df,
                       gene_col = "gene_name",
                       sample_col = "sample_name") {

  if(gene_col %in% names(df)) {
    genes <- unlist(df[,gene_col])
    df_t <- t(df[,names(df) != gene_col])
    samples <- rownames(df_t)
    df_out <- cbind(samples, as.data.frame(df_t))
    names(df_out) <- c(sample_col,genes)
    rownames(df_out) <- NULL

    print("(╯°□°）╯︵ ┻━┻")

    df_out

  } else if(sample_col %in% names(df)) {

    samples <- unlist(df[,sample_col])
    df_t <- t(df[,names(df) != sample_col])
    genes <- rownames(df_t)
    df_out <- cbind(genes, as.data.frame(df_t))
    names(df_out) <- c(gene_col, samples)
    rownames(df_out) <- NULL

    print("(╯°□°）╯︵ ┻━┻")

    df_out

  } else {
    print(paste("No column named",gene_col,"or",sample_col,"found."))
  }

}
