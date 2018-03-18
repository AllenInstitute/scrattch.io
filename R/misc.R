collapse_along <- function(x,
                           collapse = "/") {

  out <- vector(length = length(x))

   for(i in 1:length(x)) {
    out[i] <- paste(x[1:i], collapse = collapse)
  }

  out
}

flip_table <- function(df, gene_col = "gene", id_col = "sample_id") {

  if(gene_col %in% names(df)) {
    genes <- unlist(df[,gene_col])
    df_t <- t(df[,names(df) != gene_col])
    samples <- rownames(df_t)
    df_out <- cbind(samples, as.data.frame(df_t))
    names(df_out) <- c(id_col,genes)
    rownames(df_out) <- NULL

    print("(╯°□°）╯︵ ┻━┻")

    df_out

  } else if(id_col %in% names(df)) {

    samples <- unlist(df[,id_col])
    df_t <- t(df[,names(df) != id_col])
    genes <- rownames(df_t)
    df_out <- cbind(genes, as.data.frame(df_t))
    names(df_out) <- c(gene_col, samples)
    rownames(df_out) <- NULL

    print("(╯°□°）╯︵ ┻━┻")

    df_out

  } else {
    print(paste("No column named",gene_col,"or",id_col,"found."))
  }

}
