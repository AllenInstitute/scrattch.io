% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/write_10x.R
\name{write_dgCMatrix_h5}
\alias{write_dgCMatrix_h5}
\title{Write a dgCMatrix to an h5 file similar to cellRanger format}
\usage{
write_dgCMatrix_h5(mat, cols_are = "gene_names", h5_target,
  ref_name = "mm10-1.2.0_premrna", gene_ids = NULL)
}
\arguments{
\item{mat}{a dgCMatrix to write.}

\item{cols_are}{Whether columns are "gene_names" or "sample_names". If "gene_names", mat will be transposed to match 10X conventions.}

\item{h5_target}{The target .h5 file for output.}

\item{ref_name}{Reference name for storing the data}

\item{gene_ids}{If available, ENSEMBL IDs for each gene}
}
\value{
an .h5 file with these mappings from dgCMatrix -> .h5:
\itemize{
  \item colnames(mat) -> /ref_name/barcodes (after transposition if cols_are == "gene_names")
  \item rownames(mat) -> /ref_name/gene_names (after transposition if cols_are == "gene_names")
  \item mat@x -> /ref_name/data
  \item mat@i -> /ref_name/indices
  \item mat@p -> /ref_name/indptr
  \item gene_ids -> /ref_name/gene
}
}
\description{
Write a dgCMatrix to an h5 file similar to cellRanger format
}
