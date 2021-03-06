% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read_shelf.R
\name{read_shelf_sample_data}
\alias{read_shelf_sample_data}
\title{Read sample data across a set of tomes that share the same genes}
\usage{
read_shelf_sample_data(tomes, samples, regions = "exon",
  units = "counts", transform = "none", format = "data.frame")
}
\arguments{
\item{tomes}{a character vector of tome file locations}

\item{samples}{A vector of sample names to read. Required.}

\item{regions}{The gene regions to use. Can be "exon", "intron", or "both". Default = "exon".}

\item{units}{The type of values to return. Can be "counts" or "cpm". Default = "counts".}

\item{transform}{Transformation to apply to values. Can be "none", "log", "log2", "log10". Log transforms will add 1 to values before transformation. Default = "none".}

\item{format}{The format of the output. Can be "data.frame", "matrix", or "dgcMatrix" (sparse matrix). Default is "data.frame".}
}
\value{
A data.frame with sample_name as the first column and each subsequent column
containing gene expression values and named for the genes; Or a matrix with columns as genes and rows as samples.
}
\description{
Read sample data across a set of tomes that share the same genes
}
