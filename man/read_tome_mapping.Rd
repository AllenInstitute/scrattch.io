% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read_tome.R
\name{read_tome_mapping}
\alias{read_tome_mapping}
\title{Read mapping table from a tome file}
\usage{
read_tome_mapping(tome, mapping_name = NULL, columns = NULL,
  get_all = FALSE)
}
\arguments{
\item{tome}{the location of the tome file to read.}

\item{mapping_name}{the name of the mapping table to read}

\item{columns}{selected columns to read. If NULL, reads all columns. Default = NULL.}

\item{get_all}{logical, whether or not to append all other columns after the specified columns. Default = FALSE.}
}
\description{
Read mapping table from a tome file
}
