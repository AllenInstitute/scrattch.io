% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read_tome.R
\name{read_tome_anno}
\alias{read_tome_anno}
\title{Read annotations table from a tome file}
\usage{
read_tome_anno(tome, groups = NULL)
}
\arguments{
\item{tome}{The location of the tome file to read.}

\item{groups}{The groups to read - matches column names using grep. Can provide multiple with c(). If NULL, will get all columns. Default is NULL.}
}
\description{
Read annotations table from a tome file
}
