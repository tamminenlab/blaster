% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/blaster.R
\name{blast}
\alias{blast}
\title{Runs BLAST algorithm}
\usage{
blast(
  query_table,
  db_table,
  maxAccepts = 1,
  maxRejects = 16,
  minIdentity = 0.75,
  strand = "both",
  output_to_tmp_file = FALSE
)
}
\arguments{
\item{query_table}{A dataframe}

\item{db_table}{A dataframe}

\item{maxAccepts}{A number}

\item{maxRejects}{A number}

\item{minIdentity}{A number}

\item{strand}{A string}
}
\value{
A dataframe
}
\description{
Runs BLAST algorithm
}