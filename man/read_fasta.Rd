#' Read FASTA File
#'
#' Reads the contents of a nucleotide or protein FASTA file into a dataframe.
#'
#' @param filename A string specifying the name of the FASTA file to be imported.
#' @param filter An optional string specifying a sequence motif for sequence filtering.
#'   Only keeps those sequences containing this motif. Also splits the matched sequences
#'   and provides the split parts in two additional columns.
#' @param non_standard_chars A string specifying instructions for handling non-standard
#'   nucleotide or amino acid characters. Options include 'remove', 'ignore' or throw an 'error'.
#'   Defaults to 'error'.
#' @param alphabet A string specifying the query and database alphabet: 'nucleotide' or 'protein'.
#'   Defaults to 'nucleotide'.
#'
#' @return A dataframe containing FASTA ids (Id column) and sequences (Seq column).
#'   If 'filter' is specified, the split sequences are stored in additional columns Part1 and Part2.
#'
#' @examples
#' query <- system.file("extdata", "query.fasta", package = "blaster")
#' query <- read_fasta(filename = query)
#'
#' @export
read_fasta <- function(filename, filter = "", non_standard_chars = "error", alphabet = "nucleotide") {
  # function implementation
}

