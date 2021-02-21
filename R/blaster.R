#' Blaster
#' 
#' Blaster implements an efficient BLAST-like sequence
#' comparison algorithm using native R datatypes.
#' 
#' @docType package
#' @author Manu Tamminen <mavatam.@utu.fi>, Steven Schmid
#' @import Rcpp 
#' @importFrom Rcpp evalCpp
#' @useDynLib blaster
#' @name blaster
NULL


#' Runs BLAST sequence comparison algorithm.
#'
#' @param query A dataframe of the query sequences (containing Id and Seq columns)
#'              or a string specifying the FASTA file of the query sequences. 
#' @param db A dataframe of the database sequences (containing Id and Seq columns)
#'           or a string specifying the FASTA file of the database sequences.
#' @param maxAccepts A number specifying the maximum accepted hits.
#' @param maxRejects A number specifying the maximum rejected hits.
#' @param minIdentity A number specifying the minimal accepted sequence 
#'                    similarity between the query and hit sequences.
#' @param alphabet A string specifying the query and database alphabet: 
#'                 'nucleotide' or 'protein'. Defaults to 'nucleotide'.
#' @param strand A string specifying the strand to search: 'plus', 'minus' or
#'               'both'. Defaults to 'both'. Only affects nucleotide searches.
#' @param output_to_file A boolean specifying the output type. If TRUE, the
#'                       results are written into a temporary file a string
#'                       containing the file name and location is returned.
#'                       Otherwise a dataframe of the results is returned.
#'                       Defaults to FALSE.
#' @return A dataframe or a string. A dataframe is returned by default, containing
#'         the BLAST output in columns QueryId, TargetId, QueryMatchStart, QueryMatchEnd,
#'         TargetMatchStart, TargetMatchEnd, QueryMatchSeq, TargetMatchSeq, NumColumns,
#'         NumMatches, NumMismatches, NumGaps, Identity and Alignment. A string is returned
#'         if 'output_to_file' is set to TRUE. This string points to the temporary file
#'         containing the output table. 
#' @examples
#' 
#' query <- system.file("extdata", "query.fasta", package = "blaster")
#' db <- system.file("extdata", "db.fasta", package = "blaster")
#'
#' blast_table <- blast(query = query, db = db)
#'
#' query <- read_fasta(filename = query)
#' db <- read_fasta(filename = db)
#' blast_table <- blast(query = query, db = db)
#' 
#' prot <- system.file("extdata", "prot.fasta", package = "blaster")
#' prot_blast_table <- blast(query = prot, db = prot, alphabet = "protein")
#' 
#' @export
#' @importFrom utils read.csv
blast <- function(query,
                  db,
                  maxAccepts = 1,
                  maxRejects = 16,
                  minIdentity = 0.75,
                  alphabet = "nucleotide", 
                  strand = "both",
                  output_to_file = FALSE)
{
    tmp_file <- tempfile(fileext = ".csv")
    if (!output_to_file)
        on.exit(if (exists(tmp_file)) file.remove(tmp_file), add = TRUE)

    if (is.data.frame(query)) {
        query_file <- tempfile(fileext = ".fasta")
        write(with(query, paste0(">", Id, "\n", Seq)), query_file)
        query <- query_file
        on.exit(if (exists(query)) file.remove(query), add = TRUE)
    }

    if (is.data.frame(db)) {
        db_file <- tempfile(fileext = ".fasta")
        write(with(db, paste0(">", Id, "\n", Seq)), db_file)
        db <- db_file
        on.exit(if (exists(db)) file.remove(db), add = TRUE)
    }

    if (alphabet == "nucleotide")
        dna_blast(
            query,
            db,
            tmp_file,
            maxAccepts,
            maxRejects,
            minIdentity,
            strand)
    else if (alphabet == "protein")
        protein_blast(
            query,
            db,
            tmp_file,
            maxAccepts,
            maxRejects,
            minIdentity)
    else
        stop("Supported alphabet include 'nucleotide' and 'protein'.")

    if (output_to_file)
        tmp_file
    else
        read.csv(tmp_file,
                 col.names = c(
                     "QueryId", "TargetId", "QueryMatchStart",
                     "QueryMatchEnd", "TargetMatchStart",
                     "TargetMatchEnd", "QueryMatchSeq",
                     "TargetMatchSeq", "NumColumns", "NumMatches",
                     "NumMismatches", "NumGaps", "Identity", "Alignment"))
}


#' Reads the contents of nucleotide or protein FASTA file into a dataframe.
#'
#' @param filename A string specifying the name of the FASTA file to be imported.
#' @param filter An optional string specifying a sequence motif for sequence filtering.
#'               Only keeps those sequences containing this motif. Also splits the
#'               matched sequences and provides the split parts in two additional columns.
#' @param alphabet A string specifying the query and database alphabet: 
#'                 'nucleotide' or 'protein'. Defaults to 'nucleotide'.
#' @param non_standard_chars A string specifying instructions for handling non-standard
#'                           nucleotide or amino acid characters. Options include 'remove',
#'                           'ignore' or throw an 'error'. Defaults to 'error'.
#' @return A dataframe containing FASTA ids (Id column) and sequences (Seq column).
#'         If 'filter' is specified, the split sequences are stored in additional columns
#'         Part1 and Part2.
#' @examples
#' 
#' query <- system.file("extdata", "query.fasta", package = "blaster")
#'
#' query <- read_fasta(filename = query)
#'
#' @export
read_fasta <- function(filename,
                       filter = "",
                       non_standard_chars = "error",
                       alphabet = "nucleotide")
{
    if (alphabet == "nucleotide") {
        read_dna_fasta(
            filename,
            filter,
            non_standard_chars)
    } else if (alphabet == "protein") {
        read_protein_fasta(
            filename,
            filter,
            non_standard_chars)
    } else {
        stop("Supported alphabet include 'nucleotide' and 'protein'.")
    }
}
