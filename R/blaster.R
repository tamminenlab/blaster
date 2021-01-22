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


#' Creates random tmp filename
#'
#' @param length A number
#' @param suffix A string
#' @return A string; random filename
#' @examples
#' create_random_name()
#' create_random_name(length=30, suffix = ".txt")
create_random_name <- function(length = 20, suffix = ".csv")
{
    paste0(
        paste0(sample(c(letters, 1:9),
                      length),
               collapse = ""),
        suffix)
}


#' Runs BLAST algorithm
#'
#' @param query_table A dataframe
#' @param db_table A dataframe
#' @param maxAccepts A number
#' @param maxRejects A number
#' @param minIdentity A number
#' @param strand A string
#' @return A dataframe
#' @export
blast <- function(query_table,
           db_table,
           maxAccepts = 1,
           maxRejects = 16,
           minIdentity = 0.75,
           strand = "both",
           output_to_tmp_file = FALSE)
{
    tmp_file <- create_random_name()
    if (!output_to_tmp_file)
        on.exit(file.remove(tmp_file), add = TRUE)

    pre_blast(query_table,
              db_table,
          tmp_file,
          maxAccepts,
          maxRejects,
          minIdentity,
          strand)

    if (output_to_tmp_file)
        tmp_file
    else
        read.csv(tmp_file,
                 col.names = c("QueryId", "TargetId", "QueryMatchStart",
                               "QueryMatchEnd", "TargetMatchStart",
                               "TargetMatchEnd", "QueryMatchSeq",
                               "TargetMatchSeq", "NumColumns", "NumMatches",
                               "NumMismatches", "NumGaps", "Identity", "Alignment"))
}
