#include <Rcpp.h>
#include <string>
#include <fstream>

#include "FASTA/Reader.h"
#include "Alphabet/DNA.h"
#include "FileFormat.h"

using namespace Rcpp;

//' Read the contents of a Fasta file into a DataFrame
//' 
//' @param filename A string; name of the imported Fasta file.
//' @param filter A string; only include those entries which contain a sequence motif specified in this argument.
//' @return A dataframe with Id and Seq columns
//' @examples
//' 
//' query <- read_fasta(filename = "extdata/query.fasta")
//' 
//' query <- read_fasta(filename = "extdata/db.fasta", filter = "TGGTTGAGG")
//' 
//' @export
// [[Rcpp::export]]
DataFrame read_fasta(std::string filename,
                     std::string filter = "") 
{
  std::unique_ptr< SequenceReader< DNA > > dbReader( new FASTA::Reader< DNA >( filename ) );

  Sequence< DNA > seq;
  SequenceList< DNA > sequences;
  std::vector< std::string > id_split;

  while( !(dbReader->EndOfFile()) ) {
    (*dbReader) >> seq;
    sequences.push_back( std::move( seq ) );
  }

  std::string id;
  std::string sequence;
  std::string part1;
  std::string part2;
  std::vector< std::string > ids;
  std::vector< std::string > seqs;
  std::vector< std::string > parts1;
  std::vector< std::string > parts2;
  
  if (filter == "") {
    for (Sequence< DNA > i : sequences) {
      sequence = i.sequence;
      id = i.identifier;
      ids.push_back( id );
      seqs.push_back( sequence );
    }
    return DataFrame::create(Named("Id") = ids,
                             Named("Seq") = seqs);
  }
  else {
    for (Sequence< DNA > i : sequences) {
      sequence = i.sequence;
      id = i.identifier;
      if (sequence.find(filter) != std::string::npos) {
        part1 = sequence.substr(0, sequence.find(filter));
        part2 = sequence.substr(sequence.find(filter) + filter.size(),
                                sequence.size());
        ids.push_back( id );
        seqs.push_back( sequence );
        parts1.push_back( part1 );
        parts2.push_back( part2 );
      }
    }
  }
  return DataFrame::create(Named("Id") = ids,
                           Named("Seq") = seqs,
                           Named("Part1") = parts1,
                           Named("Part2") = parts2);
}
