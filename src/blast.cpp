#include <Rcpp.h>
#include <string>

#include "FASTA/Reader.h"
#include "Alphabet/DNA.h"
#include "FileFormat.h"

using namespace Rcpp;

// [[Rcpp::export]]
DataFrame blast(DataFrame seq_table) 
{
  std::vector< std::string > ids = seq_table["Id"];
  std::vector< std::string > seqs = seq_table["Seq"];

  
  std::stringstream content;

  for (int i{0}; i < ids.size(); ++i) {
    std::string id{ids[i]};
    std::string seq{seqs[i]};
    content << ">" << id << "\n" << seq << "\n";
  }

  std::istringstream iss( content.str() );

  std::unique_ptr< SequenceReader< DNA > > dbReader( new FASTA::Reader< DNA >( iss ) );
  
  Sequence< DNA > seq;
  SequenceList< DNA > sequences;

  while( !(dbReader->EndOfFile()) ) {
    (*dbReader) >> seq;
    sequences.push_back( std::move( seq ) );
  }

  std::vector< int > new_ids;
  std::vector< std::string > new_seqs;
  int counter{1};
  for (Sequence< DNA > i : sequences) {
    new_ids.push_back( counter );
    new_seqs.push_back( i.sequence );
    ++counter;
      }

  return DataFrame::create(Named("Id") = new_ids,
                           Named("Seq") = new_seqs);
}
