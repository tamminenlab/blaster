#include <Rcpp.h>
#include <string>

#include "FASTA/Reader.h"
#include "Alphabet/DNA.h"
#include "FileFormat.h"

using namespace Rcpp;

// [[Rcpp::export]]
DataFrame read_fasta(std::string filename) 
{
  std::unique_ptr< SequenceReader< DNA > > dbReader( new FASTA::Reader< DNA >( filename ) );

  Sequence< DNA > seq;
  SequenceList< DNA > sequences;

  while( !(dbReader->EndOfFile()) ) {
    (*dbReader) >> seq;
    sequences.push_back( std::move( seq ) );
  }

  std::vector< std::string > ids;
  std::vector< std::string > seqs;
  for (Sequence< DNA > i : sequences) {
    ids.push_back( i.identifier );
    seqs.push_back( i.sequence );
      }

  return DataFrame::create(Named("Id") = ids,
                           Named("Seq") = seqs);
}
