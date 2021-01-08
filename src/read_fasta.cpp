
#include <Rcpp.h>
#include <string>

#include "FASTA/Reader.h"
#include "Alphabet/DNA.h"
#include "FileFormat.h"

using namespace Rcpp;

// [[Rcpp::export]]
CharacterVector read_fasta(std::string filename, bool keep_ids = true) 
{
  auto dbReader = DetectFileFormatAndOpenReader< DNA >( filename, FileFormat::FASTA );

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

  CharacterVector z;
  z = seqs;

  if (keep_ids) z.names() = ids;
  
  return z;
}
