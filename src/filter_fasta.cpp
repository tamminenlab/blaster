
#include <Rcpp.h>
#include <string>

#include "FASTA/Reader.h"
#include "Alphabet/DNA.h"
#include "FileFormat.h"

using namespace Rcpp;

// [[Rcpp::export]]
DataFrame filter_fasta(std::string filename,
                           std::string split_sequence,
                           bool keep_ids = true) 
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
  std::vector< std::string > bcs;
  std::vector< std::string > ribos;
  int counter{1};

  for (Sequence< DNA > i : sequences) {

    std::string seq{i.sequence};
    
    if (seq.find(split_sequence) != std::string::npos) {

      std::string bc{seq.substr(0, seq.find(split_sequence))};
      std::string ribo{seq.substr(seq.find(split_sequence) + split_sequence.size(), seq.size())};

      ids.push_back(id);
      seqs.push_back(seq);
      bcs.push_back(bc);
      ribos.push_back(ribo);
    }
  }

  // CharacterVector z;
  // z = seqs;

  z = DataFrame.create(Named("Id") = ids,
                       Named("BC") = bcs,
                       Named("Ribo") = ribos,
                       Named("Sequence") = seqs);

  return z;
}
