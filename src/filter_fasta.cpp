#include <Rcpp.h>
#include <string>

#include "FASTA/Reader.h"
#include "Alphabet/DNA.h"
#include "FileFormat.h"

using namespace Rcpp;

// [[Rcpp::export]]
DataFrame filter_fasta(std::string filename,
                       std::string by,
                       bool split = true) 
{
  auto dbReader = DetectFileFormatAndOpenReader< DNA >( filename, FileFormat::FASTA );

  Sequence< DNA > seq;
  SequenceList< DNA > sequences;

  while( !(dbReader->EndOfFile()) ) {
    (*dbReader) >> seq;
    sequences.push_back( std::move( seq ) );
  }

  if (split) {
    std::vector< std::string > ids;
    std::vector< std::string > seqs;
    std::vector< std::string > bcs;
    std::vector< std::string > ribos;

    for (Sequence< DNA > i : sequences) {

      std::string seq{i.sequence};
    
      if (seq.find(by) != std::string::npos) {

        std::string bc{seq.substr(0, seq.find(by))};
        std::string ribo{seq.substr(seq.find(by) + by.size(),
                                    seq.size())};

        ids.push_back(i.identifier);
        bcs.push_back(bc);
        ribos.push_back(ribo);
      }
    }

    return DataFrame::create(Named("Id") = ids,
                             Named("Part1") = bcs,
                             Named("Part2") = ribos);
  }
  else
    {
    
      std::vector< std::string > ids;
      std::vector< std::string > seqs;

      for (Sequence< DNA > i : sequences) {

        std::string seq{i.sequence};
    
        if (seq.find(by) != std::string::npos) {
          ids.push_back(i.identifier);
          seqs.push_back(seq);
        }
      }

      return DataFrame::create(Named("Id") = ids,
                               Named("Seq") = seqs);
    }
}
