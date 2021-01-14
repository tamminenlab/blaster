#include <Rcpp.h>
#include <string>
#include <regex>

#include "FASTA/Reader.h"
#include "Alphabet/DNA.h"
#include "FileFormat.h"

using namespace Rcpp;

std::vector<std::string> split(const std::string& input,
                               const std::string& regex) {
  std::regex re(regex);
  std::sregex_token_iterator first{input.begin(), input.end(), re, -1}, last;
  return {first, last};
}

std::string get_first_element(std::string id,
                              std::string id_split_string) {
  std::string new_id;
  std::vector< std::string > id_split;
  if (id_split_string != "") {
    id_split = split(id, id_split_string);
    new_id = id_split[0];
  }
  else new_id = id;
  return new_id;
}


// [[Rcpp::export]]
DataFrame read_fasta(std::string filename,
                     std::string id_split_string = "",
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
      id = get_first_element(i.identifier,
                             id_split_string);
      ids.push_back( id );
      seqs.push_back( sequence );
    }
    return DataFrame::create(Named("Id") = ids,
                             Named("Seq") = seqs);
  }
  else {
    for (Sequence< DNA > i : sequences) {
      sequence = i.sequence;
      id = get_first_element(i.identifier,
                             id_split_string);
      if (sequence.find(filter) != std::string::npos && split) {
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
