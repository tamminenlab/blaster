#include <Rcpp.h>
#include <string>
#include <fstream>
#include <set>

#include "FASTA/Reader.h"
#include "Alphabet/DNA.h"
#include "FileFormat.h"

using namespace Rcpp;

std::string process_sequence(std::string sequence,
                             std::string non_standard_chars,
                             std::unordered_set<char> alphabet)
{
  std::string str_acc = "";

  if (non_standard_chars == "remove") {
    for (int i = 0; i < sequence.size(); ++i) 
      if (!(alphabet.find(sequence[i]) == alphabet.end()))
        str_acc = str_acc + sequence[i];
  } else if (non_standard_chars == "ignore") {
    str_acc = sequence;
  } else if (non_standard_chars == "error")  {
    for (int i = 0; i < sequence.size(); ++i) 
      if (alphabet.find(sequence[i]) == alphabet.end())
        stop("Non-standard characters in the file!");
      str_acc = sequence;
  } else {
    stop("Argument 'non_standard_chars' must be 'remove', 'ignore' or 'error'.");
  }
  return str_acc;
}

// Read the contents of a nucleotide Fasta file into a DataFrame
// [[Rcpp::export]]
DataFrame read_dna_fasta(std::string filename,
                         std::string filter,
                         std::string non_standard_chars) 
{
  std::ifstream f(filename);
  if (!f.good())
    stop("File does not exist.");
  f.close(); 
  
  std::unique_ptr< SequenceReader< DNA > > dbReader( new FASTA::Reader< DNA >( filename ) );

  Sequence< DNA > seq;
  SequenceList< DNA > sequences;

  while( !(dbReader->EndOfFile()) ) {
    (*dbReader) >> seq;
    sequences.push_back( std::move( seq ) );
  }

  std::string part1;
  std::string part2;
  std::string sequence;
  std::vector< std::string > ids;
  std::vector< std::string > seqs;
  std::vector< std::string > parts1;
  std::vector< std::string > parts2;
  std::unordered_set<char> nucleotides {'A', 'T', 'C', 'G', 'U', 'W', 'S', 'M',
                                        'K', 'R', 'Y', 'B', 'D', 'H', 'V', 'N'};
  
  if (filter == "") {
    for (Sequence< DNA > i : sequences) {
      ids.push_back( i.identifier );
      seqs.push_back( process_sequence(i.sequence, non_standard_chars, nucleotides) );
    }
    return DataFrame::create(Named("Id") = ids,
                             Named("Seq") = seqs);
  }
  else {
    for (Sequence< DNA > i : sequences) {
      sequence = process_sequence(i.sequence, non_standard_chars, nucleotides);
      if (sequence.find(filter) != std::string::npos) {
        part1 = sequence.substr(0, sequence.find(filter));
        part2 = sequence.substr(sequence.find(filter) + filter.size(),
                                sequence.size());
        ids.push_back( i.identifier );
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

// Read the contents of a protein Fasta file into a DataFrame
// [[Rcpp::export]]
DataFrame read_protein_fasta(std::string filename,
                         std::string filter,
                         std::string non_standard_chars) 
{
  std::ifstream f(filename);
  if (!f.good())
    stop("File does not exist.");
  f.close(); 
  
  std::unique_ptr< SequenceReader< Protein > > dbReader( new FASTA::Reader< Protein >( filename ) );

  Sequence< Protein > seq;
  SequenceList< Protein > sequences;

  while( !(dbReader->EndOfFile()) ) {
    (*dbReader) >> seq;
    sequences.push_back( std::move( seq ) );
  }

  std::string part1;
  std::string part2;
  std::string sequence;
  std::vector< std::string > ids;
  std::vector< std::string > seqs;
  std::vector< std::string > parts1;
  std::vector< std::string > parts2;
  std::unordered_set<char> aminoacids {'A', 'R', 'N', 'D', 'C', 'E', 'Q', 'G', 'H',
                                       'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W',
                                       'Y', 'V'};
  
  if (filter == "") {
    for (Sequence< Protein > i : sequences) {
      ids.push_back( i.identifier );
      seqs.push_back( process_sequence(i.sequence, non_standard_chars, aminoacids) );
    }
    return DataFrame::create(Named("Id") = ids,
                             Named("Seq") = seqs);
  }
  else {
    for (Sequence< Protein > i : sequences) {
      sequence = process_sequence(i.sequence, non_standard_chars, aminoacids);
      if (sequence.find(filter) != std::string::npos) {
        part1 = sequence.substr(0, sequence.find(filter));
        part2 = sequence.substr(sequence.find(filter) + filter.size(),
                                sequence.size());
        ids.push_back( i.identifier );
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
