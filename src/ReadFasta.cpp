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
    for (int i = 0; i < sequence.size(); ++i) {
      if (alphabet.find(sequence[i]) == alphabet.end())
        continue;
      if (sequence[i] != '\r')
        str_acc = str_acc + sequence[i];
    }
  } else if (non_standard_chars == "ignore") {
    for (int i = 0; i < sequence.size(); ++i) {
      if (sequence[i] != '\r')
        str_acc = str_acc + sequence[i];
    }
  } else if (non_standard_chars == "error") {
    for (int i = 0; i < sequence.size(); ++i) {
      if (alphabet.find(sequence[i]) == alphabet.end())
        stop("Non-standard characters in the file!");
      if (sequence[i] != '\r')
        str_acc = str_acc + sequence[i];
    }
  } else {
    stop("Argument 'non_standard_chars' must be 'remove', 'ignore' or 'error'.");
  }
  return str_acc;
}

std::string process_id(std::string seq_id)
{
  std::string str_acc = "";
  for (int i = 0; i < seq_id.size(); ++i) {
    if (seq_id[i] != '\r')
      str_acc = str_acc + seq_id[i];
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
  std::unordered_set<char> nucleotides {'A', 'T', 'C', 'G', 'W', 'S', 'M', '\r',
                                        'K', 'R', 'Y', 'B', 'D', 'H', 'V', 'N'};
  
  if (filter == "") {
    for (Sequence< DNA > i : sequences) {
      sequence = process_sequence(i.sequence, non_standard_chars, nucleotides);
      id = process_id(i.identifier);
      ids.push_back( id );
      seqs.push_back( sequence );
    }
    return DataFrame::create(Named("Id") = ids,
                             Named("Seq") = seqs);
  }
  else {
    for (Sequence< DNA > i : sequences) {
      sequence = process_sequence(i.sequence, non_standard_chars, nucleotides);
      id = process_id(i.identifier);
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
  std::unordered_set<char> aminoacids {'A', 'R', 'N', 'D', 'C', 'E', 'Q', 'G', 'H',
                                       'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W',
                                       'Y', 'V', '\r'};

  
  if (filter == "") {
    for (Sequence< Protein > i : sequences) {
      sequence = process_sequence(i.sequence, non_standard_chars, aminoacids);
      id = process_id(i.identifier);
      ids.push_back( id );
      seqs.push_back( sequence );
    }
    return DataFrame::create(Named("Id") = ids,
                             Named("Seq") = seqs);
  }
  else {
    for (Sequence< Protein > i : sequences) {
      sequence = process_sequence(i.sequence, non_standard_chars, aminoacids);
      id = process_id(i.identifier);
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
