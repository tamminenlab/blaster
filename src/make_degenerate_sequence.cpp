#include <string>
#include <vector>
#include <map>
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
StringVector make_degenerate_sequence(StringVector sequences,
                                      int sequence_length,
                                      double cutoff = 0.1)
{
  std::vector< std::vector< char > > aln_matrix;

  for (int i = 0 ; i < sequences.size() ; ++i) {
    if (sequences[i].size() == sequence_length ) {
      std::vector<char> chars(sequences[i].begin(), sequences[i].end());
      aln_matrix.push_back(chars);
    }
  }

  if (aln_matrix.size() == 0)
    stop("Zero sequences of the specified length provided.");

  std::vector< std::vector< char> > trans_vec(aln_matrix[0].size(),
                                              std::vector< char >());

  for (int i = 0; i < aln_matrix.size(); i++)
    for (int j = 0; j < aln_matrix[i].size(); j++)
      trans_vec[j].push_back(aln_matrix[i][j]);

  std::map< std::string, char > deg_nucleotides {
    {"A", 'A'},
    {"T", 'T'},
    {"C", 'C'},
    {"G", 'G'},
    {"AT", 'W'},
    {"CG", 'S'},
    {"CG", 'S'},
    {"AC", 'M'},
    {"GT", 'K'},
    {"AG", 'R'},
    {"CT", 'Y'},
    {"CGT", 'B'},
    {"AGT", 'D'},
    {"ACG", 'H'},
    {"ACT", 'V'},
    {"ACGT", 'N'}
  };

  double aln_len(trans_vec[0].size());
  std::vector< char > deg_seq;
  for (std::vector< char > x : trans_vec) {

    std::map< char, int > char_hash;
    for (char y : x)
      char_hash[y]++;

    std::vector< char > aln_acc;
    for (auto i : char_hash)
      if (i.second / aln_len > cutoff)
        aln_acc.push_back(i.first);

    std::sort(aln_acc.begin(), aln_acc.end());
    std::string aln_sorted(aln_acc.begin(), aln_acc.end());
    deg_seq.push_back(deg_nucleotides[aln_sorted]);
  }

  std::string deg_string(deg_seq.begin(), deg_seq.end());

  return CharacterVector::create(deg_string);

}
