#include <Rcpp.h>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <regex>
#include <set>

#include "Utils.h"

using namespace Rcpp;

// [[Rcpp::export]]
List process_blast_table(std::string filename) 
{
  std::fstream testFile(filename);    
  std::string line;
  std::vector< std::string > split_line;
  std::vector< std::string > split_query;
  std::set< std::string > primer_set;
  std::set< std::string > seq_set;

  while(getline(testFile, line)){
    split_line = split(line, ",");
    std::string target{split_line[1] + "_" + split_line[4]};
    std::string new_primer_item{split_line[0] + "," + target};
    if (split_line[0] != target)
      primer_set.insert(new_primer_item);

    split_query = split(split_line[0], "_");
    std::string new_seq_item{split_query[0] + "," + split_line[1]};
    if (split_query[0] != split_line[1])
      seq_set.insert(new_seq_item);
  }
  
  std::set< std::string >::iterator primer_it = --primer_set.end();
  primer_set.erase(primer_it); 

  std::vector < std::string > primer_left;
  std::vector < std::string > primer_right;

  for (std::set< std::string >::iterator primer_it = primer_set.begin(); primer_it != primer_set.end(); ++primer_it) {
    line = *primer_it;
    split_line = split(line, ",");
    primer_left.push_back(split_line[0]);
    primer_right.push_back(split_line[1]);
  }

  std::set< std::string >::iterator seq_it = --seq_set.end();
  seq_set.erase(seq_it); 

  std::vector < std::string > seq_left;
  std::vector < std::string > seq_right;

  for (std::set< std::string >::iterator seq_it = seq_set.begin(); seq_it != seq_set.end(); ++seq_it) {
    line = *seq_it;
    split_line = split(line, ",");
    seq_left.push_back(split_line[0]);
    seq_right.push_back(split_line[1]);
  }
  
  DataFrame primer_table = DataFrame::create(Named("Left") = primer_left,
                                             Named("Right") = primer_right);

  DataFrame seq_table = DataFrame::create(Named("Left") = seq_left,
                                          Named("Right") = seq_right);

  return List::create(primer_table, seq_table);

}
