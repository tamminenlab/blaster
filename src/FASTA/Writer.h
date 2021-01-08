#pragma once

#include "../SequenceWriter.h"

namespace FASTA {

template < typename Alphabet >
class Writer : public SequenceWriter< Alphabet > {
public:
  using SequenceWriter< Alphabet >::SequenceWriter;

  static const size_t MaxLineLength = 60;

  Writer< Alphabet >& operator<<( const Sequence< Alphabet >& seq ) {
    auto& out = this->mOutput;
    out << '>' << seq.identifier << std::endl;
    for( size_t i = 0; i < seq.Length(); i += MaxLineLength ) {
      auto sub = seq.Subsequence( i, MaxLineLength );
      out << sub.sequence << std::endl;
    }
    return *this;
  }
};

} // namespace FASTA
