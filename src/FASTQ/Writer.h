#pragma once

#include "../SequenceWriter.h"

namespace FASTQ {

template < typename Alphabet >
class Writer : public SequenceWriter< Alphabet > {
public:
  using SequenceWriter< Alphabet >::SequenceWriter;

  Writer< Alphabet >& operator<<( const Sequence< Alphabet >& seq ) {
    auto& out = this->mOutput;
    out << '@' << seq.identifier << std::endl;
    out << seq.sequence << std::endl;
    out << '+' << std::endl;
    out << seq.quality << std::endl;
    return *this;
  }
};

} // namespace FASTQ
