#pragma once

#include "../SequenceReader.h"

namespace FASTQ {

template < typename Alphabet >
class Reader : public SequenceReader< Alphabet > {
public:
  using SequenceReader< Alphabet >::SequenceReader;

  Reader< Alphabet >& operator>>( Sequence< Alphabet >& seq ) {
    ( *mTextReader ) >> seq.identifier;
    ( *mTextReader ) >> seq.sequence;
    ( *mTextReader ) >> seq.quality; // skip plusline
    ( *mTextReader ) >> seq.quality;

    // delete '>'
    seq.identifier.erase( seq.identifier.begin(), seq.identifier.begin() + 1 );

    UpcaseString( seq.sequence ); // atc -> ATC
    UpcaseString( seq.quality );

    return *this;
  }

private:
  using SequenceReader< Alphabet >::mTextReader;
};

} // namespace FASTQ
