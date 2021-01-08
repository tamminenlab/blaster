#pragma once

#include "../SequenceReader.h"

namespace FASTA {

template < typename Alphabet >
class Reader : public SequenceReader< Alphabet > {
public:
  using SequenceReader< Alphabet >::SequenceReader;

  Reader< Alphabet >& operator>>( Sequence< Alphabet >& seq ) {
    std::string identifier, sequence;
    if( mLastLine.empty() ) {
      ( *mTextReader ) >> identifier;
    } else {
      identifier = mLastLine;
    }

    std::string line;
    while( !SequenceReader< Alphabet >::EndOfFile() ) {
      ( *mTextReader ) >> line;

      if( line[ 0 ] == '>' ) {
        mLastLine = line;
        break;
      }

      sequence += line;
    }

    UpcaseString( sequence );
    seq = Sequence< Alphabet >( identifier.substr( 1 ), sequence );

    return *this;
  }

private:
  using SequenceReader< Alphabet >::mTextReader;

  std::string mLastLine;
};

} // namespace FASTA
