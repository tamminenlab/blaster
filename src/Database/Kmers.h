#pragma once

#include "../Sequence.h"
#include "../Utils.h"

#include <functional>

using Kmer = uint32_t;
const Kmer AmbiguousKmer = ( Kmer )-1;

template< typename Alphabet >
class Kmers {
public:
  using Callback = const std::function< void( const Kmer, const size_t ) >;

  Kmers( const Sequence< Alphabet >& ref, const size_t length ) : mRef( ref ) {
    mLength = std::min( { length, mRef.Length(), sizeof( Kmer ) * 8 / BitMapPolicy< Alphabet >::NumBits } );
  }

  void ForEach( const Callback& block ) const {
    const char* ptr = mRef.sequence.data();

    auto bitIndex = []( const size_t pos ) {
      return ( pos * BitMapPolicy< Alphabet >::NumBits ) % ( sizeof( Kmer ) * 8 );
    };

    auto bitMapNucleotide = []( const char base ) {
      return BitMapPolicy< Alphabet >::BitMap( base );
    };

    // First kmer
    size_t lastAmbigIndex = ( size_t ) -1;
    Kmer   kmer           = 0;
    for( size_t k = 0; k < mLength; k++ ) {
      int8_t val = bitMapNucleotide( *ptr );
      if( val < 0 ) {
        lastAmbigIndex = k;
      } else {
        kmer |= ( val << bitIndex( k ) );
      }
      ptr++;
    }

    if( lastAmbigIndex == ( size_t ) -1 ) {
      block( kmer, 0 );
    } else {
      block( AmbiguousKmer, 0 );
    }

    // For each consecutive kmer, shift window by one
    size_t maxFrame = mRef.Length() - mLength;
    for( size_t frame = 1; frame <= maxFrame; frame++, ptr++ ) {
      kmer >>= BitMapPolicy< Alphabet >::NumBits;
      int8_t val = bitMapNucleotide( *ptr );
      if( val < 0 ) {
        lastAmbigIndex = frame + mLength - 1;
      } else {
        kmer |= ( val << bitIndex( mLength - 1 ) );
      }

      if( lastAmbigIndex == ( size_t ) -1 || frame > lastAmbigIndex ) {
        block( kmer, frame );
      } else {
        block( AmbiguousKmer, frame );
      }
    }
  }

  size_t Count() const {
    return mRef.Length() - mLength + 1;
  }

private:
  size_t                      mLength;
  const Sequence< Alphabet >& mRef;
};
