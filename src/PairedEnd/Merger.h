#pragma once

#include "../Sequence.h"
#include "../FASTQ/QScore.h"
#include "../Utils.h"

#include <cassert>
#include <cfloat>
#include <iostream>

namespace PairedEnd {
// TODO: Find good defaults! Seem to stringent
static const int    MERGER_DEFAULT_MIN_OVERLAP  = 16; // bases
static const double MERGER_DEFAULT_MIN_IDENTITY = 0.9;

template < typename Alphabet >
class Merger {
public:
  Merger( const int    minOverlap  = MERGER_DEFAULT_MIN_OVERLAP,
          const double minIdentity = MERGER_DEFAULT_MIN_IDENTITY );
  bool Merge( const Sequence< Alphabet >& fwd, const Sequence< Alphabet >& rev,
              Sequence< Alphabet >* merged ) const;

private:
  int    mMinOverlap;
  double mMinIdentity;

  typedef struct {
    size_t length;
    size_t pos1;
    size_t pos2;
  } OverlapInfo;

  double ComputeOverlapScore( const char* sequence1, const char* sequence2,
                              const char* quality1, const char* quality2,
                              const size_t len ) const;
  bool FindBestOverlap( const Sequence< Alphabet >& sequence1,
                        const Sequence< Alphabet >& sequence2,
                        OverlapInfo*                overlap ) const;
  bool IsStaggered( const OverlapInfo& overlap ) const;
  void PrintOverlap( const Sequence< Alphabet >& seq1,
                     const Sequence< Alphabet >& seq2,
                     const OverlapInfo&          overlap ) const;
};

/*
 * Implementation
 */
template < typename A >
Merger< A >::Merger( const int minOverlap, const double minIdentity )
    : mMinOverlap( minOverlap ), mMinIdentity( minIdentity ) {}

template < typename A >
bool Merger< A >::Merge( const Sequence< A >& fwd, const Sequence< A >& rev,
                         Sequence< A >* mergedOut ) const {

  using FASTQ::QScore;

  OverlapInfo overlap;

  assert( fwd.quality.length() == fwd.sequence.length() );
  assert( rev.quality.length() == rev.sequence.length() );

  Sequence< A > seq1 = fwd;
  Sequence< A > seq2 = rev.Reverse().Complement();

  if( !FindBestOverlap( seq1, seq2, &overlap ) )
    return false;

  Sequence< A > merged = seq1.Subsequence( 0, overlap.length );

  for( int i = 0; i < overlap.length; i++ ) {
    char s1 = seq1.sequence[ overlap.pos1 + i ];
    char s2 = seq2.sequence[ overlap.pos2 + i ];

    int q1 = seq1.quality[ overlap.pos1 + i ] - FASTQ::Q_MIN_ASCII_BASE;
    int q2 = seq2.quality[ overlap.pos2 + i ] - FASTQ::Q_MIN_ASCII_BASE;

    if( q1 >= q2 ) {
      // Call X as merged base
      merged.sequence[ i ] = s1;
    } else {
      // Call Y as merged base
      merged.sequence[ i ] = s2;
    }

    if( MatchPolicy< A >::Match( s1, s2 ) ) {
      merged.quality[ i ] =
        FASTQ::Q_MIN_ASCII_BASE +
        QScore::Instance().CalculatePosteriorScoreForMatch( q1, q2 );
    } else {
      merged.quality[ i ] =
        FASTQ::Q_MIN_ASCII_BASE +
        QScore::Instance().CalculatePosteriorScoreForMismatch( q1, q2 );
    }
  }

  // AAA
  //  BBB -> AMMB (nonstaggered)
  //
  //  AAA
  // BBB  -> MM (staggered)
  if( !IsStaggered( overlap ) ) {
    Sequence< A > leftOverhang  = seq1.Subsequence( 0, overlap.pos1 );
    Sequence< A > rightOverhang = seq2.Subsequence( overlap.length );

    merged = leftOverhang + merged + rightOverhang;
  }

  /* std::cout << seq1.identifier << std::endl; */
  /* PrintOverlap( seq1, seq2, overlap ); */

  /* std::cout << merged.sequence << std::endl; */
  /* std::cout << merged.quality << std::endl; */
  /* std::cout << string( merged.sequence.length(), '-' ) << std::endl; */

  *mergedOut = std::move( merged );
  return true;
}

template < typename A >
double Merger< A >::ComputeOverlapScore( const char*  sequence1,
                                         const char*  sequence2,
                                         const char*  quality1,
                                         const char*  quality2,
                                         const size_t len ) const {
  using FASTQ::QScore;

  double score = 0.0;

  size_t numMismatches = 0;
  size_t maxMismatches = len - size_t( len * mMinIdentity );

  for( int i = 0; i < len; i++ ) {
    if( MatchPolicy< A >::Match( sequence1[ i ], sequence2[ i ] ) ) {
      score +=
        ( 1.0 - QScore::Instance().CalculatePosteriorErrorProbabilityForMatch(
                  quality1[ i ], quality2[ i ] ) );
    } else {
      score -=
        ( 1.0 -
          QScore::Instance().CalculatePosteriorErrorProbabilityForMismatch(
            quality1[ i ], quality2[ i ] ) );

      numMismatches++;
      if( numMismatches > maxMismatches ) {
        score = DBL_MIN;
        break;
      }
    }
  }

  return score;
}

/*
 * Slide REV along FWD, start at the right to the left
 * We care about staggered reads, so don't stop at i=3
 *
 * i = 0
 * AAA
 *    BB
 *
 * i = 1
 * AAA
 *   BB
 *
 * ...
 *
 * i = 3
 * AAA
 * BB
 *
 * ...
 *
 * i = 5 AAA
 * BB
 *
 */
template < typename A >
bool Merger< A >::FindBestOverlap( const Sequence< A >& seq1,
                                   const Sequence< A >& seq2,
                                   OverlapInfo*         overlap ) const {
  int len1 = seq1.Length();
  int len2 = seq2.Length();

  overlap->length = 0;
  overlap->pos1   = 0;
  overlap->pos2   = 0;

  double bestScore = DBL_MIN;

  // Slide base by base, finding best overlap
  for( int i = 0; i <= len1 + len2; i++ ) {
    int pos1 = std::max( len1 - i, 0 );
    int pos2 = std::max( i - len1, 0 );

    int length = std::min( len2 - pos2, i );
    if( length < mMinOverlap )
      continue;

    double score = ComputeOverlapScore(
      seq1.sequence.c_str() + pos1, seq2.sequence.c_str() + pos2,
      seq1.quality.c_str() + pos1, seq2.quality.c_str() + pos2, length );

    if( score > bestScore ) {
      bestScore       = score;
      overlap->length = length;
      overlap->pos1   = pos1;
      overlap->pos2   = pos2;
    }
  }

  return ( bestScore > DBL_MIN );
}

template < typename A >
bool Merger< A >::IsStaggered( const OverlapInfo& overlap ) const {
  return overlap.pos2 > 0;
}

template < typename A >
void Merger< A >::PrintOverlap( const Sequence< A >& seq1,
                                const Sequence< A >& seq2,
                                const OverlapInfo&   overlap ) const {
  std::cout << std::string( overlap.pos2, ' ' ) << seq1.quality << std::endl;
  std::cout << std::string( overlap.pos2, ' ' ) << seq1.sequence << std::endl;
  std::cout << std::string( std::max( overlap.pos1, overlap.pos2 ), ' ' );
  for( int i = 0; i < overlap.length; i++ ) {
    if( seq1.sequence[ overlap.pos1 + i ] ==
        seq2.sequence[ overlap.pos2 + i ] ) {
      std::cout << '|';
    } else {
      std::cout << ' ';
    }
  }
  std::cout << std::endl;
  std::cout << std::string( overlap.pos1, ' ' ) << seq2.sequence << std::endl;
  std::cout << std::string( overlap.pos1, ' ' ) << seq2.quality << std::endl;
}

} // namespace PairedEnd
