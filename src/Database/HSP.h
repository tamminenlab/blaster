#pragma once

#include "../Alignment/Cigar.h"

#include <cassert>
#include <cmath>

// High-scoring segment pair
// HSP: first and last character in sequence (i.e. seq[a1] - seq[a2])
class HSP {
public:
  size_t a1, a2;
  size_t b1, b2;
  int    score;
  Cigar  cigar;

  HSP( const size_t a1, const size_t a2, const size_t b1, const size_t b2,
       const int score = 0 )
      : a1( a1 ), a2( a2 ), b1( b1 ), b2( b2 ), score( score ) {
    assert( a2 >= a1 && b2 >= b1 );
  }

  size_t Length() const {
    return std::max( a2 - a1, b2 - b1 ) + 1;
  }

  int Score() const {
    return score;
  }

  bool IsOverlapping( const HSP& other ) const {
    return ( a1 <= other.a2 && other.a1 <= a2 )     // overlap in A direction
           || ( b1 <= other.b2 && other.b1 <= b2 ); // overlap in B direction
  }

  size_t DistanceTo( const HSP& other ) const {
    size_t dx = ( a1 > other.a2 ? a1 - other.a2 : other.a1 - a2 );
    dx        = dx > 0 ? dx - 1 : 0;

    size_t dy = ( b1 > other.b2 ? b1 - other.b2 : other.b1 - b2 );
    dy        = dy > 0 ? dy - 1 : 0;

    return sqrt( dx * dx + dy * dy );
  }

  bool operator<( const HSP& other ) const {
    return Score() < other.Score();
  }
};
