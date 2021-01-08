#pragma once

template < typename Alphabet >
struct BitMapPolicy {
  static const size_t NumBits = 0;
  inline static int8_t BitMap( const char ch ) {
    return -1;
  }
};

template < typename Alphabet >
struct ComplementPolicy {
  inline static char Complement( const char ch ) {
    return ch;
  }
};

template < typename Alphabet >
struct MatchPolicy {
  inline static bool Match( const char chA, const char chB ) {
    return chA == chB;
  }
};

template < typename Alphabet >
struct ScorePolicy {
  inline static int8_t Score( const char chA, const char chB ) {
    return MatchPolicy< Alphabet >::Match( chA, chB ) ? 1 : -1;
  }
};
