#pragma once

#include "../Alphabet.h"

struct DNA {
  typedef char CharType;

  enum class Strand { Plus, Minus, Both };
};

using RNA = DNA;

template <>
struct BitMapPolicy< DNA > {
  static const size_t NumBits = 2;

  inline static int8_t BitMap( const char base ) {
    switch( base ) {
      case 'A': return 0b00;
      case 'C': return 0b01;
      case 'U':
      case 'T': return 0b10;
      case 'G': return 0b11;
      default: return -1; // ambiguity
    }
  }
};

template <>
struct ComplementPolicy< DNA > {
  inline static char Complement( const char nuc ) {
    switch( nuc ) {
      case 'A': return 'T';
      case 'G': return 'C';
      case 'C': return 'G';
      case 'T': return 'A';
      case 'U': return 'A';

      case 'Y': return 'R';
      case 'R': return 'Y';
      case 'W': return 'W';
      case 'S': return 'S';
      case 'K': return 'M';
      case 'M': return 'K';

      case 'D': return 'H';
      case 'V': return 'B';
      case 'H': return 'D';
      case 'B': return 'V';
      case 'N': return 'N';
    }

    return nuc;
  }
};

template <>
struct ScorePolicy< DNA > {
  inline static int8_t Score( const char nucA, const char nucB ) {
    static const int ScoreMatrixSize = 26; // 'A'...'Z', waste some space for faster lookup
    static const int8_t ScoreMatrix[ ScoreMatrixSize ][ ScoreMatrixSize ] = {
      {  2, -4, -4,  2,  0,  0, -4,  2,  0,  0, -4,  0,  2,  2,  0,  0,  0,  2, -4, -4, -4,  2,  2,  0, -4,  0 },
      { -4,  2,  2,  2,  0,  0,  2,  2,  0,  0,  2,  0,  2,  2,  0,  0,  0,  2,  2,  2,  2,  2,  2,  0,  2,  0 },
      { -4,  2,  2, -4,  0,  0, -4,  2,  0,  0, -4,  0,  2,  2,  0,  0,  0, -4,  2, -4, -4,  2, -4,  0,  2,  0 },
      {  2,  2, -4,  2,  0,  0,  2,  2,  0,  0,  2,  0,  2,  2,  0,  0,  0,  2,  2,  2,  2,  2,  2,  0,  2,  0 },
      {  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 },
      {  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 },
      { -4,  2, -4,  2,  0,  0,  2, -4,  0,  0,  2,  0, -4,  2,  0,  0,  0,  2,  2, -4, -4,  2, -4,  0, -4,  0 },
      {  2,  2,  2,  2,  0,  0, -4,  2,  0,  0,  2,  0,  2,  2,  0,  0,  0,  2,  2,  2,  2,  2,  2,  0,  2,  0 },
      {  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 },
      {  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 },
      { -4,  2, -4,  2,  0,  0,  2,  2,  0,  0,  2,  0, -4,  2,  0,  0,  0,  2,  2,  2,  2,  2,  2,  0,  2,  0 },
      {  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 },
      {  2,  2,  2,  2,  0,  0, -4,  2,  0,  0, -4,  0,  2,  2,  0,  0,  0,  2,  2, -4, -4,  2,  2,  0,  2,  0 },
      {  2,  2,  2,  2,  0,  0,  2,  2,  0,  0,  2,  0,  2,  2,  0,  0,  0,  2,  2,  2,  2,  2,  2,  0,  2,  0 },
      {  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 },
      {  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 },
      {  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 },
      {  2,  2, -4,  2,  0,  0,  2,  2,  0,  0,  2,  0,  2,  2,  0,  0,  0,  2,  2, -4, -4,  2,  2,  0, -4,  0 },
      { -4,  2,  2,  2,  0,  0,  2,  2,  0,  0,  2,  0,  2,  2,  0,  0,  0,  2,  2, -4, -4,  2, -4,  0,  2,  0 },
      { -4,  2, -4,  2,  0,  0, -4,  2,  0,  0,  2,  0, -4,  2,  0,  0,  0, -4, -4,  2,  2, -4,  2,  0,  2,  0 },
      { -4,  2, -4,  2,  0,  0, -4,  2,  0,  0,  2,  0, -4,  2,  0,  0,  0, -4, -4,  2,  2, -4,  2,  0,  2,  0 },
      {  2,  2,  2,  2,  0,  0,  2,  2,  0,  0,  2,  0,  2,  2,  0,  0,  0,  2,  2, -4, -4,  2,  2,  0,  2,  0 },
      {  2,  2, -4,  2,  0,  0, -4,  2,  0,  0,  2,  0,  2,  2,  0,  0,  0,  2, -4,  2,  2,  2,  2,  0,  2,  0 },
      {  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 },
      { -4,  2,  2,  2,  0,  0, -4,  2,  0,  0,  2,  0,  2,  2,  0,  0,  0, -4,  2,  2,  2,  2,  2,  0,  2,  0 },
      {  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 }
    };

    return ScoreMatrix[ nucA - 'A' ][ nucB - 'A' ];
  }
};

template <>
struct MatchPolicy < DNA > {
  inline static bool Match( const char nucA, const char nucB ) {
    return ScorePolicy< DNA >::Score( nucA, nucB ) > 0;
  }
};
