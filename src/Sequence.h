#pragma once

#include <cassert>
#include <deque>
#include <iostream>
#include <string>
#include <algorithm>

#include "Utils.h"

#include <cassert>
#include <ctype.h>
#include <iostream>

#include "Alphabet.h"
#include "FASTQ/QScore.h"

template < typename Alphabet >
class Sequence {
public:
  Sequence();
  Sequence( const Sequence< Alphabet >& sequence );
  Sequence( Sequence< Alphabet >&& sequence );
  Sequence< Alphabet >& operator=( const Sequence< Alphabet >& other );
  Sequence( const std::string& sequence );
  Sequence( const char* sequence );
  Sequence( const std::string&                                      identifier,
            const std::basic_string< typename Alphabet::CharType >& sequence );
  Sequence( const std::string&                                      identifier,
            const std::basic_string< typename Alphabet::CharType >& sequence,
            const std::string&                                      quality );

  size_t Length() const;

  Sequence< Alphabet >
  Subsequence( const size_t pos,
               const size_t len = std::string::npos ) const;

  Sequence< Alphabet > operator+( const Sequence< Alphabet >& other ) const;
  bool                 operator==( const Sequence< Alphabet >& other ) const;
  bool                 operator!=( const Sequence< Alphabet >& other ) const;

  // Inline for faster lookup
  inline char& operator[]( const size_t index ) {
    assert( index >= 0 && index < sequence.size() );
    return sequence[ index ];
  }

  inline char operator[]( const size_t index ) const {
    assert( index >= 0 && index < sequence.size() );
    return sequence[ index ];
  }

  Sequence< Alphabet > Complement() const;
  Sequence< Alphabet > Reverse() const;

  float NumExpectedErrors() const;

  std::string identifier;
  std::string quality;

  std::basic_string< typename Alphabet::CharType > sequence;
};

template < typename Alphabet >
static std::ostream& operator<<( std::ostream&               os,
                                 const Sequence< Alphabet >& seq ) {
  if( !seq.identifier.empty() )
    os << ">" << seq.identifier << std::endl;
  if( !seq.sequence.empty() )
    os << " " << seq.sequence << std::endl;
  if( !seq.quality.empty() )
    os << " " << seq.quality << std::endl;
  return os;
}

template < typename Alphabet >
using SequenceList = std::deque< Sequence< Alphabet > >;

/*
 * Implementation
 */
template < typename A >
Sequence< A >::Sequence() : Sequence( "", "", "" ) {}

template < typename A >
Sequence< A >::Sequence( const Sequence& sequence )
    : sequence( sequence.sequence ), identifier( sequence.identifier ),
      quality( sequence.quality ) {}

template < typename A >
Sequence< A >::Sequence( Sequence< A >&& sequence )
    : sequence( std::move( sequence.sequence ) ),
      identifier( std::move( sequence.identifier ) ),
      quality( std::move( sequence.quality ) ) {}

template < typename A >
Sequence< A >& Sequence< A >::operator=( const Sequence< A >& other ) {
  sequence   = other.sequence;
  identifier = other.identifier;
  quality    = other.quality;
  return *this;
}

template < typename A >
Sequence< A >::Sequence( const std::string& sequence )
    : Sequence( "", sequence, "" ) {}

template < typename A >
Sequence< A >::Sequence( const char* sequence )
    : Sequence( "", sequence, "" ) {}

template < typename A >
Sequence< A >::Sequence(
  const std::string&                               identifier,
  const std::basic_string< typename A::CharType >& sequence )
    : Sequence( identifier, sequence, "" ) {}

template < typename A >
Sequence< A >::Sequence(
  const std::string&                               identifier,
  const std::basic_string< typename A::CharType >& sequence,
  const std::string&                               quality )
    : identifier( identifier ), sequence( sequence ), quality( quality ) {}

template < typename A >
size_t Sequence< A >::Length() const {
  return sequence.length();
}

template < typename A >
Sequence< A > Sequence< A >::Subsequence( const size_t pos,
                                          const size_t len_ ) const {
  size_t len = len_;
  if( len == std::string::npos ) {
    len = Length() - pos;
  }

  return Sequence( identifier,
                   pos < sequence.length() ? sequence.substr( pos, len ) : "",
                   pos < quality.length() ? quality.substr( pos, len ) : "" );
}

template < typename A >
Sequence< A > Sequence< A >::operator+( const Sequence< A >& other ) const {
  return Sequence< A >( identifier, sequence + other.sequence,
                        quality + other.quality );
}

template < typename A >
bool Sequence< A >::operator==( const Sequence< A >& other ) const {
  return !( *this != other );
}

template < typename A >
bool Sequence< A >::operator!=( const Sequence< A >& other ) const {
  if( Length() != other.Length() )
    return true;

  auto tit = ( *this ).sequence.begin();
  auto oit = other.sequence.begin();
  while( tit != ( *this ).sequence.end() && oit != other.sequence.end() ) {
    if( !MatchPolicy< A >::Match( *tit, *oit ) )
      return true;

    ++tit;
    ++oit;
  }

  return false;
}

template < typename A >
Sequence< A > Sequence< A >::Reverse() const {
  Sequence rev = *this;
  // Reverse sequence and quality
  std::reverse( rev.sequence.begin(), rev.sequence.end() );
  std::reverse( rev.quality.begin(), rev.quality.end() );
  return rev;
}

template < typename A >
Sequence< A > Sequence< A >::Complement() const {
  Sequence complement = *this;

  for( char& ch : complement.sequence ) {
    ch = ComplementPolicy< A >::Complement( ch );
  }

  return complement;
}

template < typename A >
float Sequence< A >::NumExpectedErrors() const {
  if( quality.empty() )
    return 0.0f;

  float numExpectedErrors = 0.0f;
  for( auto &q : quality ) {
    numExpectedErrors += FASTQ::QScore::Instance().AsciiToProbability( q );
  }
  return numExpectedErrors;
}
