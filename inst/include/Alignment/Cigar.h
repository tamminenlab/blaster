#pragma once

#include "../Utils.h"

#include <deque>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <vector>
#include <algorithm>

enum class CigarOp : char {
  Unknown   = ' ',
  Match     = '=',
  Mismatch  = 'X',
  Deletion  = 'D',
  Insertion = 'I',
};
using CigarOps = std::vector< CigarOp >;

class CigarEntry {
public:
  int     count = 0;
  CigarOp op    = CigarOp::Unknown;

  CigarEntry() {}
  CigarEntry( const int count, const CigarOp op ) : count( count ), op( op ) {}

  bool operator==( const CigarEntry& other ) const {
    return count == other.count && op == other.op;
  }
};

class Cigar : public std::deque< CigarEntry > {
public:
  Cigar() {}

  Cigar( const char* str ) : Cigar( std::string( str ) ) {}

  Cigar( const std::string& str ) {
    // Separate "3M11C" into "3 M 11 M";
    std::string sep;
    bool        lastNumeric = true;
    for( auto ch : str ) {
      bool isNumeric = ch >= '0' && ch <= '9';
      if( isNumeric != lastNumeric ) {
        sep += ' ';
      }
      sep += ch;
      lastNumeric = isNumeric;
    }

    // Now parse separated str
    std::istringstream iss( sep );
    char               op;
    int                count;

    while( 1 ) {
      iss >> count;
      iss >> op;

      if( iss.eof() )
        break;

      Add( { count, ( CigarOp ) op } );
    }
  }

  Cigar operator+( const Cigar& other ) const {
    Cigar ce = *this;
    for( auto& c : other )
      ce.Add( c );
    return ce;
  }

  Cigar& operator+=( const Cigar& other ) {
    for( auto& c : other )
      Add( c );
    return *this;
  }

  void Clear() {
    clear();
  }

  void Reverse() {
    std::reverse( begin(), end() );
  }

  void Add( const CigarOp& op ) {
    Add( CigarEntry( 1, op ) );
  }

  void Add( const CigarEntry& entry ) {
    if( entry.count == 0 )
      return;

    if( entry.op == CigarOp::Unknown )
      return;

    if( empty() ) {
      push_back( entry );
    } else {
      auto& last = *rbegin();
      if( last.op == entry.op ) {
        // merge
        last.count += entry.count;
      } else {
        push_back( entry );
      }
    }
  }

  float Identity() const {
    size_t cols    = 0;
    size_t matches = 0;

    for( const CigarEntry& c : *this ) {
      // Don't count terminal gaps towards identity calculation
      if( &c == &( *this->cbegin() ) &&
          ( c.op == CigarOp::Insertion || c.op == CigarOp::Deletion ) )
        continue;
      if( &c == &( *this->crbegin() ) &&
          ( c.op == CigarOp::Insertion || c.op == CigarOp::Deletion ) )
        continue;

      cols += c.count;
      if( c.op == CigarOp::Match )
        matches += c.count;
    }

    return cols > 0 ? float( matches ) / float( cols ) : 0.0f;
  }

  std::string ToString() const {
    std::stringstream ss;
    for( auto& c : *this ) {
      ss << c.count << ( char ) c.op;
    }
    return ss.str();
  }
};

static std::ostream& operator<<( std::ostream& os, const Cigar& cigar ) {
  return ( os << cigar.ToString() );
}
