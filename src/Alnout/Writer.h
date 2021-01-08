#pragma once

#include "../Database/HitWriter.h"

#include "../Alphabet/DNA.h"
#include "../Alphabet/Protein.h"

#include <fstream>

namespace Alnout {

template < typename Alphabet >
class Writer : public HitWriter< Alphabet > {
public:
  using HitWriter< Alphabet >::HitWriter;

  HitWriter< Alphabet >&
  operator<<( const QueryHitsPair< Alphabet >& queryWithHits ) {
    const auto& query = queryWithHits.first;
    const auto& hits  = queryWithHits.second;

    auto& out = this->mOutput;

    // Output with fixed precision (sticky)
    out << std::setiosflags( std::ios::fixed );

    out << "Query >" << query.identifier << std::endl;
    out << " %Id   TLen  Target" << std::endl;
    for( auto& hit : hits ) {
      out << std::setprecision( 0 ) << std::setw( 3 )
              << ( hit.alignment.Identity() * 100.0 ) << '%' << std::setw( 7 )
              << hit.target.Length() << "  " << hit.target.identifier
              << std::endl;
    }
    out << std::endl;

    for( const auto& hit : hits ) {
      auto queryLen  = std::to_string( query.Length() );
      auto targetLen = std::to_string( hit.target.Length() );
      auto maxLen    = std::max( queryLen.size(), targetLen.size() );

      out << " Query" << std::setw( maxLen + 1 )
              << std::to_string( query.Length() ) << Unit() << " >"
              << query.identifier << std::endl;
      out << "Target" << std::setw( maxLen + 1 )
              << std::to_string( hit.target.Length() ) << Unit() << " >"
              << hit.target.identifier << std::endl;

      size_t numCols, numMatches, numGaps;
      auto   lines = ExtractAlignmentLines(
        QueryForAlignment( hit, query ), TargetForAlignment( hit, query ), hit.alignment,
        &numCols, &numMatches, &numGaps );

      out << std::endl;
      for( auto& line : lines ) {
        auto padLen = std::max( { std::to_string( lines.back().qs ).size(),
                                  std::to_string( lines.back().ts ).size(),
                                  std::to_string( lines.back().qe ).size(),
                                  std::to_string( lines.back().te ).size() } );

        auto qstrand = QueryStrand( hit );
        auto tstrand = TargetStrand( hit );

        if( !qstrand.empty() )
          qstrand += " ";

        if( !tstrand.empty() )
          tstrand += " ";

        out << "Qry " << std::setw( padLen )
                << QueryPos( line.qs, query, hit )
                << " " << qstrand << line.q << " "
                << QueryPos( line.qe, query, hit ) << std::endl;

        out << std::string( 5 + padLen + qstrand.size(), ' ' ) << line.a
                << std::endl;

        out << "Tgt " << std::setw( padLen )
                << TargetPos( line.ts, query, hit )
                << " " << tstrand << line.t << " "
                << TargetPos( line.te, query, hit ) << std::endl;

        out << std::endl;
      }

      float identity  = float( numMatches ) / float( numCols );
      float gapsRatio = float( numGaps ) / float( numCols );
      out << numCols << " cols, " << numMatches << " ids ("
              << std::setprecision( 1 ) << ( 100.0f * identity ) << "%), "
              << numGaps << " gaps (" << std::setprecision( 1 )
              << ( 100.0f * gapsRatio ) << "%)" << std::endl;

      out << std::endl;
    }
    return *this;
  }

private:
  static const size_t MaxLineLength = 60;

  static inline char        MatchSymbol( const char A, const char B );
  static inline std::string Unit();

  static inline std::string QueryStrand( const Hit< Alphabet >& hit ) {
    return "";
  }
  static inline std::string TargetStrand( const Hit< Alphabet >& hit ) {
    return "";
  }

  static inline Sequence< Alphabet >
  QueryForAlignment( const Hit< Alphabet >&      hit,
                     const Sequence< Alphabet >& query ) {
    return query;
  }

  static inline Sequence< Alphabet >
  TargetForAlignment( const Hit< Alphabet >&      hit,
                      const Sequence< Alphabet >& query ) {
    return hit.target;
  }

  static inline size_t QueryPos( const size_t                pos,
                                 const Sequence< Alphabet >& query,
                                 const Hit< Alphabet >&      hit ) {
    return pos;
  }

  static inline size_t TargetPos( const size_t                pos,
                                  const Sequence< Alphabet >& query,
                                  const Hit< Alphabet >&      hit ) {
    return pos;
  }

  using AlignmentLine = struct {
    size_t      qs, qe;
    std::string q;

    size_t      ts, te;
    std::string t;

    std::string a;
  };
  using AlignmentLines = std::deque< AlignmentLine >;

  static AlignmentLines ExtractAlignmentLines(
    const Sequence< Alphabet >& query, const Sequence< Alphabet >& target,
    const Cigar& alignment, size_t* outNumCols = NULL,
    size_t* outNumMatches = NULL, size_t* outNumGaps = NULL ) {
    Cigar cigar = alignment;

    size_t queryStart  = 0;
    size_t targetStart = 0;

    // Dont take left terminal gap into account
    if( !cigar.empty() ) {
      const auto& fce = cigar.front();
      if( fce.op == CigarOp::Deletion ) {
        targetStart = fce.count;
        cigar.pop_front();
      } else if( fce.op == CigarOp::Insertion ) {
        queryStart = fce.count;
        cigar.pop_front();
      }
    }

    // Don't take right terminal gap into account
    if( !cigar.empty() ) {
      const auto& bce = cigar.back();
      if( bce.op == CigarOp::Deletion ) {
        cigar.pop_back();
      } else if( bce.op == CigarOp::Insertion ) {
        cigar.pop_back();
      }
    }

    bool   match;
    size_t numMatches = 0;
    size_t numCols    = 0;
    size_t numGaps    = 0;

    size_t qcount = queryStart;
    size_t tcount = targetStart;

    AlignmentLine line;
    line.qs = queryStart + 1;
    line.ts = targetStart + 1;

    AlignmentLines lines;

    for( auto& c : cigar ) {
      for( int i = 0; i < c.count; i++ ) {
        switch( c.op ) {
          case CigarOp::Insertion:
            line.t += '-';
            line.q += query[ qcount++ ];
            line.a += ' ';
            numGaps++;
            break;

          case CigarOp::Deletion:
            line.q += '-';
            line.t += target[ tcount++ ];
            line.a += ' ';
            numGaps++;
            break;

          case CigarOp::Match: // specific for match
            numMatches++;
          case CigarOp::Mismatch: // match and mismatch
            line.q += query[ qcount++ ];
            line.t += target[ tcount++ ];
            {
              const char a = line.q.back(), b = line.t.back();
              line.a += MatchSymbol( a, b );
            }
            break;

          default:
            break;
        }

        numCols++;
        if( numCols % MaxLineLength == 0 ) {
          line.qe = qcount;
          line.te = tcount;
          lines.push_back( line );

          // Start new line
          line    = AlignmentLine();
          line.qs = qcount + 1;
          line.ts = tcount + 1;
        }
      }
    }

    if( !line.a.empty() ) {
      line.qe = qcount;
      line.te = tcount;
      lines.push_back( line );
    }

    if( outNumCols )
      *outNumCols = numCols;

    if( outNumMatches )
      *outNumMatches = numMatches;

    if( outNumGaps )
      *outNumGaps = numGaps;

    return lines;
  }

};

// DNA specializations
template <>
inline char Writer< DNA >::MatchSymbol( const char A, const char B ) {
  if( A == B ) {
    return '|';
  }
  if( MatchPolicy< DNA >::Match( A, B ) ) {
    return '+';
  }
  return ' ';
}

template <>
inline std::string Writer< DNA >::Unit() {
  return "nt";
}

template <>
inline std::string Writer< DNA >::QueryStrand( const Hit< DNA >& hit ) {
  return hit.strand == DNA::Strand::Minus ? "-" : "+";
}

template <>
inline std::string Writer< DNA >::TargetStrand( const Hit< DNA >& hit ) {
  // target is always reference (plus strand)
  return "+";
}

template <>
inline size_t Writer< DNA >::QueryPos( const size_t           pos,
                                       const Sequence< DNA >& query,
                                       const Hit< DNA >&      hit ) {
  return hit.strand == DNA::Strand::Minus ? ( query.Length() - pos + 1 ) : pos;
}

template <>
inline size_t Writer< DNA >::TargetPos( const size_t           pos,
                                        const Sequence< DNA >& query,
                                        const Hit< DNA >&      hit ) {
  // target is always reference (plus strand)
  return pos;
}

template <>
inline Sequence< DNA >
Writer< DNA >::QueryForAlignment( const Hit< DNA >&      hit,
                                  const Sequence< DNA >& query ) {
  return hit.strand == DNA::Strand::Minus ? query.Reverse().Complement()
                                          : query;
}

template <>
inline Sequence< DNA >
Writer< DNA >::TargetForAlignment( const Hit< DNA >&      hit,
                                   const Sequence< DNA >& query ) {
  // target is always reference (plus strand)
  return hit.target;
}

// Protein specializations
template <>
inline char Writer< Protein >::MatchSymbol( const char A, const char B ) {
  auto score = ScorePolicy< Protein >::Score( A, B );
  if( score >= 4 )
    return '|';
  if( score >= 2 )
    return ':';
  if( score > 0 )
    return '.';
  return ' ';
}

template <>
inline std::string Writer< Protein >::Unit() {
  return "aa";
}

} // namespace Alnout
