#pragma once

#include "../Database/HitWriter.h"

#include <fstream>
#include <array>

namespace CSV {

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

    // Each hit gets a line
    for( const auto& hit : hits ) {
      // Compute stuff
      // ----
      Cigar cigar = hit.alignment;

      size_t qs = 0, qe = query.Length() - 1;
      size_t ts = 0, te = hit.target.Length() - 1;

      // Dont take left terminal gap into account
      if( !cigar.empty() ) {
        const auto& fce = cigar.front();
        if( fce.op == CigarOp::Deletion ) {
          ts += fce.count;
          cigar.pop_front();
        } else if( fce.op == CigarOp::Insertion ) {
          qs += fce.count;
          cigar.pop_front();
        }
      }

      // Don't take right terminal gap into account
      if( !cigar.empty() ) {
        const auto& bce = cigar.back();
        if( bce.op == CigarOp::Deletion ) {
          te -= bce.count;
          cigar.pop_back();
        } else if( bce.op == CigarOp::Insertion ) {
          qe -= bce.count;
          cigar.pop_back();
        }
      }

      Sequence< Alphabet > targetMatchSeq = hit.target.Subsequence( ts, te - ts + 1 );
      Sequence< Alphabet > queryMatchSeq;
      if( IsHitOnOtherStrand( hit ) ) {
        // Minus strand -> Reverse complemented query has been hit
        // (Alignment refers to the reverse complemented query)
        queryMatchSeq =
          query.Reverse().Complement().Subsequence( qs, qe - qs + 1 );
        // Then encode this information in queryMatchStart and queryMatchEnd
        // (queryMatchStart > queryMatchEnd)
        qs = query.Length() - qs - 1;
        qe = query.Length() - qe - 1;
      } else {
         queryMatchSeq = query.Subsequence( qs, qe - qs + 1 );
      }

      size_t numMatches = 0, numMismatches = 0, numColumns = 0, numGaps = 0;
      for( auto& c : cigar ) {
        for( int i = 0; i < c.count; i++ ) {
          numColumns++;
          switch( c.op ) {
            case CigarOp::Insertion: numGaps++; break;
            case CigarOp::Deletion: numGaps++; break;
            case CigarOp::Match: numMatches++; break;
            case CigarOp::Mismatch: numMismatches++; break;
            default: break;
          }
        }
      }

      float identity  = float( numMatches ) / float( numColumns );

      // Write stuff
      // ----

      // QueryId
      out << EscapeStringForCSV( query.identifier ) << ",";

      // TargetId
      out << EscapeStringForCSV( hit.target.identifier ) << ",";

      // QueryMatchStart
      out << qs + 1 << ",";

      // QueryMatchEnd
      out << qe + 1 << ",";

      // TargetMatchStart
      out << ts + 1 << ",";

      // TargetMatchEnd
      out << te + 1 << ",";

      // QueryMatchSeq
      out << EscapeStringForCSV( queryMatchSeq.sequence ) << ",";

      // TargetMatchSeq
      out << EscapeStringForCSV( targetMatchSeq.sequence ) << ",";

      // NumColumns, NumMatches, NumMismatches, NumGaps
      out << numColumns << "," << numMatches << "," << numMismatches << ","
          << numGaps << ",";

      // Identity
      out << std::setprecision( 3 ) << identity << ",";

      // Alignment
      out << cigar.ToString();

      out << std::endl;
    }

    return *this;
  }

private:
  std::string EscapeStringForCSV( const std::string& value ) {
    std::string ret = value;

    static const std::array< char, 4 > escapingRequiredChars = { ',', '\"', '\r', '\n' };

    bool needsEscaping =
      value.empty() ||
      std::any_of(
        escapingRequiredChars.begin(), escapingRequiredChars.end(),
        [&]( const char ch ) { return value.find( ch ) != std::string::npos; } );

    if( needsEscaping ) {
      ret.insert( 0, 1, '\"' );
      ret += '\"';
    }

    return ret;
  }

  static inline bool IsHitOnOtherStrand( const Hit< Alphabet >& hit ) {
    return false;
  }
};

template <>
inline bool Writer< DNA >::IsHitOnOtherStrand( const Hit< DNA >& hit ) {
  return hit.strand == DNA::Strand::Minus;
}

} // namespace CSV
