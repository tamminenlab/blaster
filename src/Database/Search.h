#pragma once

#include "../Alignment/Cigar.h"
#include "../Database.h"
#include "../Sequence.h"

#include "../Alphabet/DNA.h"

#include <deque>
#include <vector>

struct BaseSearchParams {
  int   maxAccepts  = 1;
  int   maxRejects  = 16;
  float minIdentity = 0.75f;
};

template < typename Alphabet >
struct SearchParams : public BaseSearchParams {};

template <>
struct SearchParams< DNA > : public BaseSearchParams {
  DNA::Strand strand = DNA::Strand::Plus;
};

template < typename Alphabet >
struct Hit {
  Sequence< Alphabet > target;
  Cigar                alignment;
};

template <>
struct Hit< DNA > {
  Sequence< DNA > target;
  Cigar           alignment;
  DNA::Strand     strand;
};

template < typename Alphabet >
using HitList = std::deque< Hit< Alphabet > >;

template < typename Alphabet >
using QueryHitsPair = std::pair< Sequence< Alphabet >, HitList< Alphabet > >;

template < typename Alphabet >
using SearchForHitsCallback =
  std::function< void( const Sequence< Alphabet >&, const Cigar& ) >;

template < typename Alphabet >
class Search {
public:
  Search( const Database< Alphabet >&     db,
          const SearchParams< Alphabet >& params )
      : mDB( db ), mParams( params ) {}

  inline HitList< Alphabet > Query( const Sequence< Alphabet >& query ) {
    HitList< Alphabet > hits;

    SearchForHits(
      query, [&]( const Sequence< Alphabet >& target, const Cigar& alignment ) {
        hits.push_back( { target, alignment } );
      } );

    return hits;
  }

protected:
  virtual void
  SearchForHits( const Sequence< Alphabet >&              query,
                 const SearchForHitsCallback< Alphabet >& callback ) = 0;

  const Database< Alphabet >&     mDB;
  const SearchParams< Alphabet >& mParams;
};

/*
 * For DNA, allow strand specification
 */
template <>
inline HitList< DNA > Search< DNA >::Query( const Sequence< DNA >& query ) {
  HitList< DNA > hits;

  auto strand = mParams.strand;

  if( strand == DNA::Strand::Plus || strand == DNA::Strand::Both ) {
    SearchForHits(
      query, [&]( const Sequence< DNA >& target, const Cigar& alignment ) {
        hits.push_back( { target, alignment, DNA::Strand::Plus } );
      } );
  }

  if( strand == DNA::Strand::Minus || strand == DNA::Strand::Both ) {
    SearchForHits(
      query.Reverse().Complement(),
      [&]( const Sequence< DNA >& target, const Cigar& alignment ) {
        hits.push_back( { target, alignment, DNA::Strand::Minus } );
      } );
  }

  return hits;
}
