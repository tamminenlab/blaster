#pragma once

#include <deque>
#include <vector>

#include "Sequence.h"
#include "Utils.h"

#include "Database/HSP.h"
#include "Database/Highscore.h"
#include "Database/Kmers.h"

#include "Alphabet.h"

using SequenceId = uint32_t; // SequenceId

template < typename Alphabet >
class Database {
public:
  enum ProgressType { StatsCollection, Indexing };
  using OnProgressCallback =
    std::function< void( ProgressType, const size_t, const size_t ) >;

  Database( const size_t kmerLength );

  void SetProgressCallback( const OnProgressCallback& progressCallback );
  void Initialize( const SequenceList< Alphabet >& sequences );

  size_t NumSequences() const;
  size_t KmerLength() const;
  size_t MaxUniqueKmers() const;

  const Sequence< Alphabet >& GetSequenceById( const SequenceId& seqId ) const;

  bool GetKmersForSequenceId( const SequenceId& seqId, const Kmer** kmers,
                              size_t* numKmers ) const;
  bool GetSequenceIdsIncludingKmer( const Kmer& kmer, const SequenceId** seqIds,
                                    size_t* numSeqIds ) const;

private:
  OnProgressCallback mProgressCallback;
  SequenceList< Alphabet > mSequences;

  std::vector< Kmer >   mKmers;

  size_t mKmerLength;
  size_t mMaxUniqueKmers;

  std::vector< SequenceId > mSequenceIds;
  std::vector< size_t >     mSequenceIdsOffsetByKmer;
  std::vector< size_t >     mSequenceIdsCountByKmer;

  std::vector< size_t > mKmerOffsetBySequenceId;
  std::vector< size_t > mKmerCountBySequenceId;

};

/*
 * Implementation
 */
template < typename A >
Database< A >::Database( const size_t kmerLength )
  :  mProgressCallback( []( ProgressType, const size_t, const size_t ) {} ),
     mKmerLength( kmerLength ),
     mMaxUniqueKmers( 1 << ( BitMapPolicy< A >::NumBits * mKmerLength ) )
{
  assert( BitMapPolicy< A >::NumBits * mKmerLength <= sizeof( Kmer ) * 8 );
}

template < typename A >
void Database< A >::SetProgressCallback(
  const OnProgressCallback& progressCallback ) {
  mProgressCallback = progressCallback;
}

template < typename A >
void Database< A >::Initialize( const SequenceList< A >& sequences ) {
  mSequences = sequences;

  size_t totalEntries       = 0;
  size_t totalUniqueEntries = 0;

  /* std::vector< uint32_t > count( mMaxUniqueKmers ); */
  std::vector< size_t >     uniqueCount( mMaxUniqueKmers );
  std::vector< SequenceId > uniqueIndex( mMaxUniqueKmers, -1 );

  for( SequenceId seqId = 0; seqId < mSequences.size(); seqId++ ) {
    const Sequence< A >& seq = mSequences[ seqId ];

    Kmers< A > kmers( seq, mKmerLength );
    kmers.ForEach( [&]( const Kmer kmer, const size_t pos ) {
      totalEntries++;

      // Count unique words
      if( kmer == AmbiguousKmer || uniqueIndex[ kmer ] == seqId )
        return;

      uniqueIndex[ kmer ] = seqId;
      uniqueCount[ kmer ]++;
      totalUniqueEntries++;
    } );

    // Progress
    if( seqId % 512 == 0 || seqId + 1 == mSequences.size() ) {
      mProgressCallback( ProgressType::StatsCollection, seqId + 1,
                         mSequences.size() );
    }
  }

  // Calculate indices
  mSequenceIdsOffsetByKmer.reserve( mMaxUniqueKmers );
  for( size_t i = 0; i < mMaxUniqueKmers; i++ ) {
    mSequenceIdsOffsetByKmer[ i ] =
      i > 0 ? mSequenceIdsOffsetByKmer[ i - 1 ] + uniqueCount[ i - 1 ] : 0;
  }

  // Populate DB
  mSequenceIds.reserve( totalUniqueEntries );
  mKmers.reserve( totalEntries );

  // Reset to 0
  mSequenceIdsCountByKmer = std::vector< size_t >( mMaxUniqueKmers );
  mKmerCountBySequenceId  = std::vector< size_t >( mSequences.size() );
  mKmerOffsetBySequenceId = std::vector< size_t >( mSequences.size() );

  uniqueIndex = std::vector< SequenceId >( mMaxUniqueKmers, -1 );

  auto   kmersData = mKmers.data();
  size_t kmerCount = 0;

  for( SequenceId seqId = 0; seqId < mSequences.size(); seqId++ ) {
    const Sequence< A >& seq = mSequences[ seqId ];

    mKmerOffsetBySequenceId[ seqId ] = kmerCount;

    Kmers< A > kmers( seq, mKmerLength );
    kmers.ForEach( [&]( const Kmer kmer, const size_t pos ) {
      // Encode position in kmersData implicitly
      // by saving _every_ kmer
      kmersData[ kmerCount++ ] = kmer;

      if( kmer == AmbiguousKmer || uniqueIndex[ kmer ] == seqId )
        return;

      uniqueIndex[ kmer ] = seqId;

      mSequenceIds[ mSequenceIdsOffsetByKmer[ kmer ] +
                    mSequenceIdsCountByKmer[ kmer ] ] = seqId;
      mSequenceIdsCountByKmer[ kmer ]++;
    } );

    mKmerCountBySequenceId[ seqId ] =
      kmerCount - mKmerOffsetBySequenceId[ seqId ];

    // Progress
    if( seqId % 512 == 0 || seqId + 1 == mSequences.size() ) {
      mProgressCallback( ProgressType::Indexing, seqId + 1, mSequences.size() );
    }
  }
}

template < typename A >
const Sequence< A >&
Database< A >::GetSequenceById( const SequenceId& seqId ) const {
  assert( seqId < NumSequences() );
  return mSequences[ seqId ];
}

template < typename A >
size_t Database< A >::NumSequences() const {
  return mSequences.size();
}

template < typename A >
size_t Database< A >::MaxUniqueKmers() const {
  return mMaxUniqueKmers;
}

template < typename A >
size_t Database< A >::KmerLength() const {
  return mKmerLength;
}

template < typename A >
bool Database< A >::GetKmersForSequenceId( const SequenceId& seqId,
                                           const Kmer**      kmers,
                                           size_t*           numKmers ) const {
  if( seqId >= NumSequences() )
    return false;

  const auto& offset = mKmerOffsetBySequenceId[ seqId ];
  const auto& count  = mKmerCountBySequenceId[ seqId ];

  *kmers    = &mKmers[ offset ];
  *numKmers = count;
  return count > 0;
}

template < typename A >
bool Database< A >::GetSequenceIdsIncludingKmer( const Kmer&        kmer,
                                                 const SequenceId** seqIds,
                                                 size_t* numSeqIds ) const {
  if( kmer == AmbiguousKmer )
    return false;

  if( kmer >= MaxUniqueKmers() )
    return false;

  const auto& offset = mSequenceIdsOffsetByKmer[ kmer ];
  const auto& count  = mSequenceIdsCountByKmer[ kmer ];

  *seqIds    = &mSequenceIds[ offset ];
  *numSeqIds = count;
  return count > 0;
}
