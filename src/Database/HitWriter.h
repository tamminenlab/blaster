#pragma once

#include "Search.h"

#include <fstream>

template < typename Alphabet >
class HitWriter {
public:
  HitWriter( std::ostream& output ) : mOutput( output ) {}
  HitWriter( const std::string& pathToFile )
      : mFile( pathToFile ), mOutput( mFile ) {}

  virtual HitWriter< Alphabet >&
  operator<<( const QueryHitsPair< Alphabet >& queryWithHits ) = 0;

  virtual ~HitWriter() = default;

protected:
  std::ofstream mFile;
  std::ostream& mOutput;
};
