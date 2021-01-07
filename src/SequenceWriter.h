#pragma once

#include "Sequence.h"

#include <fstream>

template< typename Alphabet >
class SequenceWriter  {
public:
  SequenceWriter( std::ostream& output ) : mOutput( output ) {}
  SequenceWriter( const std::string& pathToFile )
      : mFile( pathToFile ), mOutput( mFile ) {}

  virtual SequenceWriter< Alphabet >&
  operator<<( const Sequence< Alphabet >& seq ) = 0;

  virtual ~SequenceWriter() = default;

protected:
  std::ofstream mFile;
  std::ostream& mOutput;
};
