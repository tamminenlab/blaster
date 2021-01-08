#pragma once

#include "Sequence.h"
#include "TextReader.h"
#include "Utils.h"

#include <memory>

template< typename Alphabet >
class SequenceReader {
public:
  SequenceReader( const std::string& pathToFile )
      : mTextReader( new TextFileReader( pathToFile ) ) {}

  SequenceReader( std::istream& is )
      : mTextReader( new TextStreamReader( is ) ) {}

  bool EndOfFile() const {
    return mTextReader->EndOfFile();
  }

  size_t NumBytesRead() const {
    return mTextReader->NumBytesRead();
  }

  size_t NumBytesTotal() const {
    return mTextReader->NumBytesTotal();
  }

  virtual SequenceReader< Alphabet >& operator>>( Sequence< Alphabet >& seq ) = 0;

  void Read( const size_t count, SequenceList< Alphabet >* out ) {
    Sequence< Alphabet > seq;

    for( size_t i = 0; i < count && !EndOfFile(); i++ ) {
      *this >> seq;
      out->push_back( std::move( seq ) );
    }
  }

  virtual ~SequenceReader() = default;

protected:
  std::unique_ptr< TextReader > mTextReader;
};
