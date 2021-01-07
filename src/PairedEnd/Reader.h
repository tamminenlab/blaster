#pragma once

#include "../FASTQ/Reader.h"
#include "../Sequence.h"

namespace PairedEnd {

template < typename Alphabet >
class Reader {
private:
  FASTQ::Reader< Alphabet > mFwdReader;
  FASTQ::Reader< Alphabet > mRevReader;

public:
  Reader( const std::string& pathToFwdFile, const std::string& pathToRevFile )
      : mFwdReader( pathToFwdFile ), mRevReader( pathToRevFile ) {}

  Reader( std::istream& fwd, std::istream& rev )
      : mFwdReader( fwd ), mRevReader( rev ) {}

  void Read( const int count, SequenceList< Alphabet >* fwds,
             SequenceList< Alphabet >* revs ) {
    Sequence< Alphabet > fwd, rev;

    for( int i = 0; i < count; i++ ) {
      if( !Read( &fwd, &rev ) )
        break;

      fwds->push_back( std::move( fwd ) );
      revs->push_back( std::move( rev ) );
    }
  }

  bool Read( Sequence< Alphabet >* fwd, Sequence< Alphabet >* rev ) {
    if( EndOfFile() )
      return false;

    mFwdReader >> *fwd;
    mRevReader >> *rev;

    return true;
  }

  bool EndOfFile() const {
    return mFwdReader.EndOfFile() || mRevReader.EndOfFile();
  }

  std::streampos NumBytesRead() const {
    return mFwdReader.NumBytesRead() + mRevReader.NumBytesRead();
  }

  std::streampos NumBytesTotal() const {
    return mFwdReader.NumBytesTotal() + mRevReader.NumBytesTotal();
  }
};

} // namespace PairedEnd
