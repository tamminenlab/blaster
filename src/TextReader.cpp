#include "TextReader.h"
#include "Utils.h"

#include <fcntl.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include <algorithm>
#include <cstring>

/*
 * TextStreamReader
 */
TextStreamReader::TextStreamReader( std::istream& is ) : mInput( is ) {
  mInput.seekg( 0, mInput.end );
  mTotalBytes = mInput.tellg();
  mInput.seekg( 0, mInput.beg );
}

size_t TextStreamReader::NumBytesRead() const {
  std::streampos pos = mInput.tellg();
  return pos < 0 ? mTotalBytes : pos;
}

size_t TextStreamReader::NumBytesTotal() const {
  return mTotalBytes;
}

bool TextStreamReader::EndOfFile() const {
  return !mInput || mInput.peek() == EOF;
}

inline bool IsBlank( const std::string& str ) {
  return str.empty() || std::all_of( str.begin(), str.end(), isspace );
}

void TextStreamReader::operator>>( std::string& str ) {
  do {
    getline( mInput, str );
  } while( !EndOfFile() && IsBlank( str ) );
}

/*
 * TextFileReader
 */
void TextFileReader::NextBuffer() {
#ifdef USE_ZLIB
  if( mGzFile ) {
    mBufferSize = gzread( mGzFile, mBuffer, mTotalBufferSize );
  } else
#endif
  {
    mBufferSize = read( mFd, mBuffer, mTotalBufferSize );
  }

  mBufferPos = 0;
}

TextFileReader::TextFileReader( const std::string& fileName,
                                const size_t       totalBufferSize )
    : mBufferPos( -1 ), mBufferSize( 0 ), mTotalBufferSize( totalBufferSize ),
      mBuffer( NULL ) {
  mFd = open( fileName.c_str(), O_RDONLY ); // orly?

  if( mFd != -1 ) {
#ifdef USE_ZLIB
    mGzFile = NULL;

    // Check for GZ magic number
    uint8_t magic[ 2 ] = { 0, 0 };
    read( mFd, magic, 2 );
    lseek( mFd, 0, SEEK_SET );
    if( magic[ 0 ] == 0x1F && magic[ 1 ] == 0x8B ) {
      mGzFile = gzdopen( mFd, "rb" );
    }
#endif
    mBuffer = new char[ totalBufferSize ];

    mTotalBytes = lseek( mFd, 0, SEEK_END );
    lseek( mFd, 0, SEEK_SET );

    NextBuffer();
  }
}

TextFileReader::~TextFileReader() {
  if( mBuffer ) {
    delete[] mBuffer;
  }

  if( mFd != -1 ) {
#ifdef USE_ZLIB
    if( mGzFile ) {
      gzclose( mGzFile );
    } else
#endif
    {
      close( mFd );
    }
  }
}

void TextFileReader::operator>>( std::string& str ) {
  str.clear();
ReadLine:
  while( !EndOfFile() ) {
    char* pos =
      ( char* ) memchr( mBuffer + mBufferPos, '\n', mBufferSize - mBufferPos );

    if( pos == NULL ) {
      str += std::string( mBuffer + mBufferPos, mBufferSize - mBufferPos );
      NextBuffer();
    } else {
      size_t numBytes = pos - ( mBuffer + mBufferPos );
      str += std::string( mBuffer + mBufferPos, numBytes );

      mBufferPos += numBytes + 1; // skip '\n'
      if( mBufferPos >= mBufferSize )
        NextBuffer();

      break;
    }
  }

  if( IsBlank( str ) && !EndOfFile() )
    goto ReadLine;
}

bool TextFileReader::EndOfFile() const {
  return mFd == -1 || mBufferSize <= 0;
}

size_t TextFileReader::NumBytesRead() const {
  if( EndOfFile() ) {
    return mTotalBytes;
  } else {
    return lseek( mFd, 0, SEEK_CUR );
  }
}

size_t TextFileReader::NumBytesTotal() const {
  return mTotalBytes;
}
