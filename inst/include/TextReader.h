#pragma once

#include <fstream>
#include <iostream>
#include <string>

#ifdef USE_ZLIB
#include <zlib.h>
#endif

class TextReader {
public:
  virtual size_t NumBytesRead() const  = 0;
  virtual size_t NumBytesTotal() const = 0;

  virtual bool EndOfFile() const = 0;

  virtual void operator>>( std::string& str ) = 0;

  virtual ~TextReader() = default;
};

class TextStreamReader : public TextReader {
public:
  TextStreamReader( std::istream& is );

  size_t NumBytesRead() const;
  size_t NumBytesTotal() const;

  bool EndOfFile() const;

  void operator>>( std::string& str );

private:
  std::istream&  mInput;
  std::streampos mTotalBytes;
};

/*
 * Reading huge fastq files needs to be fast
 */
class TextFileReader : public TextReader {
public:
  TextFileReader( const std::string& fileName,
                  const size_t       totalBufferSize = 32 * 1024 );
  ~TextFileReader();

  size_t NumBytesRead() const;
  size_t NumBytesTotal() const;

  bool EndOfFile() const;

  void operator>>( std::string& str );

private:
  void NextBuffer();

  int mFd;

#ifdef USE_ZLIB
  gzFile mGzFile;
#endif

  size_t mBufferPos, mBufferSize, mTotalBufferSize;
  char*  mBuffer;
  off_t  mTotalBytes;
};
