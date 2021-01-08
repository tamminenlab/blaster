#pragma once

#include <string>
#include <vector>
#include <map>

// Sequences
#include "FASTA/Writer.h"
#include "FASTQ/Writer.h"
#include "FASTA/Reader.h"
#include "FASTQ/Reader.h"

// Hits
#include "Alnout/Writer.h"
#include "CSV/Writer.h"

enum class FileFormat {
  FASTA,
  FASTQ,
  ALNOUT,
  CSV,
};

using StringList = std::vector< std::string > ;
static const std::map< FileFormat, StringList > FileFormatEndings = {
  { FileFormat::FASTA, { "fa", "fna", "fsa", "fasta" } },
  { FileFormat::FASTQ, { "fq", "fastq" } },
  { FileFormat::ALNOUT, { "aln", "alnout" } },
  { FileFormat::CSV, { "csv" } },
};

static FileFormat InferFileFormat( const std::string& filepath, const FileFormat defaultFormat ) {
  auto pos = filepath.find_last_of( "." );
  if( pos == std::string::npos )
    return defaultFormat;

  auto ext = filepath.substr( pos + 1 );
  for( auto &ff : FileFormatEndings ) {
    for( auto &ending : ff.second ) {
      if( ext == ending )
        return ff.first;
    }
  }
  return defaultFormat;
}

template < typename A >
static std::unique_ptr< SequenceWriter < A > > DetectFileFormatAndOpenWriter( const std::string &path, const FileFormat defaultFormat ) {
  switch( InferFileFormat( path, defaultFormat ) ) {
    case FileFormat::FASTQ:
      return std::unique_ptr< SequenceWriter< A > >(
        new FASTQ::Writer< A >( path ) );

    default:
      return std::unique_ptr< SequenceWriter< A > >(
        new FASTA::Writer< A >( path ) );
  }
}

template < typename A >
static std::unique_ptr< SequenceReader < A > > DetectFileFormatAndOpenReader( const std::string &path, const FileFormat defaultFormat ) {
  switch( InferFileFormat( path, defaultFormat ) ) {
    case FileFormat::FASTQ:
      return std::unique_ptr< SequenceReader< A > >(
        new FASTQ::Reader< A >( path ) );

    default:
      return std::unique_ptr< SequenceReader< A > >(
        new FASTA::Reader< A >( path ) );
  }
}

template < typename A >
static std::unique_ptr< HitWriter < A > > DetectFileFormatAndOpenHitWriter( const std::string &path, const FileFormat defaultFormat ) {
  switch( InferFileFormat( path, defaultFormat ) ) {
    case FileFormat::CSV:
      return std::unique_ptr< HitWriter< A > >(
        new CSV::Writer< A >( path ) );

    default:
      return std::unique_ptr< HitWriter< A > >(
        new Alnout::Writer< A >( path ) );
  }
}
