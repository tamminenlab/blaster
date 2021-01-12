#pragma once

#include <string>
#include <iostream>
#include <sstream>
#include <map>
#include <iomanip>
#include <cmath>
#include <chrono>

enum class UnitType { COUNTS, BYTES };

static std::string ValueWithUnit( const float value, const UnitType unit ) {
  static const std::map< UnitType, std::map< size_t, std::string > > CONVERSION = {
    {
      UnitType::COUNTS,
      {
        { 1, "" },
        { 1000, "k" },
        { 1000 * 1000, "M" },
        { 1000 * 1000 * 1000, "G" }
      },
    },
    {
      UnitType::BYTES,
      {
        { 1, " Bytes" },
        { 1024, " kB" },
        { 1024 * 1024, " MB" },
        { 1024 * 1024 * 1024, " GB" }
      }
    }
  };

  auto list = CONVERSION.find( unit  )->second;
  auto it = list.begin();
  while( it != list.end() && value > it->first * 10 ) {
    it++;
  }

  std::stringstream ss;
  float dividedValue = value;
  if( it != list.begin() ) {
    it--;
    dividedValue = floor( value / it->first );
  }

  if( it->first == 1 ) {
    ss << std::setprecision(1) << std::setiosflags( std::ios::fixed );
  }

  ss << dividedValue;
  ss << it->second;
  return ss.str();
}

class ProgressOutput {
  using clock = std::chrono::steady_clock;

  using Stage = struct {
    std::string       label;
    UnitType          unit;
    size_t            value;
    size_t            max;
    clock::time_point lastPrint;
  };

public:
  ProgressOutput() : mActiveId( -1 ) {}

  ProgressOutput& Add( const int id, const std::string& label,
                       const UnitType unit = UnitType::COUNTS ) {
    mStages.insert( { id, Stage{ label, unit, 0, 100, clock::now() } } );
    return *this;
  }

  ProgressOutput& Set( const int id, const size_t value, const size_t max ) {
    auto& stage = mStages[ id ];
    stage.value = value;
    stage.max   = max;

    if( mActiveId == id ) {
      Print( stage );
    }
    return *this;
  }

  ProgressOutput& Activate( const int id ) {
    if( mActiveId != id )
      std::cerr << std::endl;

    mActiveId = id;
    Print( mStages[ id ] );
    return *this;
  }

private:
  void Print( Stage& stage ) {
    // Make sure we don't waste perf by outputting wasteful
    auto now    = clock::now();
    auto millis = std::chrono::duration_cast< std::chrono::milliseconds >(
                    now - stage.lastPrint )
                    .count();
    if( millis < 100 && stage.value != stage.max ) {
      return;
    }
    stage.lastPrint = now;

    std::ostream& os = std::cerr;

    size_t maxLabelLen = 0;
    for( auto& s : mStages ) {
      maxLabelLen = std::max( maxLabelLen, s.second.label.size() );
    }

    std::ios::fmtflags f( os.flags() );
    // Show one decimal point
    os << std::setiosflags( std::ios::fixed ) << std::setprecision( 1 );
    os << std::right << std::setw( maxLabelLen ) << stage.label << ": ";
    os << float( stage.value ) / stage.max * 100.0 << '%';
    os << " (" << ValueWithUnit( stage.value, stage.unit ) << ")";
    os << std::string( 20, ' ' ) << "\r" << std::flush;
    os.flags( f );
  }

  int mActiveId;
  std::map< int, Stage > mStages;
};
