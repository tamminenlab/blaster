#pragma once

#include <string>
#include <cassert>
#include <ctype.h>
#include <numeric>
#include <Rcpp.h>

static void UpcaseString( std::string& str ) {
  for( auto& ch : str )
    if( ch >= 97 && ch <= 122 ) // upcase
      ch &= ~0x20;
}
