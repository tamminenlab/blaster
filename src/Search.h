#pragma once

#include <nsearch/Database/Search.h>

#include <string>

template < typename Alphabet >
extern bool DoSearch( const std::string&              queryPath,
                      const std::string&              databasePath,
                      const std::string&              outputPath,
                      const SearchParams< Alphabet >& searchParams );
