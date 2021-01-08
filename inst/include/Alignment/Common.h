#pragma once

#include "../Sequence.h"
#include "../Utils.h"

#include <limits>
inline int MinInt() { return std::numeric_limits< int >::min() / 2; }

enum class AlignmentDirection { Forward, Reverse };
