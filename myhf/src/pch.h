#pragma once

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <format>
#include <functional>
#include <iomanip>
#include <memory>
#include <numbers>
#include <span>
#include <string>
#include <vector>

#include <Eigen\eigen>

// gcem provides constexpr math functions 
// (Note: They may be slower than std math functions and should only be used at compile time)
#include "gcem.hpp"