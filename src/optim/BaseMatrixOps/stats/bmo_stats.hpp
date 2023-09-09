/*################################################################################
  ##
  ##   Copyright (C) 2011-2022 Keith O'Hara
  ##
  ##   This file is part of the StatsLib C++ library.
  ##
  ##   Licensed under the Apache License, Version 2.0 (the "License");
  ##   you may not use this file except in compliance with the License.
  ##   You may obtain a copy of the License at
  ##
  ##       http://www.apache.org/licenses/LICENSE-2.0
  ##
  ##   Unless required by applicable law or agreed to in writing, software
  ##   distributed under the License is distributed on an "AS IS" BASIS,
  ##   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
  ##   See the License for the specific language governing permissions and
  ##   limitations under the License.
  ##
  ################################################################################*/

#ifndef BMO_STATS_INCLUDES
#define BMO_STATS_INCLUDES

#ifndef BMO_RNG_ENGINE_TYPE
    #define BMO_RNG_ENGINE_TYPE std::mt19937_64;
#endif

namespace bmo
{
namespace stats
{

using rand_engine_t = BMO_RNG_ENGINE_TYPE;

template<typename T>
using return_t = typename std::conditional<std::is_integral<T>::value,double,T>::type;

template<typename ...T>
using common_t = typename std::common_type<T...>::type;

template<typename ...T>
using common_return_t = return_t<common_t<T...>>;

//

#include "runif.hpp"
#include "rnorm.hpp"

#include "rind.hpp"

}
}

#endif
