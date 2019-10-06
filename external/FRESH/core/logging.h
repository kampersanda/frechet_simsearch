/*
 * Copyright 2018 Matteo Ceccarello
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *  http://www.apache.org/licenses/LICENSE-2.0
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 */
#pragma once

#ifndef _CURVEDIST_LOGGING_H
#define _CURVEDIST_LOGGING_H

#ifndef NDEBUG
#define PRINTD(x)                                       \
  do { std::stringstream s;                             \
    s << std::setprecision(2) << std::fixed << x;       \
    std::string msg = s.str();                          \
    std::cerr << msg << std::endl;                      \
  } while (false);

#define PRINTVAR(x) PRINTD("[" << #x << "] " << x)

#else
#define PRINTD(x)
#define PRINTVAR(x)
#endif

#endif // _CURVEDIST_LOGGING_H
