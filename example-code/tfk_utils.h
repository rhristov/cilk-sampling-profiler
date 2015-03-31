// Copyright Tim Kaler 2013

#ifndef TFK_UTILS_H_
#define TFK_UTILS_H_

#include <sys/time.h>
#include <sys/resource.h>
#include <stdio.h>
#include <sstream>
#include <string>

// Simple timer for benchmarks
double tfk_get_time() {
    struct timeval t;
    struct timezone tzp;
    gettimeofday(&t, &tzp);
    return t.tv_sec + t.tv_usec*1e-6;
}

// Simple [0,1] RNG
double tfkRand(double fMin, double fMax) {
    double f = static_cast<double>(rand()) / RAND_MAX;
    return fMin + f * (fMax - fMin);
}

template <typename T>
  std::string NumberToString(T Number) {
     std::ostringstream ss;
     ss << Number;
     return ss.str();
  }

template <typename T>
  T StringToNumber(std::string Text) {
     std::istringstream ss(Text);
     T result;
     return ss >> result ? result : 0;
  }


#endif  // TFK_UTILS_H_
