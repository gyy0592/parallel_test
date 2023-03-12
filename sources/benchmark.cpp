// benchmark.cpp
#include "benchmark.h"

void benchmark1(int iterations) {
  // declare a variable
  int a = 0;

  // perform some operations on the variable
  for (int i = 0; i < iterations; i++) {
    a++;
    a--;
    a=a+2;
    a=a*10086;
    a=a/10086;
    a=a-2;
  }
}