#ifndef __INTERVAL_BISECTION_H
#define __INTERVAL_BISECTION_H

#include <cmath>

// Creating a function template
// Trying to find an x such that |g(x) - y| < epsilon,
// starting with the interval (m, n).
template<typename T>
double interval_bisection(double y_target,  // Target y value
                          double m,         // Left interval value
                          double n,         // Right interval value
                          double epsilon,   // Tolerance
                          T g) {            // Function object of type T, named g

  // Create the initial x mid-point value
  // Find the mapped y value of g(x)
  double x = 0.5 * (m + n);
  double y = g(x);

  // While the difference between y and the y_target
  // value is greater than epsilon, keep subdividing
  // the interval into successively smaller halves
  // and re-evaluate the new y.
  do {
    if (y < y_target) {
      m = x;
    }

    if (y > y_target) {
      n = x;
    }

    x = 0.5 * (m + n);
    y = g(x);
  } while (fabs(y-y_target) > epsilon);

  return x;
}

#endif