#ifndef __NEWTON_RAPHSON_H
#define __NEWTON_RAPHSON_H

#include <cmath>

template<typename T, 
	double (T::*g)(double) const, 
	double (T::*g_prime)(double) const>

double newton_raphson(double y_target, double init, 
	double epsilon, const T& root_func){

	double y = (root_func.*g)(init);
	double x = init;

	while(fabs(y-y_target) > epsilon) {
		double d_x = (root_func.*g_prime)(x);
		x+=(y_target - y)/d_x;
		y = (root_func.*g)(x);
	}
	return x;
}

#endif