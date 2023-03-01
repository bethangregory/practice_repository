/*
 * Interpolator.hpp
 *
 *  Created on: Oct 19, 2020
 *      Author: rodney
 */

#ifndef INTERPOLATOR_HPP_
#define INTERPOLATOR_HPP_

#include <vector>
#include <cmath>
using namespace std;

class Interpolator {
public:
	Interpolator(const vector<double> &x, const vector<double> &y);
	virtual ~Interpolator();

	// return linearly interpolated y value for given x value
	// if x outside boundaries returns either highest or lowest stored y value
	double linterp(double x);

	// return log-linearly interpolated y value for given x (use if y is on a logarithmic scale in input file)
	// if x outside boundaries returns either highest or lowest stored y value
	double loglinterp(double x);

private:
	int size;
	vector<double> x_data;
	vector<double> y_data;
	vector<double> log_y_data;
	vector<double> gradients;
	vector<double> log_y_gradients;
};

#endif /* INTERPOLATOR_HPP_ */
