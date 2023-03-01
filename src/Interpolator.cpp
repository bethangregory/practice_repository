/*
 * Interpolator.cpp
 *
 *  Created on: Oct 19, 2020
 *      Author: rodney
 */

#include "Interpolator.hpp"

Interpolator::Interpolator(const vector<double> &x, const vector<double> &y) {

	// locally store x and y data
	size = x.size();
	x_data.resize(size);
	y_data.resize(size);
	log_y_data.resize(size);
	gradients.resize(size-1);
	log_y_gradients.resize(size-1);
	for (int i=0; i<size; i++)
	{
		x_data[i] = x[i];
		y_data[i] = y[i];
		log_y_data[i] = log10(y[i]);
	}

	// store all possible gradients for these x_data and y_data to speed interp calculation
	for (int i=0; i<size-1; i++)
	{
		gradients[i] = (y_data[i+1] - y_data[i]) / (x_data[i+1] - x_data[i]);
		log_y_gradients[i] = (log_y_data[i+1] - log_y_data[i]) / (x_data[i+1] - x_data[i]);
	}
}

Interpolator::~Interpolator() {

}


// return linearly interpolated y value for given x value
// if x outside boundaries returns either highest or lowest stored y value
double Interpolator::linterp(double x)
{
	double val = 0.0;

	if (x <= x_data[0])
	{
		val = y_data[0];
	}
	else if (x >= x_data.back())
	{
		val = y_data.back();
	}
	else
	{
		// locate index of first stored x_data value that is greater than our argument x value
		// this search method could probably be improved to be more efficient, but this is okay for now
		int i = lower_bound(x_data.begin(), x_data.end(), x) - x_data.begin();
		val = y_data[i-1] + gradients[i-1] * (x - x_data[i-1]);   // linear interpolation
	}

	return val;
}

// return log-linearly interpolated y value for given x (use if y is on a logarithmic scale in input file)
// if x outside boundaries returns either highest or lowest stored y value
double Interpolator::loglinterp(double x)
{
	double val = 0.0;

	if (x <= x_data[0])
	{
		val = pow(10.0, log_y_data[0]);
	}
	else if (x >= x_data.back())
	{
		val = pow(10.0, log_y_data.back());
	}
	else
	{
		// locate index of first stored x_data value that is greater than our argument x value
		// this search method could probably be improved to be more efficient, but this is okay for now
		int i = lower_bound(x_data.begin(), x_data.end(), x) - x_data.begin();
		val = log_y_data[i-1] + log_y_gradients[i-1] * (x - x_data[i-1]);   // linear interpolation
		val = pow(10.0, val);
	}

	return val;
}
