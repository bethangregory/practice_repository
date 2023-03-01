/*
 * Common_Functions.hpp
 *
 *  Created on: Aug 21, 2020
 *      Author: rodney
 */

#ifndef COMMON_FUNCTIONS_HPP_
#define COMMON_FUNCTIONS_HPP_

#include <vector>
#include <memory>
#include <iostream>
#include <sstream>
#include <fstream>
#include <algorithm>
#include <random>
#include <chrono>
#include <cmath>
using namespace std;

namespace constants {
	extern const double pi;       // pi [unitless]
	extern const double twopi;    // 2*pi [unitless]
	extern const double k_b;      // Boltzmann's Constant [erg/K]
	extern const double c;        // Speed of Light in Vacuum [cm/s]
	extern const double G;        // Gravitational Constant [cm^3/g/s^2]
	extern const double amu;      // Atomic Mass Unit [g]
	extern const double m_e;      // Electron Mass [g]
	extern const double q_e;      // Elementary Charge [C]
	extern const double jev;      // Joules/Electron Volt [unitless]
	extern const double ergev;    // ergs/Electron Volt [unitless]
}

namespace common {
	// functions for importing double type data from 2, 3, 4, 5, or 6-column csv files
	void import_csv(string filename, vector<double> &col1, vector<double> &col2);
	void import_csv(string filename, vector<double> &col1, vector<double> &col2, vector<double> &col3);
	void import_csv(string filename, vector<double> &col1, vector<double> &col2, vector<double> &col3, vector<double> &col4);
	void import_csv(string filename, vector<double> &col1, vector<double> &col2, vector<double> &col3, vector<double> &col4, vector<double> &col5);
	void import_csv(string filename, vector<double> &col1, vector<double> &col2, vector<double> &col3, vector<double> &col4, vector<double> &col5, vector<double> &col6);

	// returns interpolated value at x from parallel arrays (x_data, y_data)
	double interpolate(vector<double> &x_data, vector<double> &y_data, double x);

	// returns interpolated value at x from arrays x_data, y_data, when y_data is log-scaled
	double interpolate_logy(vector<double> &x_data, vector<double> &y_data, double x);

	// returns uniformly distributed random number from interval [0, 1)
	double get_rand();

	// returns uniformly distributed random integer between lower and upper (inclusive)
	int get_rand_int(int lower, int upper);
};

#endif /* COMMON_FUNCTIONS_HPP_ */
