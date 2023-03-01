/*
 * Planet.hpp
 *
 *  Created on: May 29, 2020
 *      Author: rodney
 */

#ifndef PLANET_HPP_
#define PLANET_HPP_

#include "Common_Functions.hpp"

class Planet {
public:
	Planet();
	virtual ~Planet();
	void init(); // initialize with defaults (Venus mass/radius)
	void init(double m, double r);
	double get_mass();
	double get_radius();
	double get_k_g();

private:
	double mass;    // [g] mass of planet
	double radius;  // [cm] radius of planet
	double k_g;     // [cm^3/s^2] planet's gravitational constant (-G*mass)
};

#endif /* PLANET_HPP_ */
