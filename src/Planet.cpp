/*
 * Planet.cpp
 *
 *  Created on: May 29, 2020
 *      Author: rodney
 */

#include "Planet.hpp"

Planet::Planet()
{
	mass = 0.0;
	radius = 0.0;
	k_g = 0.0;
}

Planet::~Planet()
{

}

void Planet::init()
{
	mass = 4.8675e27;              // Venus mass by default [g]
	radius = 6.0518e8;             // Venus radius by default [cm]
	k_g = -(mass * constants::G);  // planetary gravitational constant [cm^3/s^2]
}

void Planet::init(double m, double r)
{
	mass = m;
	radius = r;
	k_g = -(mass * constants::G);
}

double Planet::get_mass()
{
	return mass;
}

double Planet::get_radius()
{
	return radius;
}

double Planet::get_k_g()
{
	return k_g;
}
