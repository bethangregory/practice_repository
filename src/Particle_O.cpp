/*
 * Particle_O.cpp
 *
 *  Created on: Jun 24, 2020
 *      Author: rodney
 */

#include "Particle_O.hpp"

const double Particle_O::mass = 15.9994*constants::amu;
const string Particle_O::name = "O";

Particle_O::Particle_O() {

}

Particle_O::~Particle_O() {

}

double Particle_O::get_mass() const
{
	return mass;
}

string Particle_O::get_name() const
{
	return name;
}
