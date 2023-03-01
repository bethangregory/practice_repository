/*
 * Particle_CO.cpp
 *
 *  Created on: Jun 23, 2020
 *      Author: rodney
 */

#include "Particle_CO.hpp"

const double Particle_CO::mass = 28.0101*constants::amu;
const string Particle_CO::name = "CO";

Particle_CO::Particle_CO() {

}

Particle_CO::~Particle_CO() {

}

double Particle_CO::get_mass() const
{
	return mass;
}

string Particle_CO::get_name() const
{
	return name;
}
