/*
 * Particle_H.cpp
 *
 *  Created on: Jun 19, 2020
 *      Author: rodney
 */

#include "Particle_H.hpp"

const double Particle_H::mass = 1.00794*constants::amu;
const string Particle_H::name = "H";

Particle_H::Particle_H() {

}

Particle_H::~Particle_H() {
}

double Particle_H::get_mass() const
{
	return mass;
}

string Particle_H::get_name() const
{
	return name;
}
