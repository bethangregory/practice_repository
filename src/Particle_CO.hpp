/*
 * Particle_CO.hpp
 *
 *  Created on: Jun 23, 2020
 *      Author: rodney
 */

#ifndef PARTICLE_CO_HPP_
#define PARTICLE_CO_HPP_

#include "Particle.hpp"

class Particle_CO: public Particle {
public:
	Particle_CO();
	virtual ~Particle_CO();
	static const double mass;
	static const string name;
	double get_mass() const;
	string get_name() const;
};

#endif /* PARTICLE_CO_HPP_ */
