/*
 * Particle_O.hpp
 *
 *  Created on: Jun 24, 2020
 *      Author: rodney
 */

#ifndef PARTICLE_O_HPP_
#define PARTICLE_O_HPP_

#include "Particle.hpp"

class Particle_O: public Particle {
public:
	Particle_O();
	virtual ~Particle_O();
	static const double mass;
	static const string name;
	double get_mass() const;
	string get_name() const;
};

#endif /* PARTICLE_O_HPP_ */
