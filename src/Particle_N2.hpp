/*
 * Particle_N2.hpp
 *
 *  Created on: Jun 23, 2020
 *      Author: rodney
 */

#ifndef PARTICLE_N2_HPP_
#define PARTICLE_N2_HPP_

#include "Particle.hpp"

class Particle_N2: public Particle {
public:
	Particle_N2();
	virtual ~Particle_N2();
	static const double mass;
	static const string name;
	double get_mass() const;
	string get_name() const;
};

#endif /* PARTICLE_N2_HPP_ */
