/*
 * Distribution_MB.cpp
 *
 *  Created on: Jul 21, 2020
 *      Author: rodney
 */

#include "Distribution_MB.hpp"

Distribution_MB::Distribution_MB(Planet my_p, double ref_h, double ref_T)
	: Distribution(my_p, ref_h, ref_T) {

}

Distribution_MB::~Distribution_MB() {

}

void Distribution_MB::init(shared_ptr<Particle> p)
{
	double phi = constants::twopi*(common::get_rand());
	double u = 2.0*common::get_rand() - 1;
	double v_avg = sqrt(constants::k_b*ref_temp/p->get_mass());

	double x = ref_radius*sqrt(1-(u*u))*cos(phi);
	double y = ref_radius*sqrt(1-(u*u))*sin(phi);
	double z = ref_radius*u;

	// Hemispherical Adjustment For Dayside Photochemical Process
	if (x > 0)
	{
		x = -x;
	}

	double v[] = {0.0, 0.0, 0.0};
	gen_mb(v_avg, v);

	p->init_particle(x, y, z, v[0], v[1], v[2]);
}

void Distribution_MB::init_vonly(shared_ptr<Particle> p, double v_avg)
{
	double v[] = {0.0, 0.0, 0.0};
	gen_mb(v_avg, v);

	p->init_particle_vonly(v[0], v[1], v[2]);
}

double Distribution_MB::get_global_rate()
{
	return 0.0;
}
