/*
 * Distribution.cpp
 *
 *  Created on: Jul 18, 2020
 *      Author: rodney
 */

#include "Distribution.hpp"

Distribution::Distribution(Planet my_p, double ref_h, double ref_T) {
	my_planet = my_p;
	ref_height = ref_h;
	ref_radius = my_planet.get_radius() + ref_height;
	ref_temp = ref_T;
}

Distribution::~Distribution() {

}

void Distribution::gen_mb(double vavg, double v_in[])
{
	double randnum1 = common::get_rand();
	double randnum2 = common::get_rand();
	double randnum3 = common::get_rand();
	double randnum4 = common::get_rand();

	randnum1 = vavg*sqrt(-2.0*log(1.0-randnum1));
	randnum2 = constants::twopi*randnum2;
	randnum3 = vavg*sqrt(-2.0*log(1.0-randnum3));
	randnum4 = constants::twopi*randnum4;

	v_in[0] = randnum1*cos(randnum2);
	v_in[1] = randnum1*sin(randnum2);
	v_in[2] = randnum3*cos(randnum4);
}
