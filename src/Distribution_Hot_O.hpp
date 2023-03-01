/*
 * Distribution_Hot_O.hpp
 *
 *  Created on: Jul 22, 2020
 *      Author: rodney
 */

#ifndef DISTRIBUTION_HOT_O_HPP_
#define DISTRIBUTION_HOT_O_HPP_

#include "Distribution.hpp"

class Distribution_Hot_O: public Distribution {
public:
	Distribution_Hot_O(Planet my_p, double ref_h, double ref_T);
	virtual ~Distribution_Hot_O();
	void init(shared_ptr<Particle> p);
	double get_global_rate();

private:
	// these 4 are from Justin's original Fortran code; may still need at some point
	double DR160;  // [cm^3/s] dissociative recombination rate at 200km
	double H_DR;   // [cm] dissociative recombination scale height
	double T_ion;  // [K] ion temperature
	double T_e;    // [K] electron temperature

	// everything below here is needed for the new process
	double m_O2plus;  // [g] O2+ ion mass
	double m_O;    // [g] neutral O mass

	double O2plus_DR_rate_coeff;  // [cm^3/s] rate coefficient for O2+ + e -> O* + O*
	double global_rate;    // [s^-1] set by chosen production method; calling function needs to divide this by 2 to get hemispherical rate
	string source;         // which source to use when initializing particles (currently 'O2plus_DR' is only option)
	vector<vector<double>> O2plus_profile;  // stores the imported O2plus density profile
	vector<vector<double>> electron_profile; // stores the imported electron density profile
	vector<vector<double>> temp_profile;  // stores the imported temperature profiles (4-column csv expected: altitude(cm), neutral_temp(K), ion_temp(K), electron_temp(K))
	vector<vector<double>> O2plus_DR_CDF; // stores an inverse CDF made using the O2plus_DR production rate profile

	// scans O2plus_DR_CDF for new particle radius
	double get_new_radius_O2plus_DR();

	// init particle using O2plus_DR mechanism
	void init_O2plus_DR_particle(shared_ptr<Particle> p);

	// init particle using Justin's original method
	void init_old_way(shared_ptr<Particle> p);

	// generate O2plus_DR_CDF for given altitude range using imported density/temp profiles
	void make_O2plus_DR_CDF(double lower_alt, double upper_alt);

	// calculate rotational energy of O2plus ion
	double E_rot(double B, double T);
};

#endif /* DISTRIBUTION_HOT_O_HPP_ */
