/*
 * Distribution_Hot_H.hpp
 *
 *  Created on: Jul 31, 2020
 *      Author: rodney
 */

#ifndef DISTRIBUTION_HOT_H_HPP_
#define DISTRIBUTION_HOT_H_HPP_

#include "Distribution.hpp"

class Distribution_Hot_H: public Distribution {
public:
	Distribution_Hot_H(Planet my_p, double ref_h, double ref_T);
	virtual ~Distribution_Hot_H();
	void init(shared_ptr<Particle> p);
	double get_global_rate();

private:
	double m_H;                // [g] neutral H mass
	double m_Hplus;            // [g] H+ ion mass
	double m_HCOplus;          // [g] HCO+ ion mass
	double m_CO;               // [g] CO mass
  	double H_Hplus_rate_coeff; // [cm^3/s] rate coefficient for H+ + H -> H* + H+
	double HCOplus_DR_rate_coeff;  // [cm^3/s] rate coefficient for HCO+ + e -> H* + CO
        double any_mechanism_energy; // [eV] energy given to H atom when 'any_mechanism_prob' source chose (producing H at single altitude)
        int any_mechanism_alt_bin; // altitude bin in which H atoms are produced when source is 'any_mechanism_prob'
	double global_rate;    // [s^-1] set by chosen production method; calling function needs to divide this by 2 to get hemispherical rate
	string source;         // which source to use when initializing particles ('H_Hplus' or 'HCOplus_DR')
        vector<vector<double>> H_profile;
	vector<vector<double>> Hplus_profile;
	vector<vector<double>> HCOplus_profile;
	vector<vector<double>> electron_profile;
	vector<vector<double>> temp_profile;
	vector<vector<double>> H_Hplus_CDF;
	vector<vector<double>> HCOplus_DR_CDF;
        vector<vector<double>> any_mechanism_prob_CDF;

	// scans HCOplus_DR_CDF for new particle radius
	double get_new_radius_HCOplus_DR();

	// scans H_Hplus_CDF for new particle radius
	double get_new_radius_H_Hplus();

        // scans any_mechanism_prob_CDF for new particle radius
        double get_new_radius_any_mechanism_prob();

        // init particle using H_Hplus mechanism
	void init_H_Hplus_particle(shared_ptr<Particle> p);

	// init particle using HCOplus_DR mechanism
	void init_HCOplus_DR_particle(shared_ptr<Particle> p);

  	// init particle using any_mechanism_prob (to estimate escape probability)
        void init_any_mechanism_prob_particle(shared_ptr<Particle> p);

	// generate HCOplus_DR_CDF for given altitude range using imported density/temp profiles
	void make_HCOplus_DR_CDF(double lower_alt, double upper_alt);

	// generate H_Hplus_CDF for given altitude range using imported density/temp profiles
	void make_H_Hplus_CDF(double lower_alt, double upper_alt);

  	// generate any_mechanism_prob_CDF for given altitude range using imported density/temp profiles (to estimate escape probability)
        void make_any_mechanism_prob_CDF(double lower_alt, double upper_alt);
};

#endif /* DISTRIBUTION_HOT_H_HPP_ */
