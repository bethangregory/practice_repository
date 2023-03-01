/*
 * Distribution_Hot_H.cpp
 *
 *  Created on: Jul 31, 2020
 *      Author: rodney
 */

#include "Distribution_Hot_H.hpp"

Distribution_Hot_H::Distribution_Hot_H(Planet my_p, double ref_h, double ref_T)
	: Distribution(my_p, ref_h, ref_T) {

	m_H = 1.00794*constants::amu;
	m_Hplus = 1.00728*constants::amu;
	m_HCOplus = 29.0175*constants::amu;
	m_CO = 28.0101*constants::amu;
	source = "";
	global_rate = 0.0;

	temp_profile.resize(4);
	H_profile.resize(2);
	Hplus_profile.resize(2);
	HCOplus_profile.resize(2);
	electron_profile.resize(2);
	H_Hplus_CDF.resize(2);
	HCOplus_DR_CDF.resize(2);
	any_mechanism_prob_CDF.resize(2);

	// import parameters from configuration file
	// file must be named 'Hot_H.cfg' and be in same directory as executable
	double profile_bottom = 0.0;
	double profile_top = 0.0;
	string temp_prof_filename = "";
	string H_prof_filename = "";
	string Hplus_prof_filename = "";
	string HCOplus_prof_filename = "";
	string electron_prof_filename = "";

	ifstream infile;
	infile.open("Hot_H.cfg");
	if (!infile.good())
	{
		cout << "Hot H configuration file not found!\n";
		exit(1);
	}
	string line, param, val;
	vector<string> parameters;
	vector<string> values;
	int num_params = 0;

	while (getline(infile, line))
	{
		if (line[0] == '#' || line.empty() || std::all_of(line.begin(), line.end(), ::isspace))
		{
			continue;
		}
		else
		{
			stringstream str(line);
			str >> param >> val;
			parameters.push_back(param);
			values.push_back(val);
			num_params++;
			param = "";
			val = "";
		}
	}
	infile.close();

	for (int i=0; i<num_params; i++)
	{
		if (parameters[i] == "H_Hplus_rate_coeff")
		{
			H_Hplus_rate_coeff = stod(values[i]);
		}
		else if (parameters[i] == "HCOplus_DR_rate_coeff")
		{
			HCOplus_DR_rate_coeff = stod(values[i]);
		}
		else if (parameters[i] == "source")
		{
			source = values[i];
		}
		else if (parameters[i] == "profile_bottom")
		{
			profile_bottom = stod(values[i]);
		}
		else if (parameters[i] == "profile_top")
		{
			profile_top = stod(values[i]);
		}
		else if (parameters[i] == "temp_prof_filename")
		{
			temp_prof_filename = values[i];
		}
		else if (parameters[i] == "H_prof_filename")
		{
			H_prof_filename = values[i];
		}
		else if (parameters[i] == "Hplus_prof_filename")
		{
			Hplus_prof_filename = values[i];
		}
		else if (parameters[i] == "HCOplus_prof_filename")
		{
			HCOplus_prof_filename = values[i];
		}
		else if (parameters[i] == "electron_prof_filename")
		{
			electron_prof_filename = values[i];
		}
		else if (parameters[i] == "any_mechanism_energy")
		{
		  any_mechanism_energy = stod(values[i]);
		}
		else if (parameters[i] == "any_mechanism_alt_bin")
		{
		  any_mechanism_alt_bin = stoi(values[i]);
		}
		    
	}
	
	common::import_csv(temp_prof_filename, temp_profile[0], temp_profile[1], temp_profile[2], temp_profile[3]);
	common::import_csv(H_prof_filename, H_profile[0], H_profile[1]);
	common::import_csv(Hplus_prof_filename, Hplus_profile[0], Hplus_profile[1]);
	common::import_csv(HCOplus_prof_filename, HCOplus_profile[0], HCOplus_profile[1]);
	common::import_csv(electron_prof_filename, electron_profile[0], electron_profile[1]);

	if (source == "H_Hplus")
	{
		make_H_Hplus_CDF(profile_bottom, profile_top);
	}
	else if (source == "HCOplus_DR")
	{
		make_HCOplus_DR_CDF(profile_bottom, profile_top);
	}
	else if (source == "any_mechanism_prob")
	{
		make_any_mechanism_prob_CDF(profile_bottom, profile_top);
	}
}

Distribution_Hot_H::~Distribution_Hot_H() {

}

void Distribution_Hot_H::init(shared_ptr<Particle> p)
{
	if (source == "H_Hplus")
	{
		init_H_Hplus_particle(p);
	}
	else if (source == "HCOplus_DR")
	{
		init_HCOplus_DR_particle(p);
	}
	else if (source == "any_mechanism_prob")
	{
	         init_any_mechanism_prob_particle(p);
	}
}

// init particle using H_Hplus mechanism
void Distribution_Hot_H::init_H_Hplus_particle(shared_ptr<Particle> p)
{
	// altitude distribution for hot H
	double r = get_new_radius_H_Hplus();
	double alt = r - my_planet.get_radius();
	double temp_ion = common::interpolate_logy(temp_profile[0], temp_profile[2], alt);
	double temp_neut = common::interpolate_logy(temp_profile[0], temp_profile[1], alt);

	double phi = constants::twopi*(common::get_rand());
	double u = 2.0*common::get_rand() - 1.0;
	double x = r*sqrt(1-(u*u))*cos(phi);
	double y = r*sqrt(1-(u*u))*sin(phi);
	double z = r*u;

	// Hemispherical Adjustment For Dayside Photochemical Process
	if (x < 0)
	{
		x = -x;
	}

	double ion_vavg = sqrt(constants::k_b*temp_ion/m_Hplus);  // average thermal H+ velocity
	double neut_vavg = sqrt(constants::k_b*temp_neut/m_H);    // average thermal H velocity
	double v_ion[] = {0.0, 0.0, 0.0};
	double v_neut[] = {0.0, 0.0, 0.0};
	gen_mb(ion_vavg, v_ion);
	gen_mb(neut_vavg, v_neut);

	double e = 0.0;
	Matrix<double, 3, 1> p1_v = {v_neut[0], v_neut[1], v_neut[2]};
	Matrix<double, 3, 1> p2_v = {v_ion[0], v_ion[1], v_ion[2]};
	Matrix<double, 3, 1> vcm;
	vcm = (m_H*p1_v.array() + m_Hplus*p2_v.array()) / (m_H + m_Hplus);
	Matrix<double, 3, 1> p1_vcm = p1_v.array() - vcm.array();    // particle 1 c-o-m velocity
	Matrix<double, 3, 1> p2_vcm = p2_v.array() - vcm.array();    // particle 2 c-o-m velocity
	double p1_vcm_tot = sqrt(p1_vcm[0]*p1_vcm[0] + p1_vcm[1]*p1_vcm[1] + p1_vcm[2]*p1_vcm[2]);  // particle 1 c-o-m scalar velocity
	double p2_vcm_tot = sqrt(p2_vcm[0]*p2_vcm[0] + p2_vcm[1]*p2_vcm[1] + p2_vcm[2]*p2_vcm[2]);  // particle 2 c-o-m scalar velocity
	e = 0.5*m_H*p1_vcm_tot*p1_vcm_tot + 0.5*m_Hplus*p2_vcm_tot*p2_vcm_tot;

	double v = sqrt(2.0*e / (m_H + (m_H*m_H/m_Hplus)));

	// spherically isotropic velocity vector
	phi = constants::twopi*common::get_rand();
	u = 2.0*common::get_rand() - 1.0;
	double vx = v*sqrt(1-u*u)*cos(phi);
	double vy = v*sqrt(1-u*u)*sin(phi);
	double vz = v*u;

	//p->init_particle(x, y, z, v_ion[0], v_ion[1], v_ion[2]);
	p->init_particle(x, y, z, vx, vy, vz);
}

// init particle using HCOplus_DR mechanism
void Distribution_Hot_H::init_HCOplus_DR_particle(shared_ptr<Particle> p)
{
	// altitude distribution for hot H
	double r = get_new_radius_HCOplus_DR();
	double alt = r - my_planet.get_radius();
	double temp_ion = common::interpolate_logy(temp_profile[0], temp_profile[2], alt);
	//double temp_neut = common::interpolate_logy(temp_profile[0], temp_profile[1], alt);
	double temp_e = common::interpolate_logy(temp_profile[0], temp_profile[3], alt);

	double phi = constants::twopi*(common::get_rand());
	double u = 2.0*common::get_rand() - 1.0;
	double x = r*sqrt(1-(u*u))*cos(phi);
	double y = r*sqrt(1-(u*u))*sin(phi);
	double z = r*u;

	// Hemispherical Adjustment For Dayside Photochemical Process
	if (x < 0)
	{
		x = -x;
	}

	// Select HCO+ DR Electronic Channel
	double Ei = 0.0;       // excess energy from DR (eV)
	double vib_lvl = 0.0;  // vibrational level of CO after DR
	double we = 0.0;       // vibrational frequency (cm-1)
	double wexe = 0.0;     // first correction term to vibrational frequency (cm-1)

	double randnum = common::get_rand();
	if (randnum < 0.23)
	{
		Ei = 1.3;
		we = 1743.41;
		wexe = 14.36;
		double randnum2 = common::get_rand();
		if (randnum2 < 0.45)
		{
			vib_lvl = 0.0;
		}
		else if (randnum2 < 0.45+0.21)
		{
			vib_lvl = 1.0;
		}
		else if (randnum2 < 0.45+0.21+0.13)
		{
			vib_lvl = 2.0;
		}
		else if (randnum2 < 0.45+0.21+0.13+0.1)
		{
			vib_lvl = 3.0;
		}
		else if (randnum2 < 0.45+0.21+0.13+0.1+0.1)
		{
			vib_lvl = 4.0;
		}
		else
		{
			vib_lvl = 5.0;
		}
	}
	/* assuming negligible chance of going to CO(a') state, but leaving code in for possible future use
	else if (randnum < (0.23+0.385))
	{
		Ei = 0.44;
		we = 1228.6;
		wexe = 10.468;
		double randnum2 = common::get_rand();
		if (randnum2 < 0.5)
		{
			vib_lvl = 0.0;
		}
		else if (randnum2 < 0.5+0.25)
		{
			vib_lvl = 1.0;
		}
		else if (randnum2 < 0.5+0.25+0.125)
		{
			vib_lvl = 2.0;
		}
		else
		{
			vib_lvl = 3.0;
		}
	}
	*/
	else
	{
		Ei = 7.31;
		we = 2169.81358;
		wexe = 13.28831;
		/*
		double randnum2 = common::get_rand();
		if (randnum2 < 0.5)
		{
			vib_lvl = 0.0;
		}
		else if (randnum2 < 0.5+0.25)
		{
			vib_lvl = 1.0;
		}
		else if (randnum2 < 0.5+0.25+0.125)
		{
			vib_lvl = 2.0;
		}
		else if (randnum2 < 0.5+0.25+0.125+0.0625)
		{
			vib_lvl = 3.0;
		}
		else if (randnum2 < 0.5+0.25+0.125+0.0625+0.03125)
		{
			vib_lvl = 4.0;
		}
		else if (randnum2 < 0.5+0.25+0.125+0.0625+0.03125+0.015625)
		{
			vib_lvl = 5.0;
		}
		else if (randnum2 < 0.5+0.25+0.125+0.0625+0.03125+0.015625+0.0078125)
		{
			vib_lvl = 6.0;
		}
		else
		{
			vib_lvl = 7.0;
		}
		*/
		vib_lvl = 0.0;
	}

	// subtract vibrational energy from excess if vib_lvl above 0 (convert inverse cm to eV before subtracting)
	if (vib_lvl > 0.0)
	{
		Ei = Ei + ((we*0.5 - wexe*0.25) * 1.239841984e-4);
		Ei = Ei - ((we*(vib_lvl + 0.5) - (wexe*(vib_lvl + 0.5)*(vib_lvl + 0.5))) * 1.239841984e-4);
		if (Ei < 0.0)  // this only happens when considering CO(a') state but leaving in just in case...
		{
			Ei = 0.0;
		}
	}

	// convert energy to ergs
	Ei = Ei*constants::ergev;

	// Translational Energy of H Resulting from Dissociative Recombination of HCO+
	double v = sqrt(2.0*Ei / (m_H + (m_H*m_H/m_CO)));

	// spherically isotropic velocity vector
	phi = constants::twopi*common::get_rand();
	u = 2.0*common::get_rand() - 1.0;
	double vx = v*sqrt(1-u*u)*cos(phi);
	double vy = v*sqrt(1-u*u)*sin(phi);
	double vz = v*u;

	// Add initial HCO+ and e translational momentum
	double vavg = sqrt(constants::k_b*temp_ion/(m_HCOplus));  // average thermal ion velocity
	double v_ion[] = {0.0, 0.0, 0.0};
	gen_mb(vavg, v_ion);                               // sample from Maxwell-Boltzmann distribution
	vavg = sqrt(constants::k_b*temp_e/constants::m_e);    //average thermal e velocity
	double v_e[] = {0.0, 0.0, 0.0};
	gen_mb(vavg, v_e);                                 // sample from Maxwell-Boltzmann distribution

	// add it all up
	vx = vx + (m_HCOplus*v_ion[0] + constants::m_e*v_e[0]) / (m_HCOplus+constants::m_e);
	vy = vy + (m_HCOplus*v_ion[1] + constants::m_e*v_e[1]) / (m_HCOplus+constants::m_e);
	vz = vz + (m_HCOplus*v_ion[2] + constants::m_e*v_e[2]) / (m_HCOplus+constants::m_e);

	p->init_particle(x, y, z, vx, vy, vz);
}

// init particle using any_mechanism_prob (to compute  escape rate of H atoms born at a certain altitude)
void Distribution_Hot_H::init_any_mechanism_prob_particle(shared_ptr<Particle> p)
{
	// altitude distribution for hot H
        double r = get_new_radius_any_mechanism_prob(); // radius from centre of planet (cm)	
	double alt = r - my_planet.get_radius(); // altitude from surface of planet (cm)
	double temp_ion = common::interpolate_logy(temp_profile[0], temp_profile[2], alt);
	//double temp_neut = common::interpolate_logy(temp_profile[0], temp_profile[1], alt);
	double temp_e = common::interpolate_logy(temp_profile[0], temp_profile[3], alt);

	double phi = constants::twopi*(common::get_rand());
	double u = 2.0*common::get_rand() - 1.0;
	double x = r*sqrt(1-(u*u))*cos(phi);
	double y = r*sqrt(1-(u*u))*sin(phi);
	double z = r*u;

	// Hemispherical Adjustment For Dayside Photochemical Process
	if (x < 0)
	{
		x = -x;
	}

	// Select Electronic Channel
	double Ei = any_mechanism_energy;       // (eV). Assume representative excess energy of x eV for all mechanisms, with none to electronic, vibrational or rotational levels

	// convert energy to ergs
	Ei = Ei*constants::ergev;

	// Translational Energy of H Resulting from excess energy. Assume all energy goes to H, with none to a secondary particle. Note that this results in v between 1 (for large mass secondary particle) and 1.41 (for secondary particle = H) x larger than if there was a secondary particle.
	double v = sqrt(2.0*Ei /m_H); // units of cm/s

	// spherically isotropic velocity vector
	phi = constants::twopi*common::get_rand();
	u = 2.0*common::get_rand() - 1.0;
	double vx = v*sqrt(1-u*u)*cos(phi);
	double vy = v*sqrt(1-u*u)*sin(phi);
	double vz = v*u;

       	//std::cout << "alt " << alt/1e5 << "\t" << v << "\n"; // This line allows check that particles are being produced at the altitude (km) and with the velocity (cm/s) expected

	// Don't add initial translational momentum of reactants -  this is acceptable because the velocities are isotropic, so there will be a net 0 effect on v
	p->init_particle(x, y, z, vx, vy, vz);
}

// scans HCOplus_DR_CDF for new particle radius
double Distribution_Hot_H::get_new_radius_HCOplus_DR()
{
        double u = common::get_rand(); // random no between 0 and 1
	int k = 0;
	while (HCOplus_DR_CDF[0][k] < u)
	{
		k++;
	}
	return HCOplus_DR_CDF[1][k] + my_planet.get_radius();
}

// scans H_Hplus_CDF for new particle radius
double Distribution_Hot_H::get_new_radius_H_Hplus()
{
	double u = common::get_rand();
	int k = 0;
	while (H_Hplus_CDF[0][k] < u)
	{
		k++;
	}
	return H_Hplus_CDF[1][k] + my_planet.get_radius();
}

// scans any_mechanism_prob_CDF for new particle radius
double Distribution_Hot_H::get_new_radius_any_mechanism_prob()
{
	double u = common::get_rand();
	int k = 0;
	while (any_mechanism_prob_CDF[0][k] < u)
	{
		k++;
	}
	return any_mechanism_prob_CDF[1][k] + my_planet.get_radius();
}

// returns global production rate (needs to be set by chosen production method)
// calling function must divide this by 2 if only hemispherical rate is needed
double Distribution_Hot_H::get_global_rate()
{
	return global_rate;
}

// generate H_Hplus_CDF for given altitude range using imported density/temp profiles
void Distribution_Hot_H::make_H_Hplus_CDF(double lower_alt, double upper_alt)
{
	//ofstream outfile;
	//outfile.open("/home/rodney/Documents/coronaTest/rodney_hplh.dat");

	double bin_size = 10000.0; // [cm]
	int num_alt_bins = (int)((upper_alt - lower_alt) / bin_size);
	vector<double> H_Hplus_rate;
	H_Hplus_rate.resize(num_alt_bins);
	H_Hplus_CDF[0].resize(num_alt_bins);
	H_Hplus_CDF[1].resize(num_alt_bins);

	double bottom_Hplus = Hplus_profile[0][0];
	double top_Hplus = Hplus_profile[0].back();
	double bottom_H = H_profile[0][0];
	double top_H = H_profile[0].back();

	double local_g = (constants::G * my_planet.get_mass()) / (pow(my_planet.get_radius()+bottom_Hplus, 2.0));
	double Hplus_bottom_scaleheight = constants::k_b*common::interpolate_logy(temp_profile[0], temp_profile[2], bottom_Hplus)/(m_Hplus*local_g);
	local_g = (constants::G * my_planet.get_mass()) / (pow(my_planet.get_radius()+top_Hplus, 2.0));
	double Hplus_top_scaleheight = constants::k_b*common::interpolate_logy(temp_profile[0], temp_profile[2], top_Hplus)/(m_Hplus*local_g);
	local_g = (constants::G * my_planet.get_mass()) / (pow(my_planet.get_radius()+bottom_H, 2.0));
	double H_bottom_scaleheight = constants::k_b*common::interpolate_logy(temp_profile[0], temp_profile[1], bottom_H)/(m_H*local_g);
	local_g = (constants::G * my_planet.get_mass()) / (pow(my_planet.get_radius()+top_H, 2.0));
	double H_top_scaleheight = constants::k_b*common::interpolate_logy(temp_profile[0], temp_profile[1], top_H)/(m_H*local_g);

	double rate_sum = 0.0;
	double rate_sum_times_r_sqrd = 0.0;
	for (int i=0; i<num_alt_bins; i++)
	{
		H_Hplus_CDF[1][i] = lower_alt + bin_size*i;

		// get new H density by either interpolation or extrapolation
		double H_dens = 0.0;
		if (H_Hplus_CDF[1][i] < bottom_H)
		{
			H_dens = H_profile[1][0]*exp((bottom_H - H_Hplus_CDF[1][i])/H_bottom_scaleheight);
		}
		else if (H_Hplus_CDF[1][i] > top_H)
		{
			H_dens = H_profile[1].back()*exp((top_H - H_Hplus_CDF[1][i])/H_top_scaleheight);
		}
		else
		{
			H_dens = common::interpolate_logy(H_profile[0], H_profile[1], H_Hplus_CDF[1][i]);
		}

		// get new Hplus density by either interpolation or extrapolation
		double Hplus_dens = 0.0;
		if (H_Hplus_CDF[1][i] < bottom_Hplus)
		{
			Hplus_dens = Hplus_profile[1][0]*exp((bottom_Hplus - H_Hplus_CDF[1][i])/Hplus_bottom_scaleheight);
		}
		else if (H_Hplus_CDF[1][i] > top_Hplus)
		{
			Hplus_dens = Hplus_profile[1].back()*exp((top_Hplus - H_Hplus_CDF[1][i])/Hplus_top_scaleheight);
		}
		else
		{
			Hplus_dens = common::interpolate_logy(Hplus_profile[0], Hplus_profile[1], H_Hplus_CDF[1][i]);
		}

		double Ti = common::interpolate_logy(temp_profile[0], temp_profile[2], H_Hplus_CDF[1][i]);
		H_Hplus_rate[i] = H_Hplus_rate_coeff * sqrt(Ti) * H_dens * Hplus_dens;
		//outfile << H_Hplus_CDF[1][i]*1e-5 << "\t" << H_Hplus_rate[i] << "\n";
		rate_sum = rate_sum + H_Hplus_rate[i];
		rate_sum_times_r_sqrd = rate_sum_times_r_sqrd + (H_Hplus_rate[i] * 4.0 * constants::pi * pow(my_planet.get_radius()+(lower_alt+bin_size*i), 2.0));
	}
	//outfile.close();
	for (int i=0; i<num_alt_bins; i++)
	{
		if (i == 0)
		{
			H_Hplus_CDF[0][i] = H_Hplus_rate[i] / rate_sum;
		}
		else
		{
			H_Hplus_CDF[0][i] = (H_Hplus_rate[i] / rate_sum) + H_Hplus_CDF[0][i-1];
		}
	}
	//double global_rate_H_Hplus = rate_sum*bin_size * 4.0 * constants::pi * pow(my_planet.get_radius()+upper_alt, 2.0);
	double global_rate_H_Hplus = rate_sum_times_r_sqrd*bin_size;
	cout << "Global hot H production rate from H+ + H:\n" << global_rate_H_Hplus << " per second\n";
	global_rate = global_rate_H_Hplus;
}

// generate HCOplus_DR_CDF for given altitude range using imported density/temp profiles
void Distribution_Hot_H::make_HCOplus_DR_CDF(double lower_alt, double upper_alt)
{
	//ofstream outfile;
	//outfile.open("/home/rodney/Documents/coronaTest/rodney_hcopl_lsa.dat");

	double bin_size = 10000.0; // [cm]
	int num_alt_bins = (int)((upper_alt - lower_alt) / bin_size);
	vector<double> HCOplus_DR_rate;
	HCOplus_DR_rate.resize(num_alt_bins);
	HCOplus_DR_CDF[0].resize(num_alt_bins);
	HCOplus_DR_CDF[1].resize(num_alt_bins);

	double bottom_HCOplus = HCOplus_profile[0][0];
	double top_HCOplus = HCOplus_profile[0].back();
	double bottom_e = electron_profile[0][0];
	double top_e = electron_profile[0].back();

	double local_g = (constants::G * my_planet.get_mass()) / (pow(my_planet.get_radius()+bottom_HCOplus, 2.0));
	double HCOplus_bottom_scaleheight = constants::k_b*common::interpolate_logy(temp_profile[0], temp_profile[2], bottom_HCOplus)/(m_HCOplus*local_g);
	local_g = (constants::G * my_planet.get_mass()) / (pow(my_planet.get_radius()+top_HCOplus, 2.0));
	double HCOplus_top_scaleheight = constants::k_b*common::interpolate_logy(temp_profile[0], temp_profile[2], top_HCOplus)/(m_HCOplus*local_g);
	local_g = (constants::G * my_planet.get_mass()) / (pow(my_planet.get_radius()+bottom_e, 2.0));
	double e_bottom_scaleheight = constants::k_b*common::interpolate_logy(temp_profile[0], temp_profile[3], bottom_e)/(constants::m_e*local_g);
	local_g = (constants::G * my_planet.get_mass()) / (pow(my_planet.get_radius()+top_e, 2.0));
	double e_top_scaleheight = constants::k_b*common::interpolate_logy(temp_profile[0], temp_profile[3], top_e)/(constants::m_e*local_g);

	double rate_sum = 0.0;
	double rate_sum_times_r_sqrd = 0.0;
	for (int i=0; i<num_alt_bins; i++)
	{
		HCOplus_DR_CDF[1][i] = lower_alt + bin_size*i;

		// get new e density by either interpolation or extrapolation
		double e_dens = 0.0;
		if (HCOplus_DR_CDF[1][i] < bottom_e)
		{
			e_dens = electron_profile[1][0]*exp((bottom_e - HCOplus_DR_CDF[1][i])/e_bottom_scaleheight);
		}
		else if (HCOplus_DR_CDF[1][i] > top_e)
		{
			e_dens = electron_profile[1].back()*exp((top_e - HCOplus_DR_CDF[1][i])/e_top_scaleheight);
		}
		else
		{
			e_dens = common::interpolate_logy(electron_profile[0], electron_profile[1], HCOplus_DR_CDF[1][i]);
		}

		// get new HCOplus density by either interpolation or extrapolation
		double HCOplus_dens = 0.0;
		if (HCOplus_DR_CDF[1][i] < bottom_HCOplus)
		{
			HCOplus_dens = HCOplus_profile[1][0]*exp((bottom_HCOplus - HCOplus_DR_CDF[1][i])/HCOplus_bottom_scaleheight);
		}
		else if (HCOplus_DR_CDF[1][i] > top_HCOplus)
		{
			HCOplus_dens = HCOplus_profile[1].back()*exp((top_HCOplus - HCOplus_DR_CDF[1][i])/HCOplus_top_scaleheight);
		}
		else
		{
			HCOplus_dens = common::interpolate_logy(HCOplus_profile[0], HCOplus_profile[1], HCOplus_DR_CDF[1][i]);
		}

		// calculate HCOplus_DR_rate at current alt using rate coefficient from Fox 2015
		double Te = common::interpolate_logy(temp_profile[0], temp_profile[3], HCOplus_DR_CDF[1][i]);
		if (Te <= 300.0)
		{
			HCOplus_DR_rate[i] = HCOplus_DR_rate_coeff * pow((Te/300.0), -1.25) * e_dens * HCOplus_dens;
		}
		else
		{
			HCOplus_DR_rate[i] = HCOplus_DR_rate_coeff * pow((Te/300.0), -1.0) * e_dens * HCOplus_dens;
		}
		//outfile << HCOplus_DR_CDF[1][i]*1e-5 << "\t" << HCOplus_DR_rate[i] << "\n";
		rate_sum = rate_sum + HCOplus_DR_rate[i];
		rate_sum_times_r_sqrd = rate_sum_times_r_sqrd + (HCOplus_DR_rate[i] * 4.0 * constants::pi * pow(my_planet.get_radius()+(lower_alt+(bin_size*i)), 2.0));
	}
	//outfile.close();

	for (int i=0; i<num_alt_bins; i++)
	{
		if (i == 0)
		{
			HCOplus_DR_CDF[0][i] = HCOplus_DR_rate[i] / rate_sum;
		}
		else
		{
			HCOplus_DR_CDF[0][i] = (HCOplus_DR_rate[i] / rate_sum) + HCOplus_DR_CDF[0][i-1];
		}
	}
	//double global_rate_HCOplus_DR = rate_sum*bin_size * 4.0 * constants::pi * pow(my_planet.get_radius()+upper_alt, 2.0);
	double global_rate_HCOplus_DR = rate_sum_times_r_sqrd*bin_size;
	cout << "Global hot H production rate from HCO+ DR:\n" << global_rate_HCOplus_DR << " per second\n";
	global_rate = global_rate_HCOplus_DR;
}

// generate any_mechanism_prob_CDF for given altitude range - particular production rates are ignored, and all test particles are produced in a certain altitude bin
void Distribution_Hot_H::make_any_mechanism_prob_CDF(double lower_alt, double upper_alt)
{
	double bin_size = 10000.0; // [cm]
	int num_alt_bins = (int)((upper_alt - lower_alt) / bin_size);

	any_mechanism_prob_CDF[0].resize(num_alt_bins);
	any_mechanism_prob_CDF[1].resize(num_alt_bins);

	for (int i=0; i<num_alt_bins; i++)
	  		
	{
	        any_mechanism_prob_CDF[1][i] = lower_alt + bin_size*i;
		// the section below produces all test particles in a single bin
                if (i == any_mechanism_alt_bin) // any_mechanism_alt_bin is the altitude bin at which particles are produced (set in Hot_H.cfg)
		{
		        any_mechanism_prob_CDF[0][i] = 1;
		}
		else
		{
		        any_mechanism_prob_CDF[0][i] = any_mechanism_prob_CDF[0][i-1] + 0.0; 
		}
	}
}
