/*
 * Distribution_Hot_O.cpp
 *
 *  Created on: Jul 22, 2020
 *      Author: rodney
 */

#include "Distribution_Hot_O.hpp"

Distribution_Hot_O::Distribution_Hot_O(Planet my_p, double ref_h, double ref_T)
	: Distribution(my_p, ref_h, ref_T) {

	// these 4 are from Justin's original Fortran code; may still need at some point
	DR160 = 600.0;
	H_DR = 19e5;
	T_ion = 400.0;
	T_e = 1600.0;

	m_O = 15.9994*constants::amu;
	m_O2plus = 31.9983*constants::amu;
	source = "";
	O2plus_DR_rate_coeff = 0.0;
	global_rate = 0.0;

	temp_profile.resize(4);
	O2plus_profile.resize(2);
	electron_profile.resize(2);
	O2plus_DR_CDF.resize(2);

	// import parameters from configuration file
	// file must be named 'Hot_O.cfg' and be in same directory as executable
	double profile_bottom = 0.0;
	double profile_top = 0.0;
	string temp_prof_filename = "";
	string O2plus_prof_filename = "";
	string electron_prof_filename = "";

	ifstream infile;
	infile.open("Hot_O.cfg");
	if (!infile.good())
	{
		cout << "Hot O configuration file not found!\n";
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
		if (parameters[i] == "O2plus_DR_rate_coeff")
		{
			O2plus_DR_rate_coeff = stod(values[i]);
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
		else if (parameters[i] == "O2plus_prof_filename")
		{
			O2plus_prof_filename = values[i];
		}
		else if (parameters[i] == "electron_prof_filename")
		{
			electron_prof_filename = values[i];
		}
	}

	common::import_csv(temp_prof_filename, temp_profile[0], temp_profile[1], temp_profile[2], temp_profile[3]);
	common::import_csv(O2plus_prof_filename, O2plus_profile[0], O2plus_profile[1]);
	common::import_csv(electron_prof_filename, electron_profile[0], electron_profile[1]);

	if (source == "O2plus_DR")
	{
		make_O2plus_DR_CDF(profile_bottom, profile_top);
	}
	else if (source == "O2plus_DR_old_method")
	{
		global_rate = H_DR*DR160*4.0*constants::pi*pow((3557e5+H_DR), 2.0);
	}
}

Distribution_Hot_O::~Distribution_Hot_O() {

}

void Distribution_Hot_O::init(shared_ptr<Particle> p)
{
	if (source == "O2plus_DR")
	{
		init_O2plus_DR_particle(p);
	}
	else
	{
		init_old_way(p);
	}
}

// init particle using O2plus_DR mechanism
void Distribution_Hot_O::init_O2plus_DR_particle(shared_ptr<Particle> p)
{
	// altitude distribution for O2+ dissociative recombination
	double r = get_new_radius_O2plus_DR();
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

	// Select Hot O Electronic Channel ###
	// Probabilities and energies
	// Particles are activated per channel to allow testing the contribution of each population in the simulation
	double Ei = 0.0;
	double randnum = common::get_rand();
	if (randnum < 0.22)
	{
		Ei = 6.99*constants::ergev;
	}
	else if (randnum < (0.22+0.42))
	{
		Ei = 5.02*constants::ergev;
	}
	else if (randnum < (0.22+0.42+0.31))
	{
		Ei = 3.06*constants::ergev;
	}
	else
	{
		Ei = 0.84*constants::ergev;
	}

	// Thermal Rotational Energy of O2+
	Ei = Ei + E_rot(3.35967e-16, temp_ion);

	// Translational Energy per O Resulting from Dissociative Recombination of O2+
	double v = sqrt(Ei / p->get_mass());

	// spherically isotropic velocity vector
	phi = constants::twopi*common::get_rand();
	u = 2.0*common::get_rand() - 1.0;
	double vx = v*sqrt(1-u*u)*cos(phi);
	double vy = v*sqrt(1-u*u)*sin(phi);
	double vz = v*u;

	// Add initial ion and electron translational momentum
    double vavg = sqrt(constants::k_b*temp_ion/(m_O2plus));  // average thermal ion velocity
    double v_ion[] = {0.0, 0.0, 0.0};
    gen_mb(vavg, v_ion);                               // sample from Maxwell-Boltzmann distribution
    vavg = sqrt(constants::k_b*temp_e/constants::m_e);    //average thermal electron velocity
    double v_e[] = {0.0, 0.0, 0.0};
    gen_mb(vavg, v_e);                                 // sample from Maxwell-Boltzmann distribution

    // add it all up
    vx = vx + (m_O2plus*v_ion[0] + constants::m_e*v_e[0]) / (m_O2plus+constants::m_e);
    vy = vy + (m_O2plus*v_ion[1] + constants::m_e*v_e[1]) / (m_O2plus+constants::m_e);
    vz = vz + (m_O2plus*v_ion[2] + constants::m_e*v_e[2]) / (m_O2plus+constants::m_e);

    p->init_particle(x, y, z, vx, vy, vz);
}

// init particle using Justin's original method
void Distribution_Hot_O::init_old_way(shared_ptr<Particle> p)
{
	// altitude distribution for O2+ dissociative recombination
	double r = my_planet.get_radius() + 160e5 - log(common::get_rand())*H_DR;

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

	// Select Hot O Electronic Channel ###
	// Probabilities and energies
	// Particles are activated per channel to allow testing the contribution of each population in the simulation
	double Ei = 0.0;
	double randnum = common::get_rand();
	if (randnum < 0.22)
	{
		Ei = 6.99*constants::ergev;
	}
	else if (randnum < (0.22+0.42))
	{
		Ei = 5.02*constants::ergev;
	}
	else if (randnum < (0.22+0.42+0.31))
	{
		Ei = 3.06*constants::ergev;
	}
	else
	{
		Ei = 0.84*constants::ergev;
	}

	// Thermal Rotational Energy of O2+
	Ei = Ei + E_rot(3.35967e-16, T_ion);

	// Translational Energy per O Resulting from Dissociative Recombination of O2+
	double v = sqrt(Ei / p->get_mass());

	// spherically isotropic velocity vector
	phi = constants::twopi*common::get_rand();
	u = 2.0*common::get_rand() - 1.0;
	double vx = v*sqrt(1-u*u)*cos(phi);
	double vy = v*sqrt(1-u*u)*sin(phi);
	double vz = v*u;

	// Add initial ion and electron translational momentum
    double vavg = sqrt(constants::k_b*T_ion/(m_O2plus));  // average thermal ion velocity
    double v_ion[] = {0.0, 0.0, 0.0};
    gen_mb(vavg, v_ion);                               // sample from Maxwell-Boltzmann distribution
    vavg = sqrt(constants::k_b*T_e/constants::m_e);    //average thermal electron velocity
    double v_e[] = {0.0, 0.0, 0.0};
    gen_mb(vavg, v_e);                                 // sample from Maxwell-Boltzmann distribution

    // add it all up
    vx = vx + (m_O2plus*v_ion[0] + constants::m_e*v_e[0]) / (m_O2plus+constants::m_e);
    vy = vy + (m_O2plus*v_ion[1] + constants::m_e*v_e[1]) / (m_O2plus+constants::m_e);
    vz = vz + (m_O2plus*v_ion[2] + constants::m_e*v_e[2]) / (m_O2plus+constants::m_e);

    p->init_particle(x, y, z, vx, vy, vz);
}

// scans O2plus_DR_CDF for new particle radius
double Distribution_Hot_O::get_new_radius_O2plus_DR()
{
	double u = common::get_rand();
	int k = 0;
	while (O2plus_DR_CDF[0][k] < u)
	{
		k++;
	}
	return O2plus_DR_CDF[1][k] + my_planet.get_radius();
}

// Sample Rotational Energy from a Thermal Population of Rigid Rotators
double Distribution_Hot_O::E_rot(double B, double T)
{
	int j = 0;
	int nj = 200;  // highest rotational level considered

	// partition function using two terms of Euler-MacLaurin summation formula
	double Q_rot = constants::k_b*T/B + 1.0/3.0;

	double P_tot = 0.0;
	double P_j = 0.0;
	double u = common::get_rand();

	for (int i=0; i<(nj+1); i++)
	{
		P_j = (2.0*j+1.0)*exp(-j*(j+1)*B/(constants::k_b*T))/Q_rot;
		P_tot = P_tot + P_j;
		j = i;
		if (u < P_tot)
		{
			break;
		}
	}
	return j*(j+1.0)*B;
}

// returns global production rate (needs to be set by chosen production method)
// calling function must divide this by 2 if only hemispherical rate is needed
double Distribution_Hot_O::get_global_rate()
{
	return global_rate;
}

// generate O2plus_DR_CDF for given altitude range using imported density/temp profiles
void Distribution_Hot_O::make_O2plus_DR_CDF(double lower_alt, double upper_alt)
{
	ofstream outfile;
	outfile.open("./o2plus_dr_rates.dat");

	double bin_size = 10000.0; // [cm]
	int num_alt_bins = (int)((upper_alt - lower_alt) / bin_size);
	vector<double> O2plus_DR_rate;
	O2plus_DR_rate.resize(num_alt_bins);
	O2plus_DR_CDF[0].resize(num_alt_bins);
	O2plus_DR_CDF[1].resize(num_alt_bins);

	// get top and bottom densities from profiles
	double bottom_O2plus = O2plus_profile[0][0];
	double top_O2plus = O2plus_profile[0].back();
	double bottom_e = electron_profile[0][0];
	double top_e = electron_profile[0].back();

	// set scale heights for densities at top and bottom boundaries of profiles
	double local_g = (constants::G * my_planet.get_mass()) / (pow(my_planet.get_radius()+bottom_O2plus, 2.0));
	double O2plus_bottom_scaleheight = constants::k_b*common::interpolate_logy(temp_profile[0], temp_profile[2], bottom_O2plus)/(m_O2plus*local_g);
	local_g = (constants::G * my_planet.get_mass()) / (pow(my_planet.get_radius()+top_O2plus, 2.0));
	double O2plus_top_scaleheight = constants::k_b*common::interpolate_logy(temp_profile[0], temp_profile[2], top_O2plus)/(m_O2plus*local_g);
	local_g = (constants::G * my_planet.get_mass()) / (pow(my_planet.get_radius()+bottom_e, 2.0));
	double e_bottom_scaleheight = constants::k_b*common::interpolate_logy(temp_profile[0], temp_profile[3], bottom_e)/(constants::m_e*local_g);
	local_g = (constants::G * my_planet.get_mass()) / (pow(my_planet.get_radius()+top_e, 2.0));
	double e_top_scaleheight = constants::k_b*common::interpolate_logy(temp_profile[0], temp_profile[3], top_e)/(constants::m_e*local_g);

	double rate_sum = 0.0;
	double rate_sum_times_r_sqrd = 0.0;
	for (int i=0; i<num_alt_bins; i++)
	{
		O2plus_DR_CDF[1][i] = lower_alt + bin_size*i;

		// get new e density by either interpolation or extrapolation
		double e_dens = 0.0;
		if (O2plus_DR_CDF[1][i] < bottom_e)
		{
			e_dens = electron_profile[1][0]*exp((bottom_e - O2plus_DR_CDF[1][i])/e_bottom_scaleheight);
		}
		else if (O2plus_DR_CDF[1][i] > top_e)
		{
			e_dens = electron_profile[1].back()*exp((top_e - O2plus_DR_CDF[1][i])/e_top_scaleheight);
		}
		else
		{
			e_dens = common::interpolate_logy(electron_profile[0], electron_profile[1], O2plus_DR_CDF[1][i]);
		}

		// get new O2plus density by either interpolation or extrapolation
		double O2plus_dens = 0.0;
		if (O2plus_DR_CDF[1][i] < bottom_O2plus)
		{
			O2plus_dens = O2plus_profile[1][0]*exp((bottom_O2plus - O2plus_DR_CDF[1][i])/O2plus_bottom_scaleheight);
		}
		else if (O2plus_DR_CDF[1][i] > top_O2plus)
		{
			O2plus_dens = O2plus_profile[1].back()*exp((top_O2plus - O2plus_DR_CDF[1][i])/O2plus_top_scaleheight);
		}
		else
		{
			O2plus_dens = common::interpolate_logy(O2plus_profile[0], O2plus_profile[1], O2plus_DR_CDF[1][i]);
		}

		// calculate O2plus_DR_rate at current alt using rate coefficient from configuration file
		double Te = common::interpolate_logy(temp_profile[0], temp_profile[3], O2plus_DR_CDF[1][i]);

		O2plus_DR_rate[i] = 2.0 * O2plus_DR_rate_coeff * pow((300.0/Te), 0.7) * e_dens * O2plus_dens;

		outfile << O2plus_DR_CDF[1][i]*1e-5 << "\t" << O2plus_DR_rate[i] << "\n";
		rate_sum = rate_sum + O2plus_DR_rate[i];
		rate_sum_times_r_sqrd = rate_sum_times_r_sqrd + (O2plus_DR_rate[i] * 4.0 * constants::pi * pow(my_planet.get_radius()+(lower_alt+(bin_size*i)), 2.0));
	}
	outfile.close();

	for (int i=0; i<num_alt_bins; i++)
	{
		if (i == 0)
		{
			O2plus_DR_CDF[0][i] = O2plus_DR_rate[i] / rate_sum;
		}
		else
		{
			O2plus_DR_CDF[0][i] = (O2plus_DR_rate[i] / rate_sum) + O2plus_DR_CDF[0][i-1];
		}
	}

	double global_rate_O2plus_DR = rate_sum_times_r_sqrd*bin_size;
	cout << "Global hot O production rate from O2+ DR:\n" << global_rate_O2plus_DR << " per second\n";
	global_rate = global_rate_O2plus_DR;
}


