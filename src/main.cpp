/*
 * main.cpp
 *
 *  Created on: May 27, 2020
 *      Author: rodney
 */
#include <iostream>
#include <fstream>
#include <iomanip>
#include <typeinfo>
#include "Atmosphere.hpp"
using namespace std;

//subroutine to set particle types from main
shared_ptr<Particle> set_particle_type(string type)
{
	shared_ptr<Particle> p;

	if (type == "H")
	{
		p = make_shared<Particle_H>();
	}
	else if (type == "O")
	{
		p = make_shared<Particle_O>();
	}
	else if (type == "N2")
	{
		p = make_shared<Particle_N2>();
	}
	else if (type == "CO")
	{
		p = make_shared<Particle_CO>();
	}
	else if (type == "CO2")
	{
		p = make_shared<Particle_CO2>();
	}
	else
	{
		cout << "Invalid particle type specified! Please check configuration file.\n";
		exit(1);
	}
	return p;
}

int main(int argc, char* argv[])
{
	cout << "Initializing Simulation...\n";

	//initialize and read parameters from configuration file
	int num_testparts = 0;
	string part_type = "";
	string dist_type = "";
	string pos_infile = "";
	string vel_infile = "";
	string output_dir = "";
	int num_traced = 0;
	string trace_output_dir = "";
	int num_EDFs = 0;
	int EDF_alts_index = 0;
	int print_status_freq = 0;
	int output_pos_freq = 0;
	string output_pos_dir = "";
	string output_stats_dir = "";
	double profile_bottom_alt = 0.0;
	double profile_top_alt = 0.0;
	string temp_profile_filename = "";
	string neut_densities_filename = "";
	string ion_densities_filename = "";
	int timesteps = 0;
	double dt = 0.0;
	double ref_height = 0.0;
	double ref_temp = 0.0;
	double planet_mass = 0.0;
	double planet_radius = 0.0;
	double sim_lower_bound = 0.0;
	double sim_upper_bound = 0.0;
	Planet my_planet;
	vector<shared_ptr<Particle>> parts;
	shared_ptr<Distribution> dist;
	int num_bgparts = 0;
	int bg_params_index = 0;

	ifstream infile;
	infile.open("corona3d_2020.cfg");
	if (!infile.good())
	{
		cout << "Configuration file not found!\n";
		return 1;
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
		if (parameters[i] == "num_testparts")
		{
			num_testparts = stoi(values[i]);
		}
		else if (parameters[i] == "part_type")
		{
			part_type = values[i];
		}
		else if (parameters[i] == "dist_type")
		{
			dist_type = values[i];
		}
		else if (parameters[i] == "pos_infile")
		{
			pos_infile = values[i];
		}
		else if (parameters[i] == "vel_infle")
		{
			vel_infile = values[i];
		}
		else if (parameters[i] == "output_dir")
		{
			output_dir = values[i];
		}
		else if (parameters[i] == "num_traced")
		{
			num_traced = stoi(values[i]);
		}
		else if (parameters[i] == "trace_output_dir")
		{
			trace_output_dir = values[i];
		}
		else if (parameters[i] == "print_status_freq")
		{
			print_status_freq = stoi(values[i]);
		}
		else if (parameters[i] == "output_pos_freq")
		{
			output_pos_freq = stoi(values[i]);
		}
		else if (parameters[i] == "output_pos_dir")
		{
			output_pos_dir = values[i];
		}
		else if (parameters[i] == "output_stats_dir")
		{
			output_stats_dir = values[i];
		}
		else if (parameters[i] == "profile_bottom_alt")
		{
			profile_bottom_alt = stod(values[i]);
		}
		else if (parameters[i] == "profile_top_alt")
		{
			profile_top_alt = stod(values[i]);
		}
		else if (parameters[i] == "temp_profile")
		{
			temp_profile_filename = values[i];
		}
		else if (parameters[i] == "neutral_densities")
		{
			neut_densities_filename = values[i];
		}
		else if (parameters[i] == "ion_densities")
		{
			ion_densities_filename = values[i];
		}
		else if (parameters[i] == "timesteps")
		{
			timesteps = stoi(values[i]);
		}
		else if (parameters[i] == "dt")
		{
			dt = stod(values[i]);
		}
		else if (parameters[i] == "ref_height")
		{
			ref_height = stod(values[i]);
		}
		else if (parameters[i] == "ref_temp")
		{
			ref_temp = stod(values[i]);
		}
		else if (parameters[i] == "planet_mass")
		{
			planet_mass = stod(values[i]);
		}
		else if (parameters[i] == "planet_radius")
		{
			planet_radius = stod(values[i]);
		}
		else if (parameters[i] == "sim_lower_bound")
		{
			sim_lower_bound = stod(values[i]);
		}
		else if (parameters[i] == "sim_upper_bound")
		{
			sim_upper_bound = stod(values[i]);
		}
		else if (parameters[i] == "num_bgparts")
		{
			num_bgparts = stoi(values[i]);
			bg_params_index = i+1;
		}
		else if (parameters[i] == "num_EDFs")
		{
			num_EDFs = stoi(values[i]);
			EDF_alts_index = i+1;
		}
	}

	//initialize planet and test particles
	my_planet.init(planet_mass, planet_radius);
	parts.resize(num_testparts);
	for (int i=0; i<num_testparts; i++)
	{
		parts[i] = set_particle_type(part_type);
	}

	//instantiate the Distribution class to be used
	if (dist_type == "Hot_H")
	{
		dist = make_shared<Distribution_Hot_H>(my_planet, ref_height, ref_temp);
	}
	else if (dist_type == "Hot_O")
	{
		dist = make_shared<Distribution_Hot_O>(my_planet, ref_height, ref_temp);
	}
	else if (dist_type == "MB")
	{
		dist = make_shared<Distribution_MB>(my_planet, ref_height, ref_temp);
	}
	else if (dist_type == "Import")
	{
		dist = make_shared<Distribution_Import>(my_planet, ref_height, ref_temp, pos_infile, vel_infile);
	}
	else
	{
		cout << "Invalid distribution type! Please check configuration file.\n";
		return 1;
	}

	//initialize Background_Species
	string bg_config_files[num_bgparts];
	for (int i=0; i<num_bgparts; i++)
	{
		bg_config_files[i] = values[bg_params_index + i];
	}
	Background_Species bg_spec(num_bgparts, bg_config_files, my_planet, ref_temp, ref_height, temp_profile_filename, neut_densities_filename, profile_bottom_alt, profile_top_alt);

	//set up EDF altitudes to be passed to atmosphere class
	int EDF_alts[num_EDFs];
	for (int i=0; i<num_EDFs; i++)
	{
		EDF_alts[i] = stoi(values[EDF_alts_index + i]);
	}

	//set output directory paths to default if necessary
	if (output_pos_dir == "")
	{
		output_pos_dir = output_dir;
	}
	if (trace_output_dir == "")
	{
		trace_output_dir = output_dir;
	}
	if (output_stats_dir == "")
	{
		output_stats_dir = output_dir;
	}

	// initialize atmosphere and run simulation
	Atmosphere my_atmosphere(num_testparts, num_traced, trace_output_dir, my_planet, parts, dist, bg_spec, num_EDFs, EDF_alts);
	//my_atmosphere.output_velocity_distro(10000.0, output_dir + "vdist.out");
	//my_atmosphere.output_altitude_distro(100000.0, output_dir + "altdist.out");
	//my_atmosphere.output_alt_energy_distro(133e5, 0.03, output_dir + "edist.out");
	my_atmosphere.run_simulation(dt, timesteps, sim_lower_bound, sim_upper_bound, print_status_freq, output_pos_freq, output_pos_dir, output_stats_dir);
	//my_atmosphere.output_velocity_distro(10000.0, output_dir + "vdist2.out");
	//my_atmosphere.output_altitude_distro(100000.0, output_dir + "altdist2.out");

	return 0;
}
