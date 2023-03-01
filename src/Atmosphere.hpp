/*
 * Atmosphere.hpp
 *
 *  Created on: Jun 8, 2020
 *      Author: rodney
 */

#ifndef ATMOSPHERE_HPP_
#define ATMOSPHERE_HPP_

#include <vector>
#include <iomanip>
#include "Background_Species.hpp"
#include "Distribution_Hot_H.hpp"
#include "Distribution_Hot_O.hpp"
#include "Distribution_Import.hpp"
#include "Distribution_MB.hpp"
#include "Common_Functions.hpp"
using namespace std;

class Atmosphere {
public:
	Atmosphere(int n, int num_to_trace, string trace_output_dir, Planet p, vector<shared_ptr<Particle>> parts, shared_ptr<Distribution> dist, Background_Species bg, int num_EDFs, int EDF_alts[]);
	virtual ~Atmosphere();

	void output_positions(string datapath);
	void output_altitude_distro(double bin_width, string datapath);
	void output_velocity_distro(double bin_width, string datapath);
	void output_alt_energy_distro(double alt_in_cm, double e_bin_width, string datapath);
	void run_simulation(double dt, int num_steps, double lower_bound, double upper_bound, int print_status_freq, int output_pos_freq, string output_pos_dir, string output_stats_dir);

private:
	int num_parts;                      // number of particles initially spawned
	int num_traced;                     // number of particles to output trace data on
	string trace_dir;                   // directory to output particle trace data to
	int active_parts;                   // number of active particles
	Planet my_planet;                   // contains planet mass and radius
	vector<shared_ptr<Particle>> my_parts;         // array of particles to be tracked
	shared_ptr<Distribution> my_dist;              // distribution class to initialize particles
	Background_Species bg_species;      // background species used for collisions
	vector<int> traced_parts;           // indices of randomly selected trace particles

	vector<vector<int>> stats_dens_counts;  // vector for accumulating particle density counts
	vector<int> stats_coldens_counts;  // vector for accumulating integrated dayside column density counts
        vector<double> stats_angleavg_dens;  // vector for accumulating angle-averaged column density counts in x=const. plane
	vector<vector<int>> stats_dens2d_counts;  // stores a 2d grid of dayside column density counts
	int stats_num_EDFs;  // number of altitude EDFs to track; populated from corona3d_2020.cfg
	vector<int> stats_EDF_alts;  // holds list of altitudes that (in km above surface) that EDFs are tracked at
	vector<vector<vector<vector<double>>>> stats_EDFs;  // EDF counts are accumulated here
	vector<double> stats_loss_rates;  // loss rates at each EDF altitude are calculated and stored here

	// these two modules are where stats are accumulated and then output at the end of a simulation
	void update_stats(double dt, int idx);
	void output_stats(double dt, double rate, int total_parts, string output_dir);

	// output test particle trace data for selected particles
	void output_collision_data();
	void output_trace_data();
};

#endif /* ATMOSPHERE_HPP_ */
