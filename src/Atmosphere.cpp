/*
 * Atmosphere.cpp
 *
 *  Created on: Jun 8, 2020
 *      Author: rodney
 */

#include "Atmosphere.hpp"

// construct atmosphere using given parameters
Atmosphere::Atmosphere(int n, int num_to_trace, string trace_output_dir, Planet p, vector<shared_ptr<Particle>> parts, shared_ptr<Distribution> dist, Background_Species bg, int num_EDFs, int EDF_alts[])
{
	num_parts = n;                // number of test particles to track
	num_traced = num_to_trace;    // number of tracked particles to output detailed trace data for
	trace_dir = trace_output_dir;
	active_parts = num_parts;
	my_planet = p;
	my_dist = dist;
	my_parts.resize(num_parts);
	bg_species = bg;

	for (int i=0; i<num_parts; i++)
	{
		my_parts[i] = parts[i];
		my_dist->init(my_parts[i]);
	}

	//initialize stats tracking vectors
	stats_num_EDFs = num_EDFs;
	stats_EDF_alts.resize(stats_num_EDFs);
	stats_loss_rates.resize(stats_num_EDFs);
	stats_EDFs.resize(2);   // index 0 is day side EDFs, 1 is night side
	stats_EDFs[0].resize(stats_num_EDFs);
	stats_EDFs[1].resize(stats_num_EDFs);
	stats_angleavg_dens.resize(stats_num_EDFs); // vector for accumulating angle-averaged column density counts in x=const. plane

	for (int i=0; i<stats_num_EDFs; i++)
	{
		stats_EDF_alts[i] = EDF_alts[i];
		stats_loss_rates[i] = 0.0;
		stats_EDFs[0][i].resize(201);
		stats_EDFs[1][i].resize(201);

		for (int j=0; j<201; j++)
		{
			stats_EDFs[0][i][j].resize(201);
			stats_EDFs[1][i][j].resize(201);

			for(int k=0; k<201; k++)
			{
				stats_EDFs[0][i][j][k] = 0.0;
				stats_EDFs[1][i][j][k] = 0.0;
			}
		}
	}

	stats_dens_counts.resize(2);  // index 0 is day side, 1 is night side
	stats_dens_counts[0].resize(100001);
	stats_dens_counts[1].resize(100001);
	stats_coldens_counts.resize(100001);
	for (int i=0; i<100001; i++)
	{
		stats_dens_counts[0][i] = 0;
		stats_dens_counts[1][i] = 0;
		stats_coldens_counts[i] = 0;
	}

	stats_dens2d_counts.resize(1025);
	for (int i=0; i<1025; i++)
	{
		stats_dens2d_counts[i].resize(1025);
		for (int j=0; j<1025; j++)
		{
			stats_dens2d_counts[i][j] = 0;
		}
	}

	// pick trace particles if any
	if (num_traced > 0)
	{
		traced_parts.resize(num_traced);
		for (int i=0; i<num_traced; i++)
		{
			traced_parts[i] = common::get_rand_int(0, num_parts-1);
			my_parts[traced_parts[i]]->set_traced();
		}
	}
}

Atmosphere::~Atmosphere() {

}

// writes single-column output file of altitude bin counts using active particles
// bin_width is in cm; first 2 numbers in output file are bin_width and num_bins
void Atmosphere::output_altitude_distro(double bin_width, string datapath)
{
	ofstream outfile;
	outfile.open(datapath);
	double alt = 0.0;           // altitude of particle [cm]
	int nb = 0;                 // bin number

	double max_radius = 0.0;
	for (int i=0; i<num_parts; i++)
	{
		if (my_parts[i]->is_active())
		{
			if (my_parts[i]->get_radius() > max_radius)
			{
				max_radius = my_parts[i]->get_radius();
			}
		}
	}
	int num_bins = (int)((max_radius - my_planet.get_radius()) / bin_width) + 10;
	int abins[num_bins] = {0};  // array of altitude bin counts

	for (int i=0; i<num_parts; i++)
	{
		if (my_parts[i]->is_active())
		{
			alt = my_parts[i]->get_radius() - my_planet.get_radius();
			nb = (int)(alt / bin_width);
			abins[nb]++;
		}
	}

	outfile << bin_width << "\n";
	outfile << num_bins << "\n";

	for (int i=0; i<num_bins; i++)
	{
		outfile << abins[i] << "\n";
	}
	outfile.close();
}

void Atmosphere::output_collision_data()
{
	for (int i=0; i<num_traced; i++)
	{
		string filename = trace_dir + "part" + to_string(traced_parts[i]) + "_collisions.out";
		my_parts[traced_parts[i]]->dump_collision_log(filename);
	}
}

// writes 3-column output file of all current particle positions
// file is saved to location specified by datapath
void Atmosphere::output_positions(string datapath)
{
	ofstream outfile;
	outfile.open(datapath);
	for (int i=0; i<num_parts; i++)
	{
		outfile << setprecision(10) << my_parts[i]->get_x() << '\t';
		outfile << setprecision(10) << my_parts[i]->get_y() << '\t';
		outfile << setprecision(10) << my_parts[i]->get_z() << '\n';
	}
	outfile.close();
}

// output test particle trace data for selected particles
void Atmosphere::output_trace_data()
{
	for (int i=0; i<num_traced; i++)
	{
		if (my_parts[traced_parts[i]]->is_active())
		{
			ofstream position_file;
			position_file.open(trace_dir + "part" + to_string(traced_parts[i]) + "_positions.out", ios::out | ios::app);
			position_file << setprecision(10) << my_parts[traced_parts[i]]->get_x() << '\t';
			position_file << setprecision(10) << my_parts[traced_parts[i]]->get_y() << '\t';
			position_file << setprecision(10) << my_parts[traced_parts[i]]->get_z() << '\n';
			position_file.close();
		}
	}
}

// writes single-column output file of velocity bin counts using active particles
// bin_width is in cm/s; first 2 numbers in output file are bin_width and num_bins
void Atmosphere::output_velocity_distro(double bin_width, string datapath)
{
	ofstream outfile;
	outfile.open(datapath);
	double v = 0.0;             // velocity magnitude [cm/s]
	int nb = 0;                 // bin number

	double max_v = 0.0;
	for (int i=0; i<num_parts; i++)
	{
		if (my_parts[i]->is_active())
		{
			double total_v = my_parts[i]->get_total_v();
			if (total_v > max_v)
			{
				max_v = total_v;
			}
		}
	}
	int num_bins = (int)((max_v / bin_width) + 10);
	int vbins[num_bins] = {0};  // array of velocity bin counts

	for (int i=0; i<num_parts; i++)
	{
		if (my_parts[i]->is_active())
		{
			v = my_parts[i]->get_total_v();
			nb = (int)(v / bin_width);
			vbins[nb]++;
		}
	}

	outfile << bin_width << "\n";
	outfile << num_bins << "\n";

	for (int i=0; i<num_bins; i++)
	{
		outfile << vbins[i] << '\n';
	}
	outfile.close();
}

// writes single-column output file of energy bin counts using active particles at given altitude
// bin_width is in eV; first 2 numbers in output file are bin_width and num_bins
void Atmosphere::output_alt_energy_distro(double alt_in_cm, double e_bin_width, string datapath)
{
	double r = my_planet.get_radius() + alt_in_cm;
	ofstream outfile;
	outfile.open(datapath);
	double e = 0.0;             // particle energy [eV]
	int nb = 0;                 // bin number

	double max_e = 0.0;
	for (int i=0; i<num_parts; i++)
	{
		if (my_parts[i]->is_active() && my_parts[i]->get_radius() >= r && my_parts[i]->get_radius() < r + 1e5)
		{
			double total_e = my_parts[i]->get_energy_in_eV();
			if (total_e > max_e)
			{
				max_e = total_e;
			}
		}
	}
	int num_bins = (int)((max_e / e_bin_width) + 10);
	int ebins[num_bins] = {0};  // array of energy bin counts

	for (int i=0; i<num_parts; i++)
	{
		if (my_parts[i]->is_active() && my_parts[i]->get_radius() >= r && my_parts[i]->get_radius() < r + 1e5)
		{
			e = my_parts[i]->get_energy_in_eV();
			nb = (int)(e / e_bin_width);
			ebins[nb]++;
		}
	}

	outfile << e_bin_width << "\n";
	outfile << num_bins << "\n";

	for (int i=0; i<num_bins; i++)
	{
		outfile << ebins[i] << '\n';
	}
	outfile.close();
}

// iterate equation of motion and check for collisions for each active particle being tracked
// a lot of stuff in here needs to be changed to be dynamically determined at runtime
void Atmosphere::run_simulation(double dt, int num_steps, double lower_bound, double upper_bound, int print_status_freq, int output_pos_freq, string output_pos_dir, string output_stats_dir)
{
	int night_escape_count = 0;
	int day_escape_count = 0;
	double v_esc_current = 0.0;
	double v_thermal = 0.0;

	// most probable MB velocity of test particle at 200K
	//double v_mp = sqrt(2.0*constants::k_b*200.0/my_parts[0]->get_mass());

	// RMS thermal velocity of test particle at 200K
	//double v_rms = sqrt(3.0*constants::k_b*200.0/my_parts[0]->get_mass());

	// average thermal velocity of test particle at 200K
	//double v_avg = sqrt(8.0*constants::k_b*200.0/(constants::pi*my_parts[0]->get_mass()));

	// background O velocity as defined in Justin's original code
	//double v_Obg = sqrt(8.0*constants::k_b*277.6 / (constants::pi*15.9994*constants::amu));

	upper_bound = my_planet.get_radius() + upper_bound;
	lower_bound = my_planet.get_radius() + lower_bound;
	double v_esc_upper = sqrt(2.0 * constants::G * my_planet.get_mass() / upper_bound);
	double global_rate = my_dist->get_global_rate();
	double k = my_planet.get_k_g();

	vector<int> active_indices;  // list of indices for active particles
	active_indices.resize(num_parts);
	for (int i=0; i<num_parts; i++)
	{
		active_indices[i] = i;
	}

	cout << "Simulating Particle Transport...\n";

	for (int i=0; i<num_steps; i++)
	{
		if (active_parts == 0)
		{
			break;
		}

		if (print_status_freq > 0 && (i+1) % print_status_freq == 0)
		{
			double hrs = (i+1)*dt/3600.0;
			double min = (hrs - (int)hrs)*60.0;
			double sec = (min - (int)min)*60.0;
			cout << (int)hrs << "h "<< (int)min << "m " << sec << "s " << "\tActive: " << active_parts << "\tDay escape: " << day_escape_count << "\tDay fraction: " << (double)day_escape_count / (double)(num_parts) << "\tNight escape: " << night_escape_count <<  "\tNight fraction: " << (double)night_escape_count / (double)(num_parts) <<endl;
		}

		if (output_pos_freq > 0 && (i+1) % output_pos_freq == 0)
		{
			output_positions(output_pos_dir + "positions" + to_string(i+1) + ".out");
		}

		if (num_traced > 0)
		{
			output_trace_data();
		}

		for (int j=0; j<active_parts; j++)
		{
			update_stats(dt, active_indices[j]);
			my_parts[active_indices[j]]->do_timestep(dt, k);

			if (bg_species.check_collision(my_parts[active_indices[j]], dt))
			{
				my_parts[active_indices[j]]->do_collision(bg_species.get_collision_target(), bg_species.get_collision_theta(), i*dt, my_planet.get_radius());
			}

			// escape velocity at current radius
			v_esc_current = sqrt(2.0 * constants::G * my_planet.get_mass() / my_parts[active_indices[j]]->get_radius());

			// thermalized threshold velocity; set to either v_esc_current, v_mp, v_rms, or v_avg
			// v_esc_current defined just above; others defined at beginning of run_simulation function
			v_thermal = v_esc_current;

			// deactivation criteria from Justin's original Hot O simulation code (must also uncomment v_Obg declaration above to use)
			//if (my_parts[active_indices[j]]->get_radius() < (my_planet.get_radius() + 900e5) && (my_parts[active_indices[j]]->get_total_v() + v_Obg) < sqrt(2.0*constants::G*my_planet.get_mass()*(my_parts[active_indices[j]]->get_inverse_radius()-1.0/(my_planet.get_radius()+900e5))))
			//{
			//	my_parts[active_indices[j]]->deactivate("\t\tParticle was thermalized.\n\n");
			//	active_parts--;
			//	active_indices.erase(active_indices.begin() + j);
			//	j--;
			//}

			if (my_parts[active_indices[j]]->get_total_v() < v_thermal)
			{
				my_parts[active_indices[j]]->deactivate(to_string(i*dt) + "\t\tParticle was thermalized.\n\n");
				active_parts--;
				active_indices.erase(active_indices.begin() + j);
				j--;
			}
			else if (my_parts[active_indices[j]]->get_radius() >= upper_bound && my_parts[active_indices[j]]->get_total_v() >= v_esc_upper)
			{
				if (my_parts[active_indices[j]]->get_x() > 0.0)
				{
					my_parts[active_indices[j]]->deactivate(to_string(i*dt) + "\t\tReached upper bound on day side with at least escape velocity.\n\n");
					active_parts--;
					day_escape_count++;
				}
				else
				{
					my_parts[active_indices[j]]->deactivate(to_string(i*dt) + "\t\tReached upper bound on night side with at least escape velocity.\n\n");
					active_parts--;
					night_escape_count++;
				}

				active_indices.erase(active_indices.begin() + j);
				j--;
			}
			else if (my_parts[active_indices[j]]->get_radius() <= lower_bound)
			{
				my_parts[active_indices[j]]->deactivate(to_string(i*dt) + "\t\tDropped below lower bound.\n\n");
				active_parts--;
				active_indices.erase(active_indices.begin() + j);
				j--;
			}
		}
	}

	if (num_traced > 0)
	{
		output_collision_data();
	}

	output_stats(dt, (global_rate / 2.0), num_parts, output_stats_dir);

	cout << "Number of collisions: " << bg_species.get_num_collisions() << endl;
	cout << "Active particles remaining: " << active_parts << endl;
	cout << "Number of day side escaped particles: " << day_escape_count << endl;
	cout << "Number of night side escaped particles: " << night_escape_count << endl;
	cout << "Total particles spawned: " << num_parts << endl;
	cout << "Day side fraction of escaped particles: " << (double)day_escape_count / (double)(num_parts) << endl;
	cout << "Night side fraction of escaped particles: " << (double)night_escape_count / (double)(num_parts) << endl;
	cout << "Global production rate: " << global_rate << endl;
	cout << "Total loss rate: " << ((double)day_escape_count / (double)(num_parts) + (double)night_escape_count / (double)(num_parts)) * (global_rate / 2.0) << endl;
}

void Atmosphere::update_stats(double dt, int i)
{
	double x = my_parts[i]->get_x();
	double y = my_parts[i]->get_y();
	double z = my_parts[i]->get_z();
	int x_index = 0;
	//int y_index = 0;
	int z_index = 0;
	int r_xz_index = 0;
	//int r_xy_index = 0;
	int r_3d_index = 0;
	double e = 0.0;
	int e_index = 0;
	//double inverse_v_r = 0.0;
	double cos_theta = 0.0;
	int cos_index = 0;

	//inverse_v_r = abs(dt / (my_parts[i]->get_radius() - my_parts[i]->get_previous_radius()));
	r_3d_index = (int)(1e-5*(my_parts[i]->get_radius()-my_planet.get_radius()));

	if (r_3d_index >= 0 && r_3d_index <= 100000)
	{
		if (x > 0.0)  // increment dayside density count
		{
			stats_dens_counts[0][r_3d_index] += 1;
		}
		else  // increment nightside density count
		{
			stats_dens_counts[1][r_3d_index] += 1;
		}
	}

	// update dayside integrated column density count for current altitude
	r_xz_index = (int)(1e-5*(sqrt(x*x + z*z) - my_planet.get_radius()));
	if ((x >= 0.0) && (r_xz_index >= 0) && (r_xz_index <= 100000)) //&& (abs(my_parts[i]->get_y()) <= 500e5))
	{
		stats_coldens_counts[r_xz_index] += 1;
	}

	x_index = (int)(1e-5*x/100.0);
	z_index = (int)(1e-5*z/100.0);
	if ((abs(x_index) <= 512) && ((abs(z_index) <= 512)))
	{
		x_index = x_index + 512;
		stats_dens2d_counts[z_index + 512][x_index] = stats_dens2d_counts[z_index + 512][x_index] + 1;
	}

	for (int j=0; j<stats_num_EDFs; j++)
	{
		if (r_3d_index == stats_EDF_alts[j])
		{
			e = my_parts[i]->get_energy_in_eV();
			e_index = (int)(20.0*e);

			cos_theta = my_parts[i]->get_cos_theta();
			cos_index = (int)(100.0*abs(cos_theta));

			if (cos_theta > 0.0)
			{
				cos_index += 100;
			}
			else
			{
				cos_index = abs(cos_index - 100);
			}

			double radial_v = abs((my_parts[i]->get_radius() - my_parts[i]->get_previous_radius()) / dt);
			if ((e_index >= 0 && e_index <= 200) && (cos_index >= 0 && cos_index <= 200))
			{
				if (x > 0.0)
				{
					stats_EDFs[0][j][e_index][cos_index] += 1;
				}
				else
				{
					stats_EDFs[1][j][e_index][cos_index] += 1;
				}
			}
			stats_loss_rates[j] = stats_loss_rates[j] + radial_v;
		}
	}
	
	// for bg angle-averaged density calculation in constant x plane (limb observation)
	for (int j=0; j<stats_num_EDFs; j++)
	  {
	    if ((int)(1e-5*(x-my_planet.get_radius())) == stats_EDF_alts[j]) // if particle in the slab where x = the chosen altitude (towards the Sun)
	      {
		if (sqrt(pow(y,2) + pow(z,2)) <= 1e5/2)
	      {
		stats_angleavg_dens[j] += 1;
		  }
		else if (sqrt(pow(y,2) + pow(z,2)) > 1e5/2)
		{
		  stats_angleavg_dens[j] += (2/constants::pi)*asin(1e5/(2*sqrt(pow(y,2) + pow(z,2))));
		}
		}
	  }
}

void Atmosphere::output_stats(double dt, double rate, int total_parts, string output_dir)
{
	double volume = 0.0;
	double surface_upper = 0.0;
	double r_in_cm = 0.0;
	double dens_day = 0.0;
	double dens_night = 0.0;
	double coldens_area = 0.0;
	double coldens_day = 0.0;
	double sum_day = 0.0;
	double sum_night = 0.0;
	ofstream dens_day_out, dens_night_out, coldens_day_out, dens2d_out, EDF_day_out, EDF_night_out, loss_rates_out, angleavg_dens_out;
	dens_day_out.open(output_dir + "density1d_day.out");
	dens_night_out.open(output_dir + "density1d_night.out");
	coldens_day_out.open(output_dir + "column_density_day.out");
	dens2d_out.open(output_dir + "density2d.out");
	angleavg_dens_out.open(output_dir + "angleavg_dens.out");

	dens_day_out << "#alt[km]\tdensity[cm-3]\n";
	dens_night_out << "#alt[km]\tdensity[cm-3]\n";
	coldens_day_out << "#alt[km]\tcol density[cm-2]\n";
	dens2d_out << "#this file contains a 1025 x 1025 grid of 2d integrated column densities for an observer viewing XZ plane from Y=infinity; each pixel represents a 100 square km area; density units are in particles per cm^2\n";

	int size = stats_dens_counts[0].size();
	for (int i=0; i<size; i++)
	{
		sum_day = 0.0;
		sum_night = 0.0;
		r_in_cm = my_planet.get_radius() + 1e5*(double)i;
		volume = 2.0*constants::pi/3.0 * (pow(r_in_cm+1e5, 3.0) - pow(r_in_cm, 3.0));

		sum_day = (double)stats_dens_counts[0][i];
		sum_night = (double)stats_dens_counts[1][i];

		dens_day = (dt*rate/(double)total_parts*sum_day) / volume;
	        dens_night = (dt*rate/(double)total_parts*sum_night) / volume;

		dens_day_out << i << "\t\t" << dens_day << "\n";
		dens_night_out << i << "\t\t" << dens_night << "\n";
	}
	dens_day_out.close();
	dens_night_out.close();

	// output altitude profile of integrated column densities
	for (int i=0; i<100001; i++)
	{
		r_in_cm = my_planet.get_radius() + 1e5*(double)i;
		coldens_area = 0.5 * constants::pi * (pow(r_in_cm+1e5, 2.0) - pow(r_in_cm, 2.0));
		coldens_day = (dt*rate/(double)total_parts*(double)stats_coldens_counts[i]) / coldens_area;
		coldens_day_out << i << "\t\t" << coldens_day << "\n";
	}
	coldens_day_out.close();

	// output 2d column density image data
	for (int i=0; i<1025; i++)
	{
		for (int j=0; j<1025; j++)
		{
			dens2d_out << (dt*rate/(double)total_parts*(double)stats_dens2d_counts[i][j]) / 1.0e14 << "\t";
		}
		dens2d_out << "\n";
	}
	dens2d_out.close();

	loss_rates_out.open(output_dir + "loss_rates.out");
	loss_rates_out << "#alt[km]\tloss rate[s-1]\n";
	for (int i=0; i<stats_num_EDFs; i++)
	{
		EDF_day_out.open(output_dir + "EDF_day_" + to_string(stats_EDF_alts[i]) + "km.out");
		EDF_night_out.open(output_dir + "EDF_night_" + to_string(stats_EDF_alts[i]) + "km.out");
		EDF_day_out << "# rows are energy, 0eV at top, 10eV at bottom; columns are cos(theta), -1 at left, 1 at right\n";
		EDF_night_out << "# rows are energy, 0eV at top, 10eV at bottom; columns are cos(theta), -1 at left, 1 at right\n";

		r_in_cm = 1e5*(double)stats_EDF_alts[i] + my_planet.get_radius();
		surface_upper = 2.0*constants::pi * (r_in_cm+1e5) * (r_in_cm+1e5);
		volume = 2.0*constants::pi/3.0 * (pow(r_in_cm+1e5, 3.0) - pow(r_in_cm, 3.0));

		angleavg_dens_out << stats_EDF_alts[i]  << "\t" << ( dt * rate * stats_angleavg_dens[i] ) /
		                                                   ((double)total_parts*1e5*1e5) << "\n";
				
		loss_rates_out << stats_EDF_alts[i] << "\t\t" << (stats_loss_rates[i] / volume) * (dt*rate/(double)total_parts) * surface_upper << "\n";

		for (int j=0; j<201; j++)
		{
			for (int k=0; k<201; k++)
			{
				EDF_day_out << ((dt*rate/(double)total_parts)*stats_EDFs[0][i][j][k]) / (volume*0.05*0.01) << "\t";
				EDF_night_out << ((dt*rate/(double)total_parts)*stats_EDFs[1][i][j][k]) / (volume*0.05*0.01) << "\t";
			}
			EDF_day_out << "\n";
			EDF_night_out << "\n";
		}
		EDF_day_out.close();
		EDF_night_out.close();
	}
	loss_rates_out.close();
	angleavg_dens_out.close();
}
