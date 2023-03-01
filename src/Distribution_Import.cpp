/*
 * Distribution_Import.cpp
 *
 *  Created on: Jul 19, 2020
 *      Author: rodney
 */

#include "Distribution_Import.hpp"

Distribution_Import::Distribution_Import(Planet my_p, double ref_h, double ref_T, string pos_file, string vel_file)
	: Distribution(my_p, ref_h, ref_T) {
	num_particles = 0;
	next_index = 0;

	ifstream pos_input, vel_input;
	pos_input.open(pos_file);
	vel_input.open(vel_file);
	string pos_line, vel_line;

	while (getline(pos_input, pos_line) && getline(vel_input, vel_line))
	{
		num_particles++;
		stringstream str_pos(pos_line);
		stringstream str_vel(vel_line);

		positions.conservativeResize(positions.rows()+1, NoChange);
		velocities.conservativeResize(velocities.rows()+1, NoChange);

		str_pos >> positions(positions.rows()-1, 0) >> positions(positions.rows()-1, 1) >> positions(positions.rows()-1, 2);
		str_vel >> velocities(velocities.rows()-1, 0) >> velocities(velocities.rows()-1, 1) >> velocities(velocities.rows()-1, 2);
	}
	pos_input.close();
	vel_input.close();
}

Distribution_Import::~Distribution_Import() {

}

void Distribution_Import::init(shared_ptr<Particle> p)
{
	if (next_index < num_particles)
	{
		p->init_particle(positions(next_index, 0), positions(next_index, 1), positions(next_index, 2), velocities(next_index, 0), velocities(next_index, 1), velocities(next_index, 2));
		next_index++;
	}
	else
	{
		cout << "End of available particles in imported distribution met!" << endl;
	}
}

double Distribution_Import::get_global_rate()
{
	return 0.0;
}
