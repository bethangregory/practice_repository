/*
 * Particle.cpp
 *
 *  Created on: May 27, 2020
 *      Author: rodney
 */

#include "Particle.hpp"
#include "Particle_CO.hpp"
#include "Particle_CO2.hpp"
#include "Particle_H.hpp"
#include "Particle_N2.hpp"
#include "Particle_O.hpp"
#include <iostream>

Particle::Particle()
{
	active = true;
	traced = false;
	radius = 0.0;
	inverse_radius = 0.0;
	previous_radius = 0.0;
	position[0] = position[1] = position[2] = 0.0;
	velocity[0] = velocity[1] = velocity[2] = 0.0;
}

Particle::~Particle()
{

}

// deactivate this particle
void Particle::deactivate(string fate)
{
	active = false;

	// record fate of particle at bottom of collision log
	if (traced)
	{
		collision_log.push_back(fate);
	}
}

// perform collision on a particle and update velocity vector
void Particle::do_collision(shared_ptr<Particle> target, double theta, double time, double planet_r)
{
	double v_before, v_after;
	double my_mass = get_mass();
	double targ_mass = target->get_mass();
	Matrix<double, 3, 1> targ_v = {target->get_vx(), target->get_vy(), target->get_vz()};
	Matrix<double, 3, 1> vcm;
	Matrix<double, 3, 3> Rrg;

	vcm = (my_mass*velocity.array() + targ_mass*targ_v.array()) / (my_mass + targ_mass);
	Matrix<double, 3, 1> v1v = velocity.array() - vcm.array();        // particle 1 c-o-m velocity
	// Matrix<double, 3, 1> v2v = targ_v.array() - vcm.array();          // particle 2 c-o-m velocity
	double v1 = sqrt(v1v[0]*v1v[0] + v1v[1]*v1v[1] + v1v[2]*v1v[2]);  // particle 1 c-o-m scalar velocity
	// double v2 = sqrt(v2v[0]*v2v[0] + v2v[1]*v2v[1] + v2v[2]*v2v[2]);  // particle 2 c-o-m scalar velocity

	// unit vector parallel to particle 1 velocity
	Matrix<double, 3, 1> r = velocity.array() / sqrt(velocity[0]*velocity[0] + velocity[1]*velocity[1] + velocity[2]*velocity[2]); // velocity is a three-component vector: vx, vy, vz, so r is the velocity vector normalised to -1 to 1

	double alpha = atan2(velocity[1], velocity[0]);
	double phi = atan2(velocity[2], sqrt(velocity[0]*velocity[0] + velocity[1]*velocity[1]));
	double gamma = constants::twopi*common::get_rand();  //  gamma random no. from uniform distribution between 0 and 2pi

	Matrix<double, 3, 1> vp;
	vp[0] = v1*cos(alpha)*cos(phi-theta);
	vp[1] = v1*sin(alpha)*cos(phi-theta);
	vp[2] = v1*sin(phi-theta);

	double Cg = cos(gamma);
	double Sg = sin(gamma);
	double Vg = 1.0-Cg;

	Rrg(0, 0) = r[0]*r[0]*Vg+Cg;
	Rrg(0, 1) = r[0]*r[1]*Vg+r[2]*Sg;
	Rrg(0, 2) = r[0]*r[2]*Vg-r[1]*Sg;

	Rrg(1, 0) = r[0]*r[1]*Vg-r[2]*Sg;
	Rrg(1, 1) = r[1]*r[1]*Vg+Cg;
	Rrg(1, 2) = r[1]*r[2]*Vg+r[0]*Sg;

	Rrg(2, 0) = r[0]*r[2]*Vg+r[1]*Sg;
	Rrg(2, 1) = r[1]*r[2]*Vg-r[0]*Sg;
	Rrg(2, 2) = r[2]*r[2]*Vg+Cg;


	Matrix<double, 3, 1> vrel1 = Rrg * vp;

	// update post-collision velocity
	if (traced)
	{
	  v_before = get_total_v()*1e-5; //  sqrt(vx^2 + vy^2 + vz^2)
	}
	velocity = vcm.array() + vrel1.array(); // here is where velocity changes

	// in case you need the updated collision partner velocity for something
	// targ_v = vcm.array() - ((my_mass / targ_mass) * vrel1.array());

	// write to collision log if traced particle
	if (traced)
	{
		v_after = get_total_v()*1e-5;
		double alt_in_km = 1e-5*(radius - planet_r);
		collision_log.push_back(to_string(time) + "\t\t" + to_string(alt_in_km) + "\t" + target->get_name() + "\t" + to_string(theta * (180.0/constants::pi)) + "\t" + to_string(v_before) + "\t" + to_string(v_after));
	}
}

void Particle::do_timestep(double dt, double k_g)
{
	previous_radius = radius;  // record current radius as new previous radius
	Array<double, 3, 1> a = {0.0, 0.0, 0.0}; // particle acceleration vector

	// calculate acceleration at current position
	double inv_r_cube = inverse_radius*inverse_radius*inverse_radius;
	a = k_g*position.array()*inv_r_cube;

	// calculate next position and update particle
	position.array() = position.array() + (velocity.array()*dt) + (0.5*a*dt*dt);
	radius = sqrt(position[0]*position[0] + position[1]*position[1] + position[2]*position[2]);
	inverse_radius = 1.0 / radius;

	// calculate acceleration at next position
	inv_r_cube = inverse_radius*inverse_radius*inverse_radius;
	a = a + k_g*position.array()*inv_r_cube;

	// calculate next velocity using acceleration at next position and update particle
	velocity.array() = velocity.array() + 0.5*a*dt;
}

// write collision log to given file
void Particle::dump_collision_log(string filename)
{
	ofstream outfile;
	outfile.open(filename);
	outfile << "#time(s)" << "\t\t" << "alt(km)" << "\t" << "targ" << "\t" << "angle(deg)" << "\t" << "v_bef(km/s)" << "\t" << "v_aft(km/s)\n";
	int num_lines = collision_log.size();
	for (int i=0; i<num_lines; i++)
	{
		outfile << collision_log[i] << "\n";
	}
	outfile.close();
}

bool Particle::is_active() const
{
	return active;
}

bool Particle::is_traced() const
{
	return traced;
}

// return cosine of angle between particle trajectory and normal
double Particle::get_cos_theta() const
{
	double cos_theta = get_radial_v() / get_total_v();

	if (cos_theta > 1.0)
	{
		cos_theta = 1.0;
	}
	else if (cos_theta < -1.0)
	{
		cos_theta = -1.0;
	}

	return cos_theta;
}

double Particle::get_energy_in_eV() const
{
	return 0.5*get_mass()*pow(get_total_v(), 2.0)/constants::ergev;
}

double Particle::get_radial_energy_in_eV(double dt) const
{
	double radial_v = abs(radius - previous_radius) / dt;
	return 0.5*get_mass()*pow(radial_v, 2.0)/constants::ergev;
}

double Particle::get_radius() const
{
	return radius;
}

double Particle::get_radial_v() const
{
	double radial_v = (velocity[0]*position[0] + velocity[1]*position[1] + velocity[2]*position[2]) / radius;
	return radial_v;
}

double Particle::get_inverse_radius() const
{
	return inverse_radius;
}

double Particle::get_previous_radius() const
{
	return previous_radius;
}

double Particle::get_x() const
{
	return position[0];
}

double Particle::get_y() const
{
	return position[1];
}

double Particle::get_z() const
{
	return position[2];
}

double Particle::get_vx() const
{
	return velocity[0];
}

double Particle::get_vy() const
{
	return velocity[1];
}

double Particle::get_vz() const
{
	return velocity[2];
}

double Particle::get_total_v() const
{
	return sqrt(velocity[0]*velocity[0] +
			    velocity[1]*velocity[1] +
				velocity[2]*velocity[2]);
}

// initialize particle using given position and velocity
void Particle::init_particle(double x, double y, double z, double vx, double vy, double vz)
{
	radius = sqrt(x*x + y*y + z*z);
	inverse_radius = 1.0/radius;
	previous_radius = radius;
	position[0] = x;
	position[1] = y;
	position[2] = z;
	velocity[0] = vx;
	velocity[1] = vy;
	velocity[2] = vz;
}

// initialize particle with velocity only to be used for collisions
void Particle::init_particle_vonly(double vx, double vy, double vz)
{
	velocity[0] = vx;
	velocity[1] = vy;
	velocity[2] = vz;
}

// initialize a single particle at given radius using Maxwell-Boltzmann avg v
void Particle::init_particle_MB(double r, double v_avg)
{
	radius = r;
	double phi = constants::twopi*(common::get_rand());
	double u = 2.0*common::get_rand() - 1;
	inverse_radius = 1.0/r;
	previous_radius = radius;
	position[0] = r*sqrt(1-(u*u))*cos(phi);
	position[1] = r*sqrt(1-(u*u))*sin(phi);
	position[2] = r*u;

	double randnum1 = common::get_rand();
	double randnum2 = common::get_rand();
	double randnum3 = common::get_rand();
	double randnum4 = common::get_rand();

	randnum1 = v_avg*sqrt(-2.0*log(1.0-randnum1));
	randnum2 = constants::twopi*randnum2;
	randnum3 = v_avg*sqrt(-2.0*log(1.0-randnum3));
	randnum4 = constants::twopi*randnum4;

	velocity[0] = randnum1*cos(randnum2);
	velocity[1] = randnum1*sin(randnum2);
	velocity[2] = randnum3*cos(randnum4);
}

void Particle::init_particle_vonly_MB(double v_avg)
{
	double randnum1 = common::get_rand();
	double randnum2 = common::get_rand();
	double randnum3 = common::get_rand();
	double randnum4 = common::get_rand();

	randnum1 = v_avg*sqrt(-2.0*log(1.0-randnum1));
	randnum2 = constants::twopi*randnum2;
	randnum3 = v_avg*sqrt(-2.0*log(1.0-randnum3));
	randnum4 = constants::twopi*randnum4;

	velocity[0] = randnum1*cos(randnum2);
	velocity[1] = randnum1*sin(randnum2);
	velocity[2] = randnum3*cos(randnum4);
}

void Particle::set_traced()
{
	traced = true;
}
