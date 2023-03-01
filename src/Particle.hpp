/*
 * Particle.hpp
 *
 *  Created on: May 27, 2020
 *      Author: rodney
 */

#ifndef PARTICLE_HPP_
#define PARTICLE_HPP_

//#include <eigen3/Eigen/Core>  // uncomment for Ubuntu
//#include </usr/local/Cellar/eigen/3.3.9/include/eigen3/Eigen/Core>  // uncomment for Mac
#include </opt/local/include/eigen3/Eigen/Core> // uncomment for Mac option 2
#include "Common_Functions.hpp"
using namespace Eigen;

class Particle {
public:
	Particle();
	virtual ~Particle();
	virtual double get_mass() const = 0;  //must be implemented in derived classes
	virtual string get_name() const = 0;  //must be implemented in derived classes

	void deactivate(string fate);
	void do_collision(shared_ptr<Particle> target, double theta, double time, double planet_r);
	void do_timestep(double dt, double k_g);
	void dump_collision_log(string filename);
	bool is_active() const;
	bool is_thermalized() const;
	bool is_traced() const;
	double get_cos_theta() const;
	double get_energy_in_eV() const;
	double get_radial_energy_in_eV(double dt) const;
	double get_radius() const;
	double get_radial_v() const;
	double get_inverse_radius() const;
	double get_previous_radius() const;
	double get_x() const;
	double get_y() const;
	double get_z() const;
	double get_vx() const;
	double get_vy() const;
	double get_vz() const;
	double get_total_v() const;
	double get_weight() const;
	void init_particle(double x, double y, double z, double vx, double vy, double vz);
	void init_particle_vonly(double vx, double vy, double vz);
	void init_particle_MB(double r, double v_avg); // init particle from MB distribution
	void init_particle_vonly_MB(double v_avg);     // init with velocity only for collision partners
	void set_traced();

protected:
	bool active;                    // flag for whether particle is active, i.e. should still be considered in the simulation
	bool traced;                    // flag for whether particle is traced through simulation
	double radius;                  // radius from center of planet	[cm]
	double inverse_radius;          // inverse radius (for computational efficiency) [cm^-1]
	double previous_radius;         // radius at previous time step; used for tracking
	Matrix<double, 3, 1> position;  // position vector [cm,cm,cm]
	Matrix<double, 3, 1> velocity;  // velocity vector [cm/s,cm/s,cm/s]
	vector<string> collision_log;   // log of collision data kept on traced particles
};

#endif /* PARTICLE_HPP_ */
