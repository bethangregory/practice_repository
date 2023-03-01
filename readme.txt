corona3d_2020 aims to be a flexible planetary exosphere modeling program, able to handle
multiple particle types with a variety of initial distribution options and background
species profiles. See the sample 'corona3d_2020.cfg' file under the 'src' folder for
details on the currently available configuration parameter options.

**********************
** QUICKSTART GUIDE **
**********************

This software has been tested on both Linux and Mac systems using the GCC compilier.
To get started, first download the source code to the desired directory from a terminal:

     cd <desired directory>
     git clone https://github.com/rdelliott76/corona3d_2020.git

Now cd into the src directory and compile the application:

     cd corona3d_2020/src/
     make

Note that you may receive a compilation error due to missing the Eigen library. If this
happens, try this on Ubuntu:

     sudo apt install libeigen3-dev

On Mac:

     brew install eigen

You may need to look up how to install the Eigen library on your particular system if
neither of these commands works. Once installed, retry make, and it should work now. The relevant
'include' line may also need to be commented in in Particle.hpp.

Once installed, modify the corona3d_2020.cfg file to have the desired parameters. To
run the program, execute the following command from within the src directory:

     ./corona3d_2020



*************************
** SOURCE CODE DETAILS **
*************************

Following is a (somewhat) detailed explanation of the source code files and program logic
to aid future developers who would like to contribute to this software. There are also
fairly extensive comments throughout all of the source code files to aid in understanding.
A general knowledge of object-oriented programming concepts in C++ will be immensely helpful.

** GENERAL PROGRAM LOGIC OVERVIEW **

The first noteable bit of code that gets executed upon running the program is in the
Common_Functions.cpp file, specifically the get_seed() function. This bit of code seeds the
random number generator. It first looks for a file named "rng_seed" within the same directory
as the executable. This file should contain a single long long integer if you would like to run
a replicateable simulation for testing purposes. If the "rng_seed" file does not exist, then
the current system clock time is used as the seed. Execution now passes to main.cpp.

The main.cpp source file begins by reading in the main configuration file. This file must be
named "corona3d_2020.cfg" and must be in the same directory as the executeable. After reading
in all of the needed paramerters, main.cpp then constructs all of the necessary objects to
run the simulation. These objects are constructed and initialized in the following order: Planet,
array of Particles (actually, an array of shared pointers to Particle child classes, usually O or H),
Distribution (actually, a shared pointer to one of the Distribution child classes, usually
Distribution_Hot_H or Distribution_Hot_O), Background_Species, Atmosphere.

Once all the necessary objects are constructed, the run_simulation method within Atmosphere.cpp is called.
This routine executes the main loop of the program. More specifically, it consists of an outer for loop
that iterates for the number of time steps specified in the main configuration file. Within this loop,
there is a nested for loop that cycles through all of the active particles in the simulation. For each time
step, for each active particle, several things happen:  First, any console output or particle trace data is
output, if necessary. Then, the radial density and EDF counts are updated. Next, the particle's position is
updated according to its current velocity vector by timestepping the 3-dimensional equation of motion using
the value of dt supplied in the main configuration file. Then, the particle's velocity vector is updated using
the acceleration due to gravity at its current altitude. Next, a collision probability check is performed using
the background species' densities at the current altitude along with the available energy-dependent collision
cross sections supplied in the configuration files. If it is determined that a collision occurred, then the
particle's velocity vector is again updated using convervation of momentum (only elastic collisions are considered
currently). Finally, a deactivation check is done using three deacitivation criteria.

