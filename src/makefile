CFLAGS=-O2 #g -O0 -Wall -Wextra

corona3d_2020: Atmosphere.o Background_Species.o Common_Functions.o Distribution_Hot_H.o Distribution_Hot_O.o Distribution_Import.o Distribution_MB.o Distribution.o Interpolator.o main.o Particle_CO.o Particle_CO2.o Particle_H.o Particle_N2.o Particle_O.o Particle.o Planet.o
	g++ $(CFLAGS) Atmosphere.o Background_Species.o Common_Functions.o Distribution_Hot_H.o Distribution_Hot_O.o Distribution_Import.o Distribution_MB.o Distribution.o Interpolator.o main.o Particle_CO.o Particle_CO2.o Particle_H.o Particle_N2.o Particle_O.o Particle.o Planet.o -o corona3d_2020

Atmosphere.o: Atmosphere.cpp
	g++ $(CFLAGS) -c Atmosphere.cpp

Background_Species.o: Background_Species.cpp
	g++ $(CFLAGS) -c Background_Species.cpp

Common_Functions.o: Common_Functions.cpp
	g++ $(CFLAGS) -c Common_Functions.cpp

Distribution_Hot_H.o: Distribution_Hot_H.cpp
	g++ $(CFLAGS) -c Distribution_Hot_H.cpp

Distribution_Hot_O.o: Distribution_Hot_O.cpp
	g++ $(CFLAGS) -c Distribution_Hot_O.cpp

Distribution_Import.o: Distribution_Import.cpp
	g++ $(CFLAGS) -c Distribution_Import.cpp

Distribution_MB.o: Distribution_MB.cpp
	g++ $(CFLAGS) -c Distribution_MB.cpp

Distribution.o: Distribution.cpp
	g++ $(CFLAGS) -c Distribution.cpp

Interpolator.o: Interpolator.cpp
	g++ $(CFLAGS) -c Interpolator.cpp

main.o: main.cpp
	g++ $(CFLAGS) -c main.cpp

Particle_CO.o: Particle_CO.cpp
	g++ $(CFLAGS) -c Particle_CO.cpp

Particle_CO2.o: Particle_CO2.cpp
	g++ $(CFLAGS) -c Particle_CO2.cpp

Particle_H.o: Particle_H.cpp
	g++ $(CFLAGS) -c Particle_H.cpp

Particle_N2.o: Particle_N2.cpp
	g++ $(CFLAGS) -c Particle_N2.cpp

Particle_O.o: Particle_O.cpp
	g++ $(CFLAGS) -c Particle_O.cpp

Particle.o: Particle.cpp
	g++ $(CFLAGS) -c Particle.cpp

Planet.o: Planet.cpp
	g++ $(CFLAGS) -c Planet.cpp

clean:
	rm *.o
	rm corona3d_2020
