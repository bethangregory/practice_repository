/*
 * Common_Functions.cpp
 *
 *  Created on: Aug 21, 2020
 *      Author: rodney
 */

#include "Common_Functions.hpp"

// function to check if custom random seed exists in local file "rng_seed"
// if file does not exist, uses system clock to generate seed
static long long get_seed()
{
	long long s = 0;
	ifstream seedfile;
	seedfile.open("rng_seed");
	if (!seedfile.good())
	{
		s = chrono::high_resolution_clock::now().time_since_epoch().count();
	}
	else
	{
		string line;
		getline(seedfile, line);
		s = stoll(line);
	}
	seedfile.close();
	return s;
}

// seed random number generator using get_seed() function above
// to access externally, must include "Common_Functions.hpp" and call
// using common::get_rand() (will return uniform real between 0 and 1)
static long long seed = get_seed();
static mt19937 rand_generator(seed);   // Mersenne Twister PRNG (apparently, pretty good)
static uniform_real_distribution<double> rand_dist(0.0, 1.0);  // dist to be used with get_rand()

namespace constants {
	const double pi    = M_PI;            // pi [unitless]
	const double twopi = 2*pi;             // 2*pi [unitless]
	const double k_b   = 1.380649e-16;     // Boltzmann's Constant [erg/K]
	const double c     = 29979245800.0;    // Speed of Light in Vacuum [cm/s]
	const double G     = 6.67430e-8;       // Gravitational Constant [cm^3/g/s^2]
	const double amu   = 1.660538782e-24;  // Atomic Mass Unit [g]
	const double m_e   = 9.10938215e-28;   // Electron Mass [g]
	const double q_e   = 1.602176487e-19;  // Elementary Charge [C]
	const double jev   = q_e;              // Joules/Electron Volt [unitless]
	const double ergev = jev*1.0e7;        // ergs/Electron Volt [unitless]
}

namespace common {

	void import_csv(string filename, vector<double> &col1, vector<double> &col2)
	{
		ifstream infile;
		infile.open(filename);
		if (!infile.good())
		{
			cout << "\"" << filename << "\" not found!\n";
			exit(1);
		}
		string line, word;
		vector<string> row;
		while (getline(infile, line))
		{
			row.clear();
			if (line[0] == '#' || line.empty() || std::all_of(line.begin(), line.end(), ::isspace))
			{
				continue;
			}
			else
			{
				stringstream str(line);
				while(getline(str, word, ','))
				{
					row.push_back(word);
				}
				col1.push_back(stod(row[0]));
				col2.push_back(stod(row[1]));
			}
		}
		infile.close();
	}

	void import_csv(string filename, vector<double> &col1, vector<double> &col2, vector<double> &col3)
	{
		ifstream infile;
		infile.open(filename);
		if (!infile.good())
		{
			cout << "\"" << filename << "\" not found!\n";
			exit(1);
		}
		string line, word;
		vector<string> row;
		while (getline(infile, line))
		{
			row.clear();
			if (line[0] == '#' || line.empty() || std::all_of(line.begin(), line.end(), ::isspace))
			{
				continue;
			}
			else
			{
				stringstream str(line);
				while(getline(str, word, ','))
				{
					row.push_back(word);
				}
				col1.push_back(stod(row[0]));
				col2.push_back(stod(row[1]));
				col3.push_back(stod(row[2]));
			}
		}
		infile.close();
	}

	void import_csv(string filename, vector<double> &col1, vector<double> &col2, vector<double> &col3, vector<double> &col4)
	{
		ifstream infile;
		infile.open(filename);
		if (!infile.good())
		{
			cout << "\"" << filename << "\" not found!\n";
			exit(1);
		}
		string line, word;
		vector<string> row;
		while (getline(infile, line))
		{
			row.clear();
			if (line[0] == '#' || line.empty() || std::all_of(line.begin(), line.end(), ::isspace))
			{
				continue;
			}
			else
			{
				stringstream str(line);
				while(getline(str, word, ','))
				{
					row.push_back(word);
				}
				col1.push_back(stod(row[0]));
				col2.push_back(stod(row[1]));
				col3.push_back(stod(row[2]));
				col4.push_back(stod(row[3]));
			}
		}
		infile.close();
	}

	void import_csv(string filename, vector<double> &col1, vector<double> &col2, vector<double> &col3, vector<double> &col4, vector<double> &col5)
	{
		ifstream infile;
		infile.open(filename);
		if (!infile.good())
		{
			cout << "\"" << filename << "\" not found!\n";
			exit(1);
		}
		string line, word;
		vector<string> row;
		while (getline(infile, line))
		{
			row.clear();
			if (line[0] == '#' || line.empty() || std::all_of(line.begin(), line.end(), ::isspace))
			{
				continue;
			}
			else
			{
				stringstream str(line);
				while(getline(str, word, ','))
				{
					row.push_back(word);
				}
				col1.push_back(stod(row[0]));
				col2.push_back(stod(row[1]));
				col3.push_back(stod(row[2]));
				col4.push_back(stod(row[3]));
				col5.push_back(stod(row[4]));
			}
		}
		infile.close();
	}

	void import_csv(string filename, vector<double> &col1, vector<double> &col2, vector<double> &col3, vector<double> &col4, vector<double> &col5, vector<double> &col6)
	{
		ifstream infile;
		infile.open(filename);
		if (!infile.good())
		{
			cout << "\"" << filename << "\" not found!\n";
			exit(1);
		}
		string line, word;
		vector<string> row;
		while (getline(infile, line))
		{
			row.clear();
			if (line[0] == '#' || line.empty() || std::all_of(line.begin(), line.end(), ::isspace))
			{
				continue;
			}
			else
			{
				stringstream str(line);
				while(getline(str, word, ','))
				{
					row.push_back(word);
				}
				col1.push_back(stod(row[0]));
				col2.push_back(stod(row[1]));
				col3.push_back(stod(row[2]));
				col4.push_back(stod(row[3]));
				col5.push_back(stod(row[4]));
				col6.push_back(stod(row[5]));
			}
		}
		infile.close();
	}

	// returns interpolated value at x from parallel arrays (x_data, y_data)
	// assumes that x_data has at least two elements, is sorted and is strictly monotonic increasing
	double interpolate(vector<double> &x_data, vector<double> &y_data, double x)
	{
		/*
		int size = x_data.size();

		int i = 0;                    // find left end of interval for interpolation
		if (x >= x_data[size - 2])    // special case: beyond right end
		{
			i = size - 2;
		}
		else
		{
			while (x > x_data[i+1]) i++;
		}
		double x_left = x_data[i], y_left = y_data[i], x_right = x_data[i+1], y_right = y_data[i+1];  // points on either side (unless beyond ends)
		*/

		if (x <= x_data[0])
		{
			return y_data[0];
		}
		else if (x >= x_data.back())
		{
			return y_data.back();
		}
		else
		{
			int i = lower_bound(x_data.begin(), x_data.end(), x) - x_data.begin();
			double x_left = x_data[i-1], x_right = x_data[i], y_left = y_data[i-1], y_right = y_data[i];
			double dydx = (y_right - y_left) / (x_right - x_left);   // gradient
			return y_left + dydx * (x - x_left);   // linear interpolation
		}
	}

	// returns interpolated value at x from arrays x_data, y_data, when y_data is log-scaled
	double interpolate_logy(vector<double> &x_data, vector<double> &y_data, double x)
	{
		if (x <= x_data[0])
		{
			return y_data[0];
		}
		else if (x >= x_data.back())
		{
			return y_data.back();
		}
		else
		{
			int i = lower_bound(x_data.begin(), x_data.end(), x) - x_data.begin();
			double x_left = x_data[i-1], x_right = x_data[i], y_left = log10(y_data[i-1]), y_right = log10(y_data[i]);
			double dydx = (y_right - y_left) / (x_right - x_left);   // gradient
			return pow(10.0, y_left + dydx * (x - x_left));   // linear interpolation
		}
	}

	// returns uniformly distributed random number from interval [0, 1)
	double get_rand()
	{
		return rand_dist(rand_generator);
	}

	// returns uniformly distributed random integer between lower and upper (inclusive)
	int get_rand_int(int lower, int upper)
	{
		uniform_int_distribution<int> dist(lower, upper);
		return dist(rand_generator);
	}
}
