// includes, system 
#include <ctime>
#include <iomanip>
#include <iostream>
#include <fstream>

// includes CUDA
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

// includes Thrust
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/reduce.h>

#include "config.h"
#include "nbody.h"
#include "nbody_exception.h"
#include "ode.h"
#include "pp_disk.h"
#include "options.h"

using namespace std;

string combine_path(string dir, string filename)
{
	if (dir.size() > 0) {
		if (*(dir.end() - 1) != '/' && *(dir.end() - 1) != '\\') {
			return dir + '/' + filename;
		}
		else {
			return dir + filename;
		}
	}
	else {
		return filename;
	}
}

string get_filename(const string& path)
{
	string result;

	if (path.size() > 0)
	{
		size_t pos = path.find_last_of("/\\");
		result = path.substr(pos + 1);
	}

	return result;
}

string get_filename_without_ext(const string& path)
{
	string result;

	if (path.size() > 0)
	{
		size_t pos = path.find_last_of("/\\");
		result = path.substr(pos + 1);
		pos = result.find_last_of('.');
		result = result.substr(0, pos);
	}

	return result;
}

string get_directory(const string& path)
{
	string result;

	if (path.size() > 0)
	{
		size_t pos = path.find_last_of("/\\");
		result = path.substr(0, pos);
	}

	return result;
}

string get_extension(const string& path)
{
	string result;

	if (path.size() > 0)
	{
		size_t pos = path.find_last_of('.');
		result = path.substr(pos + 1);
	}

	return result;
}

string get_printout_file(options& opt, int pcount)
{
	char buffer[1024];
	sprintf(buffer, "nBodies_1_127_0_0_0_0_0_RK_ppd2_pos-%.5d.txt", pcount);
	return combine_path(opt.printoutDir, string(buffer));
}

int main(int argc, const char** argv)
{
	cout << "Solaris.NBody.Cuda.Test main.cu started" << endl;

	time_t start = time(NULL);

	// Integrate the pp_disk ode
	try {
		options opt(argc, argv);

		pp_disk* ppd		= opt.create_pp_disk();
		integrator* intgr	= opt.create_integrator(ppd);

		ttt_t pp			= 0;
		ttt_t ps			= 0;
		ttt_t dt			= 0;

		ostream* positionsf = 0;
		ostream* orbelemf	= 0;
		//ostream* collisionsf= 0;
		int pcount			= 0;
		int ccount			= 0;

		string filename;
		if (!opt.printoutToFile) {
			positionsf = &cout;
			orbelemf   = &cout;
			//collisionsf = &cerr;
		}
		else {
			//collisionsf = new ofstream(combine_path(opt.printoutDir, "col.txt").c_str());
			//positionsf = new ofstream(get_printout_file(opt, pcount++).c_str());
			filename = get_filename_without_ext(opt.filename) + ".ppd.";
			string filenameWithExt = filename + get_extension(opt.filename);
			positionsf = new ofstream(combine_path(opt.printoutDir, filenameWithExt), std::ios::app);
			filenameWithExt = filename + "oe." + get_extension(opt.filename);
			orbelemf = new ofstream(combine_path(opt.printoutDir, filenameWithExt), std::ios::app);
		}

		while (ppd->t < opt.timeStop) {

			if (opt.printout) {
				if (pp >= opt.printoutPeriod) {
					pp = 0;
				}

				// Start of a print-out period, create new file if necessary
				if (pp == 0) {
					cerr << setprecision(10) << setw(16) << dt << " [d], ";
					cerr << setprecision(5) << setw(6) << (ppd->t / opt.timeStop * 100) << " %" << endl;
				}

				if (0 <= pp && pp <= opt.printoutLength) {
					if (ps >= opt.printoutStep) {
						ps = 0;
					}

					if (ps == 0) {
						// Print out positions
						ppd->copy_to_host();
						ppd->print_positions(*positionsf);
						ppd->calculate_orbelem(0);
						ppd->h_orbelem = ppd->d_orbelem;
						ppd->print_orbelem(*orbelemf);
					}
				}
			}
			dt = intgr->step();

			pp += dt;
			ps += dt;
		}

		delete ppd;
		delete intgr;
		delete positionsf;
		delete orbelemf;
	}	// end try bracket
	catch (nbody_exception& ex) {
		cerr << "Error: " << ex.what() << endl;
	}
	cout << "Total time: " << time(NULL) - start << " s" << endl;

	return 0;

	// Integrate the nbody ode
	try {
		options opt(argc, argv);

		nbody* nb			= opt.create_nbody();
		integrator* intgr	= opt.create_integrator(nb);

		ttt_t pp			= 0;
		ttt_t ps			= 0;
		ttt_t dt			= 0;

		ostream* positionsf = 0;
		//ostream* collisionsf= 0;
		int pcount			= 0;
		int ccount			= 0;

		if (!opt.printoutToFile) {
			positionsf = &cout;
			//collisionsf = &cerr;
		}
		else {
			//collisionsf = new ofstream(combine_path(opt.printoutDir, "col.txt").c_str());
			positionsf = new ofstream(get_printout_file(opt, pcount++).c_str());
		}

		while (nb->t < opt.timeStop) {

			if (opt.printout) {
				if (pp >= opt.printoutPeriod) {
					pp = 0;
				}

				// Start of a print-out period, create new file if necessary
				if (pp == 0) {
					cerr << setprecision(5) << setw(6) << (nb->t / opt.timeStop * 100) << " %" << endl;
				}

				if (0 <= pp && pp <= opt.printoutLength) {
					if (ps >= opt.printoutStep) {
						ps = 0;
					}

					if (ps == 0) {
						// Print out positions
						nb->copy_to_host();
						nb->print_positions(*positionsf);
					}
				}
			}
			dt = intgr->step();

			pp += dt;
			ps += dt;
		}

		delete nb;
		delete intgr;
		delete positionsf;
	}	// end try bracket
	catch (nbody_exception& ex) {
		cerr << "Error: " << ex.what() << endl;
	}
	cout << "Total time: " << time(NULL) - start << " s" << endl;

	return 0;
}