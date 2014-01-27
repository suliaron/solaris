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
#include "options.h"

using namespace std;

string combine_path(string dir, string path)
{
	if (dir.size() > 0) {
		if (*(dir.end() - 1) != '/' && *(dir.end() - 1) != '\\') {
			return dir + '/' + path;
		}
		else {
			return dir + path;
		}
	}
	else {
		return path;
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

void split_filename(const std::string& path)
{
	std::cout << "Splitting: " << path << '\n';
	unsigned found = path.find_last_of("/\\");
	std::cout << " dir:  " << path.substr(0,found) << " get_directory(): " << get_directory(path) << '\n';
	std::cout << " file: " << path.substr(found+1) << " get_filename(): " << get_filename(path) << '\n';
	std::cout << " file: " << path.substr(found+1) << " get_filename_without_ext(): " << get_filename_without_ext(path) << '\n';
	std::cout << " ext:  " << get_extension(path)  << " get_extension(): " << get_extension(path) << '\n';
}

string get_printout_file(options& opt, int pcount)
{
	char buffer[1024];
	sprintf(buffer, "pos-%.5d.txt", pcount);
	return combine_path(opt.printoutDir, string(buffer));
}

// -nBodies 1 1 0 0 0 0 19998 -i RK -dt 200 -t 3652500 -p 36525 36525 36525 -o C:\Work\Cuda.Solaris.TestRuns\2MStar_5MJupiter_Disc65-270_01 -f C:\Work\Cuda.Solaris.TestRuns\2MStar_5MJupiter_Disc65-270_01\Rezso.txt
// -nBodies 1 1 0 0 0 0 0 -i RKN -dt 200 -t 3652500 -p 36525 36525 36525 -o D:\Work\Projects\solaris\src\Solaris.NBody.Cuda.Test\TestRuns\2_Body -f D:\Work\Projects\solaris\src\Solaris.NBody.Cuda.Test\TestRuns\2_Body\TwoBody.txt
// -nBodies 1 1 0 0 0 0 0 -i RK -dt 200 -t 3652500 -p 36525 36525 36525 -o C:\Work\Cuda.Solaris.TestRuns\2_Body -f C:\Work\Cuda.Solaris.TestRuns\2_Body\TwoBody.txt
// -nBodies 1 0 0 0 1 0 0 -i RK -dt 40 -t 365250 -p 3652.5 3652.5 3652.5 -o C:\Work\Cuda.Solaris.TestRuns\2_Body_Gas -f C:\Work\Cuda.Solaris.TestRuns\2_Body_Gas\TwoBody.txt

int main(int argc, const char** argv)
{
	cout << "Solaris.NBody.Cuda.Test main.cu started" << endl;

	time_t start = time(NULL);

	// Integrate the nbody ode
	//try {
	//	options opt(argc, argv);

	//	nbody* nb			= opt.create_nbody();
	//	integrator* intgr	= opt.create_integrator(nb);

	//	ttt_t pp			= 0;
	//	ttt_t ps			= 0;
	//	ttt_t dt			= 0;

	//	ostream* positionsf = 0;
	//	ostream* collisionsf= 0;
		int pcount			= 0;
	//	int ccount			= 0;

	//	if (!opt.printoutToFile) {
	//		positionsf = &cout;
	//		collisionsf = &cerr;
	//	}
	//	else {
	//		//collisionsf = new ofstream(combine_path(opt.printoutDir, "col.txt").c_str());
	//		positionsf = new ofstream(get_printout_file(opt, pcount++).c_str());
	//	}

	//	while (nb->t < opt.timeStop) {

	//		if (opt.printout) {
	//			if (pp >= opt.printoutPeriod) {
	//				pp = 0;
	//			}

	//			// Start of a print-out period, create new file if necessary
	//			if (pp == 0) {
	//				cerr << setprecision(5) << setw(6) << (nb->t / opt.timeStop * 100) << " %" << endl;
	//			}

	//			if (0 <= pp && pp <= opt.printoutLength) {
	//				if (ps >= opt.printoutStep) {
	//					ps = 0;
	//				}

	//				if (ps == 0) {
	//					// Print out positions
	//					nb->copy_to_host();
	//					nb->print_positions(*positionsf);
	//				}
	//			}
	//		}
	//		dt = intgr->step();

	//		pp += dt;
	//		ps += dt;
	//	}

	//	delete nb;
	//	delete intgr;
	//}	// end try bracket
	//catch (nbody_exception& ex) {
	//	cerr << "Error: " << ex.what() << endl;
	//}
	//cout << "Total time: " << time(NULL) - start << " s" << endl;


	start = time(NULL);
	// Integrate the planets ode
	try {
		options opt(argc, argv);

		planets* pl			= opt.create_planets();
		integrator* intgr	= opt.create_integrator(pl);

		ttt_t pp			= 0;
		ttt_t ps			= 0;
		ttt_t dt			= 0;

		ostream* positionsf = 0;
		ostream* collisionsf= 0;
		int pcount			= 1;
		int ccount			= 0;

		string filename;
		if (!opt.printoutToFile) {
			positionsf = &cout;
			collisionsf = &cerr;
		}
		else {
			//collisionsf = new ofstream(combine_path(opt.printoutDir, "col.txt").c_str());
			//filename = get_filename_without_ext(opt.filename) + ".nBodies.";
			//int i = 1;
			//while (i < argc) {
			//	string p = argv[i];
			//	if (p == "-nBodies") {
			//		i++;
			//		string number;
			//		for (int k = i; k < i + 6; k++ ) {
			//			number = argv[k];
			//			filename += number + '_';
			//		}
			//		number = argv[i+6];
			//		filename += number;
			//		break;
			//	}
			//}
			//string filenameWithExt = filename + '.' + get_extension(opt.filename);
			//positionsf = new ofstream(get_printout_file(opt, pcount++).c_str());
			//positionsf = new ofstream(combine_path(opt.printoutDir, filename), std::ios::app);
		}

		while (pl->t < opt.timeStop) {

			if (opt.printout) {
				if (pp >= opt.printoutPeriod) {
					pp = 0;
				}

				// Start of a print-out period, create new file if necessary

				if (0 <= pp && pp <= opt.printoutLength) {
					if (ps >= opt.printoutStep) {
						ps = 0;
					}

					if (ps == 0) {
						// Print out positions
						pl->copy_to_host();

						if (positionsf) {
							delete positionsf;
						}
						positionsf = new ofstream(get_printout_file(opt, pcount++).c_str());
						pl->print_positions(*positionsf);
						positionsf->flush();
					}
				}
			}
			dt = intgr->step();
			cerr << setprecision(6) << setw(7) << dt << " [d]\t";
			cerr << setprecision(6) << setw(7) << (pl->t / opt.timeStop * 100) << " %" << endl;

			pp += dt;
			ps += dt;
		}

		delete pl;
		delete intgr;
	}	// end try bracket
	catch (nbody_exception& ex) {
		cerr << "Error: " << ex.what() << endl;
	}

	cout << "Total time: " << time(NULL) - start << " s" << endl;
	return 0;
}