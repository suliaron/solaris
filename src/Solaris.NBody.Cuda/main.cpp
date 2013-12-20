#include <algorithm>
#include <cstdlib>
#include <time.h>
#include <iostream>
#include <fstream>

#include "ode.h"
#include "nbody.h"
#include "options.h"
#include "integrator.h"
#include "rungekutta.h"
#include "nbody_exception.h"

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

string get_printout_file(options& opt, int pcount)
{
	char buffer[1024];
	sprintf(buffer, "pos-%.5d.dat", pcount);
	return combine_path(opt.printoutDir, string(buffer));
}

string get_collisions_file(options& opt)
{
	return combine_path(opt.printoutDir, "col.dat");
}

int main(int argc, const char** argv)
{
	time_t start = time(NULL);

	try
	{
		options opt(argc, argv);

		nbody* nb			= (nbody*)opt.create_ode();
		integrator* intgr	= opt.create_integrator(nb);

		ttt_t pp			= 0;
		ttt_t ps			= 0;
		ttt_t dt			= 0;

		ostream* positionsf = 0;
		ostream* collisionsf= 0;
		int pcount			= 0;
		int ccount			= 0;

		if (!opt.printoutToFile) {
			positionsf = &cout;
			collisionsf = &cerr;
		}
		else {
			collisionsf = new ofstream(combine_path(opt.printoutDir, "col.dat").c_str());
		}

		while (nb->t < opt.timeStop) {
			if (opt.printout) {
				if (pp >= opt.printoutPeriod) {
					pp = 0;
				}

				// Start of a print-out period, create new file if necessary
				if (pp == 0) {
					cerr << (int)(nb->t / opt.timeStop * 100) << '%' << endl;
					if (opt.printoutToFile)	{
						if (positionsf) {
							delete positionsf;
						}
						positionsf = new ofstream(get_printout_file(opt, pcount++).c_str());
					}
					ccount = nb->print_collisions(*collisionsf, ccount);
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
			nb->detect_collisions();
		}
		ccount = nb->print_collisions(*collisionsf, ccount);

		if (opt.printoutToFile) {
			if (positionsf) {
				positionsf->flush();
				delete positionsf;
			}

			if (collisionsf) {
				collisionsf->flush();
				delete collisionsf;
			}
		}

		delete nb;
		delete intgr;
	} // end try bracket
	catch (nbody_exception& ex) {
		cerr << ex.what() << endl;
		options::print_usage();
	}
	time_t end = time(NULL);

	cerr << "Total time: " << end - start << " s" << endl;

	return 0;
}
