#pragma once

#include <cstdlib>

#include "config.h"
#include "integrator.h"
#include "nbody.h"
#include "ode.h"
#include "pp_disk.h"

class gas_disk;
class number_of_bodies;

using namespace std;

class options
{
public:
	typedef enum integrator_type
			{ 
				INTEGRATOR_EULER,
				INTEGRATOR_RUNGEKUTTA2,
				INTEGRATOR_OPT_RUNGEKUTTA2,
				INTEGRATOR_RUNGEKUTTA4,
				INTEGRATOR_OPT_RUNGEKUTTA4,
				INTEGRATOR_RUNGEKUTTANYSTROM,
				INTEGRATOR_OPT_RUNGEKUTTANYSTROM
			} integrator_type_t;


public:
	int		n;						// Number of bodies
	ttt_t	timeStart;				// Start time
	ttt_t	timeStop;				// Stop time
	ttt_t	dt;						// Initial time step
	var_t	buffer_radius;			// collision buffer
	bool_t	printout;				// Printout enabled
	bool_t	printoutToFile;			// Printout to file
	ttt_t	printoutPeriod;			// Printout period
	ttt_t	printoutStep;			// Printout step size	
	ttt_t	printoutLength;			// Printout length
	string	printoutDir;			// Printout directory
	string	filename;				// Input file name

	number_of_bodies	*nBodies;
	gas_disk			*gasDisk;

private:
	integrator_type_t inttype;		// Integrator type
	bool_t adaptive;				// Adaptive step size
	var_t tolerance;				// Tolerance
	bool_t file;					// Input file supplied
	int filen;						// Number of entries in input file
	bool_t random;					// Generate random data

public:
	options(int argc, const char** argv);
	~options();

	static void print_usage();

	ode*		create_ode();
	nbody*		create_nbody();
	pp_disk*	create_pp_disk();
	integrator* create_integrator(ode* f);

private:
	void create_default_options();
	void parse_options(int argc, const char** argv);

	void initial_condition(nbody* nb);
};
