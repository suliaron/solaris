#pragma once

#include <cstdlib>

#include "config.h"
#include "euler.h"
#include "integrator.h"
#include "nbody.h"
#include "number_of_bodies.h"
#include "ode.h"
#include "planets.h"
#include "rungekutta.h"
#include "rungekuttanystrom.h"

using namespace std;

class options
{
public:
	typedef enum integrator_type
			{ 
				INTEGRATOR_EULER,
				INTEGRATOR_RUNGEKUTTA,
				INTEGRATOR_RUNGEKUTTANYSTROM
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

private:
	integrator_type_t inttype;		// Integrator type
	bool_t adaptive;				// Adaptive step size
	var_t tolerance;				// Tolerance
	bool_t file;					// Input file supplied
	int filen;						// Number of entries in input file
	bool_t random;					// Generate random data

	number_of_bodies*	bodies;

public:
	options(int argc, const char** argv);
	~options();


	static void print_usage();

	ode*		create_ode();
	nbody*		create_nbody();
	planets*	create_planets();
	integrator* create_integrator(ode* f);

private:
	void create_default_options();
	void parse_options(int argc, const char** argv);

	void initial_condition(nbody* nb);
	void initial_condition(planets* pl);
};