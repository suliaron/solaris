#pragma once

#include <cstdlib>

#include "config.h"
#include "integrator.h"
#include "euler.h"
#include "rungekutta.h"
#include "rungekuttanystrom.h"
#include "ode.h"
#include "nbody.h"

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
	int n;							// Number of bodies
	ttt_t timeStart;				// Start time
	ttt_t timeStop;					// Stop time
	ttt_t dt;						// Initial time step
	var_t buffer_radius;			// collision buffer
	bool printout;					// Printout enabled
	ttt_t printoutPeriod;			// Printout period
	ttt_t printoutStep;				// Printout step size	
	ttt_t printoutLength;			// Printout length
	bool printoutToFile;			// Printout to file
	string printoutDir;				// Printout directory

private:
	integrator_type_t inttype;		// Integrator type
	bool adaptive;					// Adaptive step size
	var_t tolerance;				// Tolerance
	bool file;						// Input file supplied
	int filen;						// Number of entries in input file
	string filename;				// Input file name
	bool random;					// Generate random data

public:
	options(int argc, const char** argv);
	~options();

	static void print_usage();

	ode* create_ode();
	integrator* create_integrator(ode * f);

private:
	void create_default_options();
	void parse_options(int argc, const char** argv);
};