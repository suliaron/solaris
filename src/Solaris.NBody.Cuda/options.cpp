#include "options.h"
#include "nbody_exception.h"

options::options(int argc, const char** argv)
{
	create_default_options();
	parse_options(argc, argv);
}

options::~options()
{
}

void options::create_default_options()
{
	n = 256;
	inttype = INTEGRATOR_EULER;
	adaptive = false;
	tolerance = 0;
	timeStart = 0;
	timeStop = 1000;
	dt = (var_t)0.1;
	buffer_radius = 3;
	printout = false;
	printoutPeriod = 100;
	printoutStep = 1;
	printoutLength = 10;
	printoutToFile = false;
	printoutDir = "";
	file = false;
	filen = 0;
	filename = "";
	random = true;
}

void options::print_usage()
{
	cout << "Usage: CudaNBody <parameterlis>" << endl;
	cout << "Parameters:" << endl;
	cout << "     -n <number>   : Number of bodies" << endl;
	cout << "     -i <type>     : Integrator type" << endl;
	cout << "                     E : Euler" << endl;
	cout << "                     RK : Runge-Kutta" << endl;
	cout << "                     RKN : Runge-Kutta-Nystrom" << endl;
	cout << "     -a            : Adaptive time step" << endl;
	cout << "     -t <number>   : Tolerance" << endl;
	cout << "     -dt <number>  : Initial time step" << endl;
	cout << "     -b <number>   : Buffer factor for collisions" << endl;
	cout << "     -i <type>     : Integrator type" << endl;
	cout << "     -p <period> <length> <stepsize>" << endl;
	cout << "                   : Print-out enabled with given parameters" << endl;
	cout << "     -f <filename> <number> : Input file, number of entries" << endl;
	cout << "     -r            : Generate random data" << endl;
}

void options::parse_options(int argc, const char** argv)
{
	int i = 1;

	while (i < argc)
	{
		string p = argv[i];

		// Number of bodies
		if (p == "-n")
		{
			i++;
			n = atoi(argv[i]);
		}
		// Integrator type
		else if (p == "-i")
		{
			i++;
			p = argv[i];
			if (p == "E")
			{
				inttype = INTEGRATOR_EULER;
			}
			else if (p == "RK")
			{
				inttype = INTEGRATOR_RUNGEKUTTA;
			}
			else if (p == "RKN")
			{
				inttype = INTEGRATOR_RUNGEKUTTANYSTROM;
			}
			else
			{
				throw nbody_exception("Invalid integrator type.");
			}
		}
		// Adaptive method
		else if (p == "-a")
		{
			adaptive = true;
			i++;
			tolerance = (var_t)atof(argv[i]);
		}
		// Time end
		else if (p == "-t")
		{
			i++;
			timeStop = (var_t)atof(argv[i]);
		}
		// Time step
		else if (p == "-dt")
		{
			i++;
			dt = (var_t)atof(argv[i]);
		}
		// Radius buffer factor
		else if (p == "-b")
		{
			i++;
			buffer_radius = (var_t)atof(argv[i]);
		}
		// Print out period
		else if (p == "-p")
		{
			printout = true;
			i++;
			printoutPeriod = (var_t)atof(argv[i]);
			i++;
			printoutLength = (var_t)atof(argv[i]);
			i++;
			printoutStep = (var_t)atof(argv[i]);
		}
		// Print-out location
		else if (p == "-o")
		{
			i++;
			printoutDir = argv[i];
			printoutToFile = true;
		}
		// Input file
		else if (p == "-f")
		{
			i++;
			filename = argv[i];
			file = true;
		}
		else if (p == "-r")
		{
			random = true;
		}
		else
		{
			throw nbody_exception("Invalid switch on command-line.");
		}

		i++;
	}
}

ode* options::create_ode()
{
	nbody* nb = new nbody(n);

	nb->t = timeStart;
	
	if (file)
	{
		nb->load(filename, n);
	}
	/*else if (random)
	{
		nb->generate_random(n);
	}*/

	nb->copy_to_device();
	return nb;
}

integrator* options::create_integrator(ode * f)
{
	integrator* intgr;

	switch (inttype)
	{
	case INTEGRATOR_EULER:
		intgr = new euler(*f, dt);
		break;
	case INTEGRATOR_RUNGEKUTTA:
		intgr = new rungekutta<4>(*f, dt, adaptive, tolerance);
		break;
	case INTEGRATOR_RUNGEKUTTANYSTROM:
		intgr = new rungekuttanystrom<9>(*f, dt, adaptive, tolerance);
		break;
	default:
		throw nbody_exception("Not implemented.");
	}

	//intgr->set_adaptive(adaptive);

	return intgr;
}

