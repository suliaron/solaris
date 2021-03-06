#include "Constants.h"
#include "gas_disk.h"
#include "options.h"
#include "nbody_exception.h"
#include "number_of_bodies.h"

#include "euler.h"
#include "midpoint.h"
#include "rk4.h"
#include "rkn76.h"
#include "rungekutta.h"
#include "rungekuttanystrom.h"

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
	n				= 256;
	nBodies			= 0;
	inttype			= INTEGRATOR_EULER;
	adaptive		= false;
	tolerance		= 0;
	timeStart		= 0;
	timeStop		= 1000;
	dt				= (var_t)0.1;
	buffer_radius	= 3;
	printout		= false;
	printoutPeriod	= 100;
	printoutStep	= 1;
	printoutLength	= 10;
	printoutToFile	= false;
	printoutDir		= "";
	file			= false;
	filen			= 0;
	filename		= "";
	random			= true;
	gasDisk			= 0;
}

void options::print_usage()
{
	cout << "Usage: CudaNBody <parameterlis>" << endl;
	cout << "Parameters:" << endl;
	cout << "     -n <number>        : Number of self-interacting bodies" << endl;
	cout << "     -nBodies <nStar> <nGP> <nRP> <nPP> <nSP> <nPl> <nTP> : Number of bodies" << endl;
	cout << "     -i <type>          : Integrator type" << endl;
	cout << "                          E : Euler" << endl;
	cout << "                          oRK2 : optimized 2nd order Runge-Kutta" << endl;
	cout << "                          oRK4 : optimized 4th order Runge-Kutta" << endl;
	cout << "                          oRKN : optimized 7(6) Runge-Kutta-Nystrom" << endl;
	cout << "                          RK2  : 2nd order Runge-Kutta" << endl;
	cout << "                          RK4  : 4th order Runge-Kutta" << endl;
	cout << "                          RKN  : Runge-Kutta-Nystrom" << endl;
	cout << "     -gas          : Embed the planets into a gas disk" << endl;
	cout << "     -a <number>   : Use adaptive time step with <number> as tolerance" << endl;
	cout << "     -t0 <number>  : Start time " << endl;
	cout << "     -t <number>   : Stop time " << endl;
	cout << "     -dt <number>  : Initial time step" << endl;
	cout << "     -b <number>   : Buffer factor for collisions" << endl;
	cout << "     -p <period> <length> <stepsize>" << endl;
	cout << "                   : Print-out enabled with given parameters" << endl;
	cout << "     -o <directory>: Output directory" << endl;
	cout << "     -f <filename> : Input file, number of entries" << endl;
	cout << "     -r            : Generate random data" << endl;
}

void options::parse_options(int argc, const char** argv)
{
	int i = 1;

	while (i < argc) {
		string p = argv[i];

		// Number of bodies
		if (p == "-n") {
			i++;
			n = atoi(argv[i]);
			if (2 > n) {
				throw nbody_exception("Number of bodies must exceed 2.");
			}
		}
		else if (p == "-nBodies") {
			i++;
			int	star				= atoi(argv[i++]);
			int	giant_planet		= atoi(argv[i++]);
			int	rocky_planet		= atoi(argv[i++]);
			int	proto_planet		= atoi(argv[i++]);
			int	super_planetesimal	= atoi(argv[i++]);
			int	planetesimal		= atoi(argv[i++]);
			int	test_particle		= atoi(argv[i]);
			this->nBodies = new number_of_bodies(star, giant_planet, rocky_planet, proto_planet, super_planetesimal, planetesimal, test_particle);
		}
		// Initialize a gas_disk object with default values
		else if (p == "-gas") {
			var2_t eta = {2.0e-3,   1.0/2.0	};
			var2_t rho = {1.0e-9, -11.0/4.0	};		// g / cm^3
			var2_t sch = {5.0e-2,   5.0/4.0	};
			var2_t tau = {2.0/3.0,  2.0		};
			rho.x	*= Constants::GramPerCm3ToSolarPerAu3; // M_sun / AU^3
			ttt_t	t0 = 0.0;
			ttt_t	t1 = 100.0 * Constants::YearToDay;
			ttt_t	timeScale = 10.0 * Constants::YearToDay;
			gas_disk::gas_decrease_t gasDecrease = gas_disk::gas_decrease_t::CONSTANT;
			gasDisk = new gas_disk(rho, sch, eta, tau, gasDecrease, t0, t1, timeScale);
		}
		// Integrator type
		else if (p == "-i") {
			i++;
			p = argv[i];
			if (p == "E") {
				inttype = INTEGRATOR_EULER;
			}
			else if (p == "RK2")	{
				inttype = INTEGRATOR_RUNGEKUTTA2;
			}
			else if (p == "oRK2")	{
				inttype = INTEGRATOR_OPT_RUNGEKUTTA2;
			}
			else if (p == "RK4")	{
				inttype = INTEGRATOR_RUNGEKUTTA4;
			}
			else if (p == "oRK4")	{
				inttype = INTEGRATOR_OPT_RUNGEKUTTA4;
			}
			else if (p == "RKN") {
				inttype = INTEGRATOR_RUNGEKUTTANYSTROM;
			}
			else if (p == "oRKN") {
				inttype = INTEGRATOR_OPT_RUNGEKUTTANYSTROM;
			}
			else {
				throw nbody_exception("Invalid integrator type.");
			}
		}
		// Adaptive method
		else if (p == "-a")	{
			adaptive = true;
			i++;
			tolerance = (var_t)atof(argv[i]);
		}
		// Time start
		else if (p == "-t0")	{
			i++;
			timeStart = (var_t)atof(argv[i]) * Constants::YearToDay;
		}
		// Time end
		else if (p == "-t")	{
			i++;
			timeStop = (var_t)atof(argv[i]) * Constants::YearToDay;
		}
		// Time step
		else if (p == "-dt") {
			i++;
			dt = (var_t)atof(argv[i]);
		}
		// Radius buffer factor
		else if (p == "-b") {
			i++;
			buffer_radius = (var_t)atof(argv[i]);
		}
		// Print out period
		else if (p == "-p")	{
			printout = true;
			i++;
			printoutPeriod = (var_t)atof(argv[i]) * Constants::YearToDay;
			i++;
			printoutLength = (var_t)atof(argv[i]) * Constants::YearToDay;
			i++;
			printoutStep = (var_t)atof(argv[i]) * Constants::YearToDay;
		}
		// Print-out location
		else if (p == "-o")	{
			i++;
			printoutDir = argv[i];
			printoutToFile = true;
		}
		// Input file
		else if (p == "-f")	{
			i++;
			filename = argv[i];
			file = true;
		}
		else if (p == "-r")	{
			random = true;
		}
		else {
			throw nbody_exception("Invalid switch on command-line.");
		}
		i++;
	}
}

void options::initial_condition(nbody* nb)
{
	vec_t*	coor = (vec_t*)nb->h_y[0].data();
	vec_t*	velo = (vec_t*)nb->h_y[1].data();
	nbody::param_t* param = (nbody::param_t*)nb->h_p.data();

	int i = 0;
	// Set the initial conditions
	{
		// Star
		param[i].mass = 1.0;		// M_sun
		param[i].radius = 0.0014;	// AU

		coor[i].x = 0.0;			// AU
		coor[i].y = 0.0;
		coor[i].z = 0.0;
		coor[i].w = 0.0;

		velo[i].x = 0.0;			// AU / day
		velo[i].y = 0.0;
		velo[i].z = 0.0;
		velo[i].w = 0.0;

		// Planet 0, 1, 2, ..., n
		for (i = 1; i < this->n; i++) {
			param[i].mass = 1.0e-7;			// M_sun
			param[i].radius = 0.000014;		// AU

			coor[i].x = 1.0 + (i-1)*0.1;	// AU
			coor[i].y = 0.0;
			coor[i].z = 0.0;
			coor[i].w = 0.0;

			// Compute the circular velocity for planet 0
			var_t	r = sqrt( SQR(coor[i].x - coor[0].x) + SQR(coor[i].y - coor[0].y) + SQR(coor[i].z - coor[0].z) );
			var_t	v = K * sqrt(param[0].mass / r);

			velo[i].x = 0.0;			// AU / day
			velo[i].y = v;
			velo[i].z = 0.0;
			velo[i].w = 0.0;
		}
	}

	// Transform the variables to the barycentric system
	{
		// Compute the total mass of the system
		var_t totalMass = 0.0;
		for (int j = 0; j < this->n; j++ ) {
			totalMass += param[j].mass;
		}

		// Position and velocity of the system's barycenter
		vec_t R0 = {0.0, 0.0, 0.0, 0.0};
		vec_t V0 = {0.0, 0.0, 0.0, 0.0};

		// Compute the position and velocity of the barycenter of the system
		for (int j = 0; j < this->n; j++ ) {
			R0.x += param[j].mass * coor[j].x;
			R0.y += param[j].mass * coor[j].y;
			R0.z += param[j].mass * coor[j].z;

			V0.x += param[j].mass * velo[j].x;
			V0.y += param[j].mass * velo[j].y;
			V0.z += param[j].mass * velo[j].z;
		}
		R0.x /= totalMass;
		R0.y /= totalMass;
		R0.z /= totalMass;

		V0.x /= totalMass;
		V0.y /= totalMass;
		V0.z /= totalMass;

		// Transform the bodies coordinates and velocities
		for (int j = 0; j < this->n; j++ ) {
			coor[j].x -= R0.x;
			coor[j].y -= R0.y;
			coor[j].z -= R0.z;

			velo[j].x -= V0.x;
			velo[j].y -= V0.y;
			velo[j].z -= V0.z;
		}
	}
}

ode* options::create_ode()
{
	nbody* nb = new nbody(n, timeStart);

	nb->t = timeStart;
	
	if (file) {
		nb->load(filename, n);
	}

	nb->copy_to_device();

	return nb;
}

nbody*	options::create_nbody()
{
	nbody*	nb = new nbody(n, timeStart);

	nb->t = timeStart;

	if (file) {
		nb->load(filename, n);
	}
	else {
		initial_condition(nb);
	}
	nb->copy_to_device();

	return nb;
}

pp_disk*	options::create_pp_disk()
{
	pp_disk *ppd = new pp_disk(nBodies, gasDisk, timeStart);

	ppd->t = timeStart;

	if (file) {
		ppd->load(filename, nBodies->total);
	}
	else {
		;
	}
	ppd->copy_to_device();

	return ppd;
}

integrator* options::create_integrator(ode* f)
{
	integrator* intgr;

	switch (inttype)
	{
	case INTEGRATOR_EULER:
		intgr = new euler(*f, dt);
		break;
	case INTEGRATOR_RUNGEKUTTA2:
		intgr = new rungekutta<2>(*f, dt, adaptive, tolerance);
		break;
	case INTEGRATOR_OPT_RUNGEKUTTA2:
		intgr = new midpoint(*f, dt, adaptive, tolerance);
		break;
	case INTEGRATOR_RUNGEKUTTA4:
		intgr = new rungekutta<4>(*f, dt, adaptive, tolerance);
		break;
	case INTEGRATOR_OPT_RUNGEKUTTA4:
		intgr = new rk4(*f, dt, adaptive, tolerance);
		break;
	case INTEGRATOR_RUNGEKUTTANYSTROM:
		intgr = new rungekuttanystrom<9>(*f, dt, adaptive, tolerance);
		break;
	case INTEGRATOR_OPT_RUNGEKUTTANYSTROM:
		intgr = new rkn76(*f, dt, adaptive, tolerance);
		break;
	default:
		throw nbody_exception("Requested integrator is not implemented.");
	}

	return intgr;
}
