// includes, system 
#include <cmath>
#include <ctime>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>

#include "Constants.h"

#define PI 3.1415926535897932384626
#define K2 0.0002959122082855911025

#define MASS_STAR		1.0			// M_sol
#define MASS_JUPITER	1.0e-3		// M_sol
#define RAD_STAR		0.005		// AU

#define MASS_SUN		1.9891E+30	// kg
#define MASS_FACTOR		5e-7		// M_sol
#define MASS_MU			log(4.0)
#define MASS_S			0.3
#define MASS_MIN		1.0e-20		// M_sol
#define MASS_MAX		1.0e-19		// M_sol

#define DIST_MIN		4.5			// AU
#define DIST_MAX		15			// AU

#define DENSITY			3000.0		// kg m-3
#define AU				149.6e9		// m

// It must be enclosed in parentheses in order to give correct results in
// the case of a division i.e. 1/SQR(x) -> 1/((x)*(x))
#define	SQR(x)			((x)*(x))
#define	CUBE(x)			((x)*(x)*(x))

using namespace std;

typedef double	var_t;
typedef int		int_t;

typedef struct orbelem
{
	//! Semimajor-axis of the body
	var_t sma;
	//! Eccentricity of the body
	var_t ecc;
	//! Inclination of the body
	var_t inc;
	//! Argument of the pericenter
	var_t peri;
	//! Longitude of the ascending node
	var_t node;
	//! Mean anomaly
	var_t mean;
} orbelem_t;

typedef struct double4 
{
	var_t	x;
	var_t	y;
	var_t	z;
	var_t	w;
} vec_t;

// Draw a number from a given distribution
var_t generate_random(var_t xmin, var_t xmax, var_t p(var_t))
{
	var_t x;
	var_t y;

	do
	{
		x = xmin + (var_t)rand() / RAND_MAX * (xmax - xmin);
		y = (var_t)rand() / RAND_MAX;
	}
	while (y > p(x));

	return x;
}

var_t pdf_mass_lognormal(var_t x)
{
	return 1.0 / sqrt(2 * PI) / MASS_S * exp(-pow(log(x) - MASS_MU, 2) / 2 / MASS_S / MASS_MU);
}

var_t pdf_distance_squared(var_t d)
{
	return d * d / DIST_MAX / DIST_MAX;
}

var_t pdf_distance_exp(var_t d)
{
	return exp(-d) * d * d;
}

var_t pdf_const(var_t x)
{
	return 1;
}

void calculate_circle_coord(var_t d, var_t phi, var_t* x, var_t* y)
{
	*x = d * cos(phi);
	*y = d * sin(phi);
}

void calculate_circle_veloc(var_t d, var_t phi, var_t* vx, var_t* vy)
{
	var_t v = sqrt(K2 * MASS_STAR / d);
	
	*vx = v * sin(phi);
	*vy = - v * cos(phi);
}

var_t calculate_radius(var_t m)
{
	var_t V = m * MASS_SUN / DENSITY;	// m3
	V /= AU * AU * AU;		// AU3
	
	return pow(3.0 / 4.0 / PI * V, 1.0 / 3.0);
}

int_t	kepler_equation_solver(var_t ecc, var_t mean, var_t eps, var_t* E)
{
	if (ecc == 0.0 || mean == 0.0 || mean == PI) {
        *E = mean;
		return 0;
    }
    *E = mean + ecc * (sin(mean)) / (1.0 - sin(mean + ecc) + sin(mean));
    var_t E1 = 0.0;
    var_t error;
    int_t step = 0;
    do {
        E1 = *E - (*E - ecc * sin(*E) - mean) / (1.0 - ecc * cos(*E));
        error = fabs(E1 - *E);
        *E = E1;
    } while (error > eps && step++ <= 15);
	if (step > 15 ) {
		return 1;
	}

	return 0;
}

int_t	calculate_phase(var_t mu, const orbelem_t* oe, vec_t* rVec, vec_t* vVec)
{
    var_t ecc = oe->ecc;
	var_t E = 0.0;
	if (kepler_equation_solver(ecc, oe->mean, 1.0e-14, &E) == 1) {
		return 1;
	}
    var_t v = 2.0 * atan(sqrt((1.0 + ecc) / (1.0 - ecc)) * tan(E / 2.0));

    var_t p = oe->sma * (1.0 - SQR(ecc));
    var_t r = p / (1.0 + ecc * cos(v));
    var_t kszi = r * cos(v);
    var_t eta = r * sin(v);
    var_t vKszi = -sqrt(mu / p) * sin(v);
    var_t vEta = sqrt(mu / p) * (ecc + cos(v));

    var_t cw = cos(oe->peri);
    var_t sw = sin(oe->peri);
    var_t cO = cos(oe->node);
    var_t sO = sin(oe->node);
    var_t ci = cos(oe->inc);
    var_t si = sin(oe->inc);

    vec_t P;
	P.x = cw * cO - sw * sO * ci;
	P.y = cw * sO + sw * cO * ci;
	P.z = sw * si;
    vec_t Q;
	Q.x = -sw * cO - cw * sO * ci;
	Q.y = -sw * sO + cw * cO * ci;
	Q.z = cw * si;

	rVec->x = kszi * P.x + eta * Q.x;
	rVec->y = kszi * P.y + eta * Q.y;
	rVec->z = kszi * P.z + eta * Q.z;

	vVec->x = vKszi * P.x + vEta * Q.x;
	vVec->y = vKszi * P.y + vEta * Q.y;
	vVec->z = vKszi * P.z + vEta * Q.z;

	return 0;
}

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

int generate_nbody(string filename, int n)
{
	var_t d;
	var_t phi;
	var_t m;

	var_t x, y, z;
	var_t vx, vy, vz;
	var_t r;

	char sep = ' ';

	std::ofstream	output;
	output.open(filename, std::ios_base::app);

	// Output central mass
	output << 0 << sep;
	output << 0.0 << sep;
	output << MASS_STAR << sep << RAD_STAR << sep;
	output << 0.0 << sep << 0.0 << sep << 0.0 << sep;
	output << 0.0 << sep << 0.0 << sep << 0.0;
	output << endl;

	// Output planets
	for (int i = 1; i < n; i ++)
	{
		d = generate_random(DIST_MIN, DIST_MAX, pdf_distance_squared);
		phi = generate_random(0, 2*PI, pdf_const);
		m = 1.0e-19; //MASS_FACTOR * generate_random(MASS_MIN, MASS_MAX, pdf_mass_lognormal);

		calculate_circle_coord(d, phi, &x, &y);
		calculate_circle_veloc(d, phi, &vx, &vy);
		r = calculate_radius(m);
		z = vz = 0;

		output << i << sep;
		output << 0 << sep;
		output << m << sep << r << sep;
		output << x << sep << y << sep << z << sep;
		output << vx << sep << vy << sep << vz;
		output << endl;
	}

	return 0;
}

int generate_nbody_Rezso(string filename, int n)
{
	var_t m, m0, m1;
	var_t r;

	char sep = ' ';

	std::ofstream	output;
	output.open(filename, std::ios_base::app);

	// Output central mass
	m0 = 2.0 * MASS_STAR;
	output << 0 << sep;
	output << 0.0 << sep;
	output << m0 << sep << 1.2 * RAD_STAR << sep;
	output << 0.0 << sep << 0.0 << sep << 0.0 << sep;
	output << 0.0 << sep << 0.0 << sep << 0.0;
	output << endl;

	// Calculate phase of giant planet
	orbelem_t oe = {64.5, 0.0, 0.0, 0.0, 0.0, 0.0};
	vec_t rVec = {0.0, 0.0, 0.0, 0.0};
	vec_t vVec = {0.0, 0.0, 0.0, 0.0};
	m1 = 5.0 * MASS_JUPITER;
	var_t mu = K2*(m0 + m1);
	int_t ret_code = calculate_phase(mu, &oe, &rVec, &vVec);
	if (ret_code == 1) {
		cerr << "Could not calculate the phase." << endl;
		return ret_code;
	}

	// Output giant planet
	output << 1 << sep;
	output << 0.0 << sep;
	output << m1 << sep << 1.0e-1 * RAD_STAR << sep;
	output << rVec.x << sep << rVec.y << sep << rVec.z << sep;
	output << vVec.x << sep << vVec.y << sep << vVec.z << sep;
	output << endl;

	srand (time(0));
	// Output planets
	for (int i = 2; i < n; i ++)
	{
		oe.sma = generate_random(70.0, 270.0, pdf_const);
		oe.ecc = 0.0;
		oe.inc = atan(0.05); // tan(i) = h/r = 5.0e-2
		oe.peri = generate_random(0.0, 2.0*PI, pdf_const);
		oe.node = generate_random(0.0, 2.0*PI, pdf_const);
		oe.mean = generate_random(0.0, 2.0*PI, pdf_const);

		var_t mu = K2*(m0 + 0.0);
		int_t ret_code = calculate_phase(mu, &oe, &rVec, &vVec);
		if (ret_code == 1) {
			cerr << "Could not calculate the phase." << endl;
			return ret_code;
		}

		m = 0.0;
		r = 0.0;

		output << i << sep;
		output << 0 << sep;
		output << m << sep << r << sep;
		output << rVec.x << sep << rVec.y << sep << rVec.z << sep;
		output << vVec.x << sep << vVec.y << sep << vVec.z << sep;
		output << endl;
	}

	return 0;
}

int generate_2_body(string filename, int n)
{
	var_t m, m0, m1;
	var_t r;

	char sep = ' ';

	std::ofstream	output;
	output.open(filename, ios::trunc);

	// Output central mass
	m0 = 1.0;
	output << 0 << sep;
	output << 0.0 << sep;
	output << m0 << sep << Constants::SolarRadiusToAu << sep;
	output << 0.0 << sep << 0.0 << sep << 0.0 << sep;
	output << 0.0 << sep << 0.0 << sep << 0.0 << endl;

	// Calculate phase of super-planetesimal
	var_t s = Constants::DegreeToRadian;
	orbelem_t oe = {5.20336301, 0.04839266, 1.3053 * s, 274.1977 * s, 100.55615 * s, 19.65053 * s};
	vec_t rVec = {0.0, 0.0, 0.0, 0.0};
	vec_t vVec = {0.0, 0.0, 0.0, 0.0};
	m1 = Constants::JupiterToSolar;
	var_t mu = Constants::Gauss2*(m0 + m1);
	int_t ret_code = calculate_phase(mu, &oe, &rVec, &vVec);
	if (ret_code == 1) {
		cerr << "Could not calculate the phase." << endl;
		return ret_code;
	}

	// Output super-planetesimal
	output << 1 << sep;
	output << 0.0 << sep;
	output << m1 << sep << 1.0e-5 * RAD_STAR << sep;
	output << setprecision(16) << rVec.x << sep << setprecision(16) << rVec.y << sep << setprecision(16) << rVec.z << sep;
	output << setprecision(16) << vVec.x << sep << setprecision(16) << vVec.y << sep << setprecision(16) << vVec.z << sep;
	output << endl;

	return 0;
}

int main(int argc, char* argv[])
{
	int n = atoi(argv[1]);
	int ret_code = 0;

	const string baseDir = "C:\\Work\\Solaris.Cuda.TestRuns";
	const string subDir = "2_Body";
	string curDir = combine_path(baseDir, subDir);
	
	srand(time(NULL));

	//ret_code = generate_nbody("E:\\Work\\VSSolutions\\solaris\\src\\Solaris.NBody.Cuda.Test\\TestRun\\16384_Body\\PlanetGen02.txt", n);
	//ret_code = generate_nbody_Rezso("E:\\Work\\VSSolutions\\solaris\\src\\Solaris.NBody.Cuda.Test\\TestRun\\Rezso\\Rezso.txt", n);
	ret_code = generate_2_body(combine_path(curDir, "TwoBody.txt"), 2);

	return ret_code;
}
