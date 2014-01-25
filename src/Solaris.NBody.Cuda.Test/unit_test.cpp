#define _CRT_SECURE_NO_WARNING

// includes, system 
#include <iostream>

// includes CUDA
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

// includes Thrust
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/reduce.h>

//#include "config.h"
//#include "planets.h"
#include "Constants.h"
#include "gas_disc.h"


using namespace std;


//! Type of time variables
typedef double		ttt_t;
//! Type of variables
typedef double		var_t;
//! Type of tuple
typedef double2		var2_t;
//! Type of vectors
typedef double4		vec_t;
//! Type of boolean variables
typedef bool		bool_t;
//! Type of integer variables
typedef int			int_t;
//! Type of integer tuples variables
typedef int2		int2_t;

typedef struct param
{
	//! Mass of body in M_sol
	var_t mass;
	//! Radius of body in AU
	var_t radius;
	//! Density of body in M_sol AU-3
	var_t density;
	//! Used for the drag force  TODO
	var_t gamma_stokes;
	//! Used for the drag force  TODO
	var_t gamma_epstein;
} planets_param_t;

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

//typedef struct gaspar
//{
//	var2_t	rho;
//	var2_t	sch;
//	var2_t	eta;
//	var2_t	tau;
//} gaspar_t;

#define K			(var_t)0.01720209895
#define K2			(var_t)0.0002959122082855911025

#define	PI			(var_t)3.1415926535897932384626
#define	TWOPI		(var_t)6.2831853071795864769253
#define	TORAD		(var_t)0.0174532925199432957692
#define TODEG		(var_t)57.295779513082320876798

// It must be enclosed in parentheses in order to give correct results in
// the case of a division i.e. 1/SQR(x) -> 1/((x)*(x))
#define	SQR(x)		((x)*(x))
#define	CUBE(x)		((x)*(x)*(x))

void shift_into_range(var_t lower, var_t upper, var_t* value)
{
    var_t range = upper - lower;
    while (upper <= *value) {
        *value -= range;
    }
    while (lower > *value) {
        *value += range;
    }
}

vec_t	cross_product(const vec_t* v, const vec_t* u)
{
	vec_t result;

	result.x = v->y*u->z - v->z*u->y;
    result.y = v->z*u->x - v->x*u->z;
    result.z = v->x*u->y - v->y*u->x;

	return result;
}

var_t	dot_product(const vec_t* v, const vec_t* u)
{
	return v->x * u->x + v->y * u->y + v->z * u->z;
}

var_t	norm2(const vec_t* v)
{
	return SQR(v->x) + SQR(v->y) + SQR(v->z);
}

var_t	norm(const vec_t* v)
{
	return sqrt(norm2(v));
}

vec_t	circular_velocity(var_t mu, const vec_t* rVec)
{
	vec_t result = {0.0, 0.0, 0.0, 0.0};

	var_t r		= sqrt(SQR(rVec->x) + SQR(rVec->y));
	var_t vc	= sqrt(mu/r);

	var_t p;
	if (rVec->x == 0.0 && rVec->y == 0.0) {
		return result;
	}
	else if (rVec->y == 0.0) {
		result.y = rVec->x > 0.0 ? vc : -vc;
	}
	else if (rVec->x == 0.0) {
		result.x = rVec->y > 0.0 ? -vc : vc;
	}
	else if (rVec->x >= rVec->y) {
		p = rVec->y / rVec->x;
		result.y = rVec->x >= 0 ? vc/sqrt(1.0 + SQR(p)) : -vc/sqrt(1.0 + SQR(p));
		result.x = -result.y*p;
	}
	else {
		p = rVec->x / rVec->y;
		result.x = rVec->y >= 0 ? -vc/sqrt(1.0 + SQR(p)) : vc/sqrt(1.0 + SQR(p));
		result.y = -result.x*p;
	}

	return result;
}

vec_t	gas_velocity(var2_t eta, var_t mu, const vec_t* rVec)
{
	vec_t result = circular_velocity(mu, rVec);
	var_t r		= sqrt(SQR(rVec->x) + SQR(rVec->y));

	var_t v		 = sqrt(1.0 - 2.0*eta.x * pow(r, eta.y));
	result.x	*= v;
	result.y	*= v;
	
	return result;
}

// TODO: implemet INNER_EDGE to get it from the input
#define INNER_EDGE 0.1 // AU
var_t	gas_density_at(const gas_disc* gasDisc, const vec_t* rVec)
{
	var_t result = 0.0;

	var_t r		= sqrt(SQR(rVec->x) + SQR(rVec->y));
	var_t h		= gasDisc->sch.x * pow(r, gasDisc->sch.y);
	var_t arg	= SQR(rVec->z/h);
	if (INNER_EDGE < r) {
		result	= gasDisc->rho.x * pow(r, gasDisc->rho.y) * exp(-arg);
	}
	else {
		var_t a	= gasDisc->rho.x * pow(INNER_EDGE, gasDisc->rho.y - 4.0);
		result	= a * SQR(SQR(r)) * exp(-arg);
	}

	return result;
}
#undef INNER_EDGE

var_t	calculate_kinetic_energy(const vec_t* vVec)
{
	return 0.5 * norm2(vVec);
}

var_t	calculate_potential_energy(var_t mu, const vec_t* rVec)
{
	return -mu / norm(rVec);
}

var_t	calculate_energy(var_t mu, const vec_t* rVec, const vec_t* vVec)
{
	return calculate_kinetic_energy(vVec) + calculate_potential_energy(mu, rVec);
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

#define	sq3	1.0e-14
int_t	calculate_sma_ecc(var_t mu, const vec_t* rVec, const vec_t* vVec, var_t* sma, var_t* ecc)
{
	// Calculate energy, h
    var_t h = calculate_energy(mu, rVec, vVec);
    if (h >= 0.0) {
        return 1;
    }

	// Calculate semi-major axis, a
    *sma = -mu / (2.0 * h);

    vec_t cVec = cross_product(rVec, vVec);
	cVec.w = norm2(&cVec);		// cVec.w = c2

	// Calculate eccentricity, e
    var_t e2 = 1.0 + 2.0 * h * cVec.w / SQR(mu);
	*ecc = fabs(e2) < sq3 ? 0.0 : sqrt(e2); 

    return 0;
}
#undef	sq3

#define	sq2 1.0e-14
#define	sq3	1.0e-14
int_t	calculate_orbelem(var_t mu, const vec_t* rVec, const vec_t* vVec, orbelem_t* oe)
{
	// Calculate energy, h
    var_t h = calculate_energy(mu, rVec, vVec);
    if (h >= 0.0) {
        return 1;
    }

	var_t r = norm(rVec);
	var_t v = norm(vVec);

	vec_t cVec	= cross_product(rVec, vVec);
	vec_t vxc	= cross_product(vVec, &cVec);
	vec_t lVec;
	lVec.x		= -mu/r * rVec->x + vxc.x;
	lVec.y		= -mu/r * rVec->y + vxc.y;
	lVec.z		= -mu/r * rVec->z + vxc.z;
	lVec.w		= norm(&lVec);

	cVec.w = norm2(&cVec);		// cVec.w = c2
    
    // Calculate eccentricity, e
	var_t ecc = 1.0 + 2.0 * h * cVec.w / SQR(mu);
	ecc = abs(ecc) < sq3 ? 0.0 : sqrt(ecc); 

	// Calculate semi-major axis, a
    var_t sma = -mu / (2.0 * h);

    // Calculate inclination, incl
	cVec.w = sqrt(cVec.w);		// cVec.w = c
    var_t cosi = cVec.z / cVec.w;
    var_t sini = sqrt(SQR(cVec.x) + SQR(cVec.y)) / cVec.w;
    var_t incl = acos(cosi);
    if (incl < sq2) {
        incl = 0.0;
    }
    
    // Calculate longitude of node, O
    var_t node = 0.0;
    if (incl != 0.0) {
		var_t tmpx = -cVec.y / (cVec.w * sini);
        var_t tmpy =  cVec.x / (cVec.w * sini);
		node = atan2(tmpy, tmpx);
		shift_into_range(0.0, 2.0*PI, &node);
    }
    
    // Calculate argument of pericenter, w
    var_t E		= 0.0;
    var_t peri	= 0.0;
    if (ecc != 0.0) {
		var_t tmpx = ( lVec.x * cos(node) + lVec.y * sin(node)) / lVec.w;
        var_t tmpy = (-lVec.x * sin(node) + lVec.y * cos(node)) /(lVec.w * cosi);
        peri = atan2(tmpy, tmpx);
        shift_into_range(0.0, 2.0*PI, &peri);

        tmpx = 1.0 / ecc * (1.0 - r / sma);
		tmpy = dot_product(rVec, vVec) / (sqrt(mu * sma) * ecc);
        E = atan2(tmpy, tmpx);
        shift_into_range(0.0, 2.0*PI, &E);
    }
    else {
        peri = 0.0;
        E = atan2(rVec->y, rVec->x);
        shift_into_range(0.0, 2.0*PI, &E);
    }
    
    // Calculate mean anomaly, M
    var_t M = E - ecc * sin(E);
    shift_into_range(0.0, 2.0*PI, &M);

	oe->sma	= sma;
	oe->ecc	= ecc;
	oe->inc	= incl;
	oe->peri= peri;
	oe->node= node;
	oe->mean= M;

	return 0;
}
#undef	sq2
#undef	sq3

var_t	orbital_period(var_t mu, var_t sma) 
{
	return TWOPI * sqrt(CUBE(sma)/mu);
}

var_t	orbital_frequency(var_t mu, var_t sma) 
{
	return 1.0 / orbital_period(mu, sma);;
}


int	unit_test_of_nbody_util()
{
	bool	succeeded = true;
	char	func_name[256];
	char	err_msg[1024];

	cout << "The unit_test_of_nbody_util() started.\n\nThe unit test of the" << endl;

	{
		bool	failed = false;
		strcpy(func_name, "shift_into_range");

		var_t value = -1.0;
		shift_into_range(0.0, 1.0, &value);
		if (0.0 != value) {
			sprintf(err_msg, "\t%30s() function failed at line %d.", func_name, __LINE__);
			cerr << err_msg << endl;
			failed = true;
		}

		value = 0.0;
		shift_into_range(0.0, 1.0, &value);
		if (0.0 != value) {
			sprintf(err_msg, "\t%30s() function failed at line %d.", func_name, __LINE__);
			cerr << err_msg << endl;
			failed = true;
		}

		value = 1.0;
		shift_into_range(0.0, 1.0, &value);
		if (0.0 != value) {
			sprintf(err_msg, "\t%30s() function failed at line %d.", func_name, __LINE__);
			cerr << err_msg << endl;
			failed = true;
		}

		value = 2.0;
		shift_into_range(0.0, 1.0, &value);
		if (0.0 != value) {
			sprintf(err_msg, "\t%30s() function failed at line %d.", func_name, __LINE__);
			cerr << err_msg << endl;
			failed = true;
		}

		value = 1.0;
		shift_into_range(2.0, 5.0, &value);
		if (4.0 != value) {
			sprintf(err_msg, "\t%30s() function failed at line %d.", func_name, __LINE__);
			cerr << err_msg << endl;
			failed = true;
		}

		if (!failed) {
			sprintf(err_msg, "\t%30s() function passed.", func_name);
			cout << err_msg << endl;
		} 
		else {
			succeeded = false;
		}
	}

	{
		bool	failed = false;
		strcpy(func_name, "cross_product");

		vec_t a = {1, 2, 3, 0};
		vec_t b = {2, 2, 1, 0};
		vec_t result = {0, 0, 0, 0};
		result = cross_product(&a, &b);
		if (-4.0 != result.x) {
			sprintf(err_msg, "\t%30s() function failed at line %d.", func_name, __LINE__);
			cerr << err_msg << endl;
			failed = true;
		}
		if (5.0 != result.y) {
			sprintf(err_msg, "\t%30s() function failed at line %d.", func_name, __LINE__);
			cerr << err_msg << endl;
			failed = true;
		}
		if (-2.0 != result.z) {
			sprintf(err_msg, "\t%30s() function failed at line %d.", func_name, __LINE__);
			cerr << err_msg << endl;
			failed = true;
		}

		if (!failed) {
			sprintf(err_msg, "\t%30s() function passed.", func_name);
			cout << err_msg << endl;
		}
		else {
			succeeded = false;
		}
	}

	{
		bool	failed = false;
		strcpy(func_name, "dot_product");

		vec_t a = {1, 2, 3, 0};
		vec_t b = {2, 2, 1, 0};

		var_t result = dot_product(&a, &b);
		if (9.0 != result) {
			sprintf(err_msg, "\t%30s() function failed at line %d.", func_name, __LINE__);
			cerr << err_msg << endl;
			failed = true;
		}

		b.x = -4.0; b.y = 5.0; b.z = -2.0;
		result = dot_product(&a, &b);
		if (0.0 != result) {
			sprintf(err_msg, "\t%30s() function failed at line %d.", func_name, __LINE__);
			cerr << err_msg << endl;
			failed = true;
		}

		if (!failed) {
			sprintf(err_msg, "\t%30s() function passed.", func_name);
			cout << err_msg << endl;
		}
		else {
			succeeded = false;
		}
	}
	
	{
		bool	failed = false;
		strcpy(func_name, "norm2");

		vec_t a = {1, 2, 3, 0};

		var_t result = norm2(&a);
		if (14.0 != result) {
			sprintf(err_msg, "\t%30s() function failed at line %d.", func_name, __LINE__);
			cerr << err_msg << endl;
			failed = true;
		}

		if (!failed) {
			sprintf(err_msg, "\t%30s() function passed.", func_name);
			cout << err_msg << endl;
		}
		else {
			succeeded = false;
		}
	}
	
	{
		bool	failed = false;
		strcpy(func_name, "norm");

		vec_t a = {1, 2, 3, 0};

		var_t result = norm(&a);
		if (3.7416573867739413855837487323165 != result) {
			sprintf(err_msg, "\t%30s() function failed at line %d.", func_name, __LINE__);
			cerr << err_msg << endl;
			failed = true;
		}

		if (!failed) {
			sprintf(err_msg, "\t%30s() function passed.", func_name);
			cout << err_msg << endl;
		}
		else {
			succeeded = false;
		}
	}
		
	{
		bool	failed = false;
		strcpy(func_name, "circular_velocity");

		vec_t result;

		var_t mu = 1.0;
		vec_t rVec = {1.0, 0.0, 0.0, 0.0};
		result = circular_velocity(mu, &rVec);
		if (0.0 != result.x || 1.0 != result.y || 0.0 != result.z) {
			sprintf(err_msg, "\t%30s() function failed at line %d.", func_name, __LINE__);
			cerr << err_msg << endl;
			failed = true;
		}
		rVec.x = 0.0; rVec.y = 1.0;
		result = circular_velocity(mu, &rVec);
		if (-1.0 != result.x || 0.0 != result.y || 0.0 != result.z) {
			sprintf(err_msg, "\t%30s() function failed at line %d.", func_name, __LINE__);
			cerr << err_msg << endl;
			failed = true;
		}
		rVec.x = -1.0; rVec.y = 0.0;
		result = circular_velocity(mu, &rVec);
		if (0.0 != result.x || -1.0 != result.y || 0.0 != result.z) {
			sprintf(err_msg, "\t%30s() function failed at line %d.", func_name, __LINE__);
			cerr << err_msg << endl;
			failed = true;
		}
		rVec.x = 0.0; rVec.y = -1.0;
		result = circular_velocity(mu, &rVec);
		if (1.0 != result.x || 0.0 != result.y || 0.0 != result.z) {
			sprintf(err_msg, "\t%30s() function failed at line %d.", func_name, __LINE__);
			cerr << err_msg << endl;
			failed = true;
		}

		mu = 0.001;
		rVec.x = 5.0; rVec.y = 0.5;
		result = circular_velocity(mu, &rVec);
		// vc = 0.01410699961161117446257797364463
		// v.y = 0.01403698925583099106511593011858
		// v.x = -0.001403698925583099106511593011858
		if (fabs(-0.0014036989255830 - result.x) > 1.0e-15 || fabs(0.01403698925583099 - result.y) > 1.0e-15 || 0.0 != result.z) {
			sprintf(err_msg, "\t%30s() function failed at line %d.", func_name, __LINE__);
			cerr << err_msg << endl;
			failed = true;
		}

		if (!failed) {
			sprintf(err_msg, "\t%30s() function passed.", func_name);
			cout << err_msg << endl;
		}
		else {
			succeeded = false;
		}
	}

	{
		bool	failed = false;
		strcpy(func_name, "gas_velocity");

		vec_t result;

		var_t mu = 1.0;
		vec_t rVec = {1.0, 0.0, 0.0, 0.0};
		var2_t eta = {1.0e-3, 0.5};
		result = gas_velocity(eta, mu, &rVec);
		if (0.0 != result.x || 0.99899949949937412368543414284205 != result.y || 0.0 != result.z) {
			sprintf(err_msg, "\t%30s() function failed at line %d.", func_name, __LINE__);
			cerr << err_msg << endl;
			failed = true;
		}

		if (!failed) {
			sprintf(err_msg, "\t%30s() function passed.", func_name);
			cout << err_msg << endl;
		}
		else {
			succeeded = false;
		}
	}

	{
		bool	failed = false;
		strcpy(func_name, "gas_density_at");

		var_t result;

		var_t mu = 1.0;
		gas_disc gasDisc;
		gasDisc.rho.x = 1.0e-10; // g/cm3
		gasDisc.rho.y = -3.0;
		gasDisc.sch.x = 5.0e-2;	// AU
		gasDisc.sch.y = 3.0/2.0;
		vec_t rVec = {0.1, 0.0, 0.0, 0.0};
		result = gas_density_at(&gasDisc, &rVec);
		if (fabs(1.0e-7 - result) > 1.0e-15) {
			sprintf(err_msg, "\t%30s() function failed at line %d.", func_name, __LINE__);
			cerr << err_msg << endl;
			failed = true;
		}

		rVec.x = 0.05; rVec.y = 0.0; rVec.z = 0.0; rVec.w = 0.0;
		result = gas_density_at(&gasDisc, &rVec);
		if (fabs(6.25e-9 - result) > 1.0e-15) {
			sprintf(err_msg, "\t%30s() function failed at line %d.", func_name, __LINE__);
			cerr << err_msg << endl;
			failed = true;
		}

		rVec.x = 1.0; rVec.y = 0.0; rVec.z = 0.0; rVec.w = 0.0;
		result = gas_density_at(&gasDisc, &rVec);
		if (1.0e-10 != result) {
			sprintf(err_msg, "\t%30s() function failed at line %d.", func_name, __LINE__);
			cerr << err_msg << endl;
			failed = true;
		}

		rVec.x = 0.0; rVec.y = 1.0; rVec.z = 0.0; rVec.w = 0.0;
		result = gas_density_at(&gasDisc, &rVec);
		if (1.0e-10 != result) {
			sprintf(err_msg, "\t%30s() function failed at line %d.", func_name, __LINE__);
			cerr << err_msg << endl;
			failed = true;
		}

		rVec.x = 1.0; rVec.y = 0.0; rVec.z = 5.0e-2; rVec.w = 0.0;
		result = gas_density_at(&gasDisc, &rVec);
		if (fabs(3.6787944117144232159552377016146e-11 - result) > 1.0e-16) {
			sprintf(err_msg, "\t%30s() function failed at line %d.", func_name, __LINE__);
			cerr << err_msg << endl;
			failed = true;
		}

		if (!failed) {
			sprintf(err_msg, "\t%30s() function passed.", func_name);
			cout << err_msg << endl;
		}
		else {
			succeeded = false;
		}
	}

	{
		bool	failed = false;
		strcpy(func_name, "calculate_potential_energy");

		var_t result;

		var_t mu = 1.0;
		vec_t rVec = {0.1, 0.0, 0.0, 0.0};
		result = calculate_potential_energy(mu, &rVec);
		if (-10.0 != result) {
			sprintf(err_msg, "\t%30s() function failed at line %d.", func_name, __LINE__);
			cerr << err_msg << endl;
			failed = true;
		}

		rVec.x = 3.0; rVec.y = 4.0; rVec.z = 0.0; rVec.w = 0.0;
		result = calculate_potential_energy(mu, &rVec);
		if (-0.2 != result) {
			sprintf(err_msg, "\t%30s() function failed at line %d.", func_name, __LINE__);
			cerr << err_msg << endl;
			failed = true;
		}

		if (!failed) {
			sprintf(err_msg, "\t%30s() function passed.", func_name);
			cout << err_msg << endl;
		}
		else {
			succeeded = false;
		}
	}

	{
		bool	failed = false;
		strcpy(func_name, "calculate_kinetic_energy");

		var_t result;

		vec_t rVec = {0.1, 0.0, 0.0, 0.0};
		result = calculate_kinetic_energy(&rVec);
		if (fabs(5.0e-3 - result) > 1.0e-16) {
			sprintf(err_msg, "\t%30s() function failed at line %d.", func_name, __LINE__);
			cerr << err_msg << endl;
			failed = true;
		}

		rVec.x = 3.0; rVec.y = 4.0; rVec.z = 0.0; rVec.w = 0.0;
		result = calculate_kinetic_energy(&rVec);
		if (fabs(12.5 - result) > 1.0e-16) {
			sprintf(err_msg, "\t%30s() function failed at line %d.", func_name, __LINE__);
			cerr << err_msg << endl;
			failed = true;
		}

		if (!failed) {
			sprintf(err_msg, "\t%30s() function passed.", func_name);
			cout << err_msg << endl;
		}
		else {
			succeeded = false;
		}
	}

	{
		bool	failed = false;
		strcpy(func_name, "kepler_equation_solver");

		var_t result = 0.0;
		var_t ecc = 0.5;
		var_t mean = 27.0 * TORAD;
		int_t ret_code = kepler_equation_solver(ecc, mean, 1.0e-10, &result);
		if (0 != ret_code || fabs(48.43417991487915 - result*TODEG) > 1.0e-10) {
			sprintf(err_msg, "\t%30s() function failed at line %d.", func_name, __LINE__);
			cerr << err_msg << endl;
			failed = true;
		}

		mean = 125.0 * TORAD;
		ret_code = kepler_equation_solver(ecc, mean, 1.0e-10, &result);
		if (0 != ret_code || fabs(142.4568442405949 - result*TODEG) > 1.0e-10) {
			sprintf(err_msg, "\t%30s() function failed at line %d.", func_name, __LINE__);
			cerr << err_msg << endl;
			failed = true;
		}

		mean = 205.0 * TORAD;
		ret_code = kepler_equation_solver(ecc, mean, 1.0e-10, &result);
		if (0 != ret_code || fabs(196.7457972979956 - result*TODEG) > 1.0e-10) {
			sprintf(err_msg, "\t%30s() function failed at line %d.", func_name, __LINE__);
			cerr << err_msg << endl;
			failed = true;
		}

		mean = 300.0 * TORAD;
		ret_code = kepler_equation_solver(ecc, mean, 1.0e-10, &result);
		if (0 != ret_code || fabs(271.36018243209764 - result*TODEG) > 1.0e-10) {
			sprintf(err_msg, "\t%30s() function failed at line %d.", func_name, __LINE__);
			cerr << err_msg << endl;
			failed = true;
		}

		mean = 350.0 * TORAD;
		ret_code = kepler_equation_solver(ecc, mean, 1.0e-10, &result);
		if (0 != ret_code || fabs(340.38113495327434 - result*TODEG) > 1.0e-10) {
			sprintf(err_msg, "\t%30s() function failed at line %d.", func_name, __LINE__);
			cerr << err_msg << endl;
			failed = true;
		}

		if (!failed) {
			sprintf(err_msg, "\t%30s() function passed.", func_name);
			cout << err_msg << endl;
		}
		else {
			succeeded = false;
		}
	}

	{
		bool	failed = false;
		strcpy(func_name, "calculate_phase");

		const var_t m0 = 1.0;
		const var_t m1 = 0.0;
		const var_t mu = K2*(m0 + m1);
		orbelem_t oe;
		oe.sma  =   1.944132;
		oe.ecc  =   0.0737774;
		oe.inc  =  22.50914 * TORAD;
		oe.peri = 124.04003 * TORAD;
		oe.node = 175.33712 * TORAD;
		oe.mean = 217.31447 * TORAD;

// Emesétől adatok:
//Epoch [JD]      a [AU]          e               i [deg]         peric [deg] Omega [deg]     tau [JD]        M [deg]         n [deg/day] 
//2.456600500E+06 1.944132000E+00 7.377740000E-02 2.250914000E+01 1.240400300E+02 1.753371200E+02 0.000000000E+00 2.173144700E+02 3.635929300E-01

//Epoch [JD]      X[AU]                  Y[AU]                  Z[AU] VX[AU/day] VY[AU/day]             VZ[AU/day]
//+2.45660050E+06 -1.824214307389804E+00 +9.081741252311457E-01 -3.136484471795240E-01 -4.959970348405931E-03 -9.644216216654788E-03 +4.150430491295794E-03

		vec_t rVec = {0.0, 0.0, 0.0, 0.0};
		vec_t vVec = {0.0, 0.0, 0.0, 0.0};

		int_t ret_code = calculate_phase(mu, &oe, &rVec, &vVec);
		if (0 != ret_code	|| fabs(-1.824214307389804 - rVec.x)  > 1.0e-14
							|| fabs(+0.9081741252311457 - rVec.y) > 1.0e-14
							|| fabs(-0.3136484471795240 - rVec.z) > 1.0e-14) {
			sprintf(err_msg, "\t%30s() function failed at line %d.", func_name, __LINE__);
			cerr << err_msg << endl;
			failed = true;
		}
		if ( 				   fabs(-4.959970348405931E-03 - vVec.x) > 1.0e-14
							|| fabs(-9.644216216654788E-03 - vVec.y) > 1.0e-14
							|| fabs(+4.150430491295794E-03 - vVec.z) > 1.0e-14) {
			sprintf(err_msg, "\t%30s() function failed at line %d.", func_name, __LINE__);
			cerr << err_msg << endl;
			failed = true;
		}
		if (!failed) {
			sprintf(err_msg, "\t%30s() function passed.", func_name);
			cout << err_msg << endl;
		}
		else {
			succeeded = false;
		}

		strcpy(func_name, "calculate_sma_ecc");
		ret_code = calculate_sma_ecc(mu, &rVec, &vVec, &oe.sma, &oe.ecc);
		if (0 != ret_code || fabs(1.944132 - oe.sma) > 1.0e-6 || fabs(0.0737774 - oe.ecc) > 1.0e-7) {
			sprintf(err_msg, "\t%30s() function failed at line %d.", func_name, __LINE__);
			cerr << err_msg << endl;
			failed = true;
		}

		if (!failed) {
			sprintf(err_msg, "\t%30s() function passed.", func_name);
			cout << err_msg << endl;
		}
		else {
			succeeded = false;
		}

		strcpy(func_name, "calculate_orbelem");
		ret_code = calculate_orbelem(mu, &rVec, &vVec, &oe);
		if (0 != ret_code	|| fabs(  1.944132000000000 - oe.sma) > 1.0e-14 
							|| fabs(  0.073777400000000 - oe.ecc) > 1.0e-14
							|| fabs( 22.509140000000000 - oe.inc * TODEG) > 1.0e-13
							|| fabs(124.040030000000000 - oe.peri* TODEG) > 1.0e-12
							|| fabs(175.337120000000000 - oe.node* TODEG) > 1.0e-12
							|| fabs(217.314470000000000 - oe.mean* TODEG) > 1.0e-12) {
			sprintf(err_msg, "\t%30s() function failed at line %d.", func_name, __LINE__);
			cerr << err_msg << endl;
			failed = true;
		}

		if (!failed) {
			sprintf(err_msg, "\t%30s() function passed.", func_name);
			cout << err_msg << endl;
		}
		else {
			succeeded = false;
		}
	}

	{
		bool	failed = false;
		strcpy(func_name, "orbital_period");

		const var_t m0 = 1.0;
		const var_t m1 = 1.0 / 3.3294605e5;
		const var_t mu = K2*(m0 + m1);

		var_t result = orbital_period(mu, 1.0);
		if (fabs(365.25635 - result) > 1.0e-5) {
			sprintf(err_msg, "\t%30s() function failed at line %d.", func_name, __LINE__);
			cerr << err_msg << endl;
			failed = true;
		}

		if (!failed) {
			sprintf(err_msg, "\t%30s() function passed.", func_name);
			cout << err_msg << endl;
		}
		else {
			succeeded = false;
		}
	}

	return succeeded ? 0 : 1;
}


int main(int argc, const char** argv)
{
	cudaError_t cudaStatus = cudaSuccess;
	int		result = 0;
	char	func_name[256];
	char	err_msg[1024];

	cout << "Solaris.NBody.Cuda.Test unit_test.cpp started\n\n";

	{
		strcpy(func_name, "unit_test_of_nbody_util");

		result = unit_test_of_nbody_util();
		if (0 == result) {
			sprintf(err_msg, "The unit test(s) of the %s() function passed.", func_name);
			cout << endl << err_msg << endl;
		}
		else {
			sprintf(err_msg, "The unit test(s) of the %s() function failed.", func_name);
			cout << endl << err_msg << endl;
		}
	}

}