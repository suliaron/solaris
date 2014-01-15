#include "nbody_util.h"

// Utility functions
__host__ __device__
void shift_into_range(var_t lower, var_t upper, var_t* value)
{
    var_t range = upper - lower;
    while (*value >= upper) {
        *value -= range;
    }
    while (*value < lower) {
        *value += range;
    }
}

// General 3D vector routines
__host__ __device__
vec_t	cross_product(const vec_t* v, const vec_t* u)
{
	vec_t result;

	result.x = v->y*u->z - v->z*u->y;
    result.y = v->z*u->x - v->x*u->z;
    result.z = v->x*u->y - v->y*u->x;

	return result;
}

__host__ __device__
var_t	dot_product(const vec_t* v, const vec_t* u)
{
	return v->x * u->x + v->y * u->y + v->z * u->z;
}

__host__ __device__
var_t	norm2(const vec_t* v)
{
	return SQR(v->x) + SQR(v->y) + SQR(v->z);
}

__host__ __device__
var_t	norm(const vec_t* v)
{
	return sqrt(norm2(v));
}

__host__ __device__ 
vec_t circular_velocity(var_t mu, var_t r, var_t alpha)
{
	vec_t	result;

	var_t v		= sqrt(mu/r);
	result.x	=-v*sin(alpha);
	result.y	= v*cos(alpha);
	result.z	= 0.0;

	return result;
}

__host__ __device__
vec_t gas_velocity(var2_t eta, var_t mu, var_t r, var_t alpha)
{
	vec_t result = circular_velocity(mu, r, alpha);

	var_t v		 = sqrt(1.0 - 2.0*eta.x * pow(r, eta.y));
	result.x	*= v;
	result.y	*= v;
	
	return result;
}

// TODO: implemet INNER_EDGE to get it from the input
#define INNER_EDGE 0.046 // AU ~ 10 R_sol
__host__ __device__
var_t	gas_density_at(const planets::gaspar_t* gaspar, var_t r, var_t z)
{
	var_t result = 0.0;

	var_t h		= gaspar->sch.x * pow(r, gaspar->sch.y);
	var_t arg	= SQR(z/h);
	if (r > INNER_EDGE) {
		result	= gaspar->rho.x * pow(r, gaspar->rho.y) * exp(-arg);
	}
	else {
		var_t a	= gaspar->rho.x * pow(INNER_EDGE, gaspar->rho.y - 4.0);
		result	= a * SQR(SQR(r)) * exp(-arg);
	}

	return result;
}
#undef INNER_EDGE

__host__ __device__
var_t	calculate_kinetic_energy(const vec_t* velo)
{
	return 0.5 * norm2(velo);
}

__host__ __device__
var_t	calculate_potential_energy(var_t mu, const vec_t* coor)
{
	return -mu / norm(coor);
}

__host__ __device__
var_t	calculate_energy(var_t mu, const vec_t* coor, const vec_t* velo)
{
	return calculate_kinetic_energy(velo) + calculate_potential_energy(mu, coor);
}

__host__ __device__
int_t kepler_equation_solver(var_t ecc, var_t mean, var_t eps, var_t* E)
{
	if (ecc == 0.0 || mean == 0.0 || mean == PI)
    {
        *E = mean;
		return 0;
    }
    *E = mean + ecc * (sin(mean)) / (1.0 - sin(mean + ecc) + sin(mean));
    var_t E1 = 0.0;
    var_t error;
    int_t step = 0;
    do {
        E1 = *E - (*E - ecc * sin(*E) - mean) / (1.0 - ecc * cos(*E));
        error = abs(E1 - *E);
        *E = E1;
        step++;
    } while (error > eps && step <= 25);
	if (step > 25 ) {
		return 1;
	}

	return 0;
}

__host__ __device__
int_t calculate_phase(var_t mu, const planets::orbelem_t* oe, vec_t* coor, vec_t* velo)
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

	coor->x = kszi * P.x + eta * Q.x;
	coor->y = kszi * P.y + eta * Q.y;
	coor->z = kszi * P.z + eta * Q.z;

	velo->x = vKszi * P.x + vEta * Q.x;
	velo->y = vKszi * P.y + vEta * Q.y;
	velo->z = vKszi * P.z + vEta * Q.z;

	return 0;
}

#define	sq3	1.0e-14
__host__ __device__
int_t	calculate_sma_ecc(var_t mu, const vec_t* coor, const vec_t* velo, var_t* sma, var_t* ecc)
{
	// Calculate energy, h
    var_t h = calculate_energy(mu, coor, velo);
    if (h >= 0.0) {
        return 1;
    }

	// Calculate semi-major axis, a
    *sma = -mu / (2.0 * h);

    vec_t cVec = cross_product(coor, velo);
	cVec.w = norm2(&cVec);		// cVec.w = c2

	// Calculate eccentricity, e
    *ecc = 1.0 + 2.0 * cVec.w * h / SQR(mu);
	*ecc = abs(*ecc) < sq3 ? 0.0 : sqrt(*ecc); 

    return 0;
}
#undef	sq3

#define	sq2 1.0e-14
#define	sq3	1.0e-14
__host__ __device__
int_t	calculate_orbelem(var_t mu, const vec_t* coor, const vec_t* velo, planets::orbelem_t* orbelem)
{
	// Calculate energy, h
    var_t h = calculate_energy(mu, coor, velo);
    if (h >= 0.0) {
        return 1;
    }

	var_t r = norm(coor);
	var_t v = norm(velo);

	vec_t cVec	= cross_product(coor, velo);
	vec_t vxc	= cross_product(velo, &cVec);
	vec_t lVec;
	lVec.x		= -mu/r * coor->x + vxc.x;
	lVec.y		= -mu/r * coor->y + vxc.y;
	lVec.z		= -mu/r * coor->z + vxc.z;
	lVec.w		= norm(&lVec);

	cVec.w = norm2(&cVec);		// cVec.w = c2
    
    // Calculate eccentricity, e
	var_t ecc = 1.0 + 2.0*cVec.w * h / SQR(mu);
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
		tmpy = dot_product(coor, velo) / (sqrt(mu * sma) * ecc);
        E = atan2(tmpy, tmpx);
        shift_into_range(0.0, 2.0*PI, &E);
    }
    else {
        peri = 0.0;
        E = atan2(coor->y, coor->x);
        shift_into_range(0.0, 2.0*PI, &E);
    }
    
    // Calculate mean anomaly, M
    var_t M = E - ecc * sin(E);
    shift_into_range(0.0, 2.0*PI, &M);

	orbelem->sma	= sma;
	orbelem->ecc	= ecc;
	orbelem->inc	= incl;
	orbelem->peri	= peri;
	orbelem->node	= node;
	orbelem->mean	= M;

	return 0;
}
#undef	sq2
#undef	sq3
