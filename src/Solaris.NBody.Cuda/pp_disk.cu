// includes system 
#include <iostream>
#include <fstream>

// includes CUDA
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

// includes Thrust
#include "thrust\device_vector.h"
#include "thrust\host_vector.h"
#include "thrust\generate.h"
#include "thrust\copy.h"
#include "thrust\transform.h"

// includes project's header
#include "config.h"
#include "interaction_bound.h"
#include "gas_disc.h"
#include "nbody_exception.h"
#include "number_of_bodies.h"
#include "pp_disk.h"

using namespace std;

#define THREADS_PER_BLOCK	256

static cudaError_t HandleError(cudaError_t cudaStatus, const char *file, int line)
{
    if (cudaSuccess != cudaStatus) {
        printf( "%s in %s at line %d\n", cudaGetErrorString( cudaStatus ), file, line );
        return cudaStatus;
    }
	return cudaStatus;
}
#define HANDLE_ERROR(cudaStatus) (HandleError(cudaStatus, __FILE__, __LINE__))

static __host__ __device__
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

static __host__ __device__
vec_t	vector_subtract(const vec_t* a, const vec_t* b)
{
	vec_t result;

	result.x = a->x - b->x;
    result.y = a->y - b->y;
    result.z = a->z - b->z;

	return result;
}

static __host__ __device__
vec_t	cross_product(const vec_t* v, const vec_t* u)
{
	vec_t result;

	result.x = v->y*u->z - v->z*u->y;
    result.y = v->z*u->x - v->x*u->z;
    result.z = v->x*u->y - v->y*u->x;

	return result;
}

static __host__ __device__
var_t	dot_product(const vec_t* v, const vec_t* u)
{
	return v->x * u->x + v->y * u->y + v->z * u->z;
}

static __host__ __device__
var_t	norm2(const vec_t* v)
{
	return SQR(v->x) + SQR(v->y) + SQR(v->z);
}

static __host__ __device__
var_t	norm(const vec_t* v)
{
	return sqrt(norm2(v));
}

static __host__ __device__
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

static __host__ __device__
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
static __host__ __device__
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

static __host__ __device__
var_t	calculate_kinetic_energy(const vec_t* vVec)
{
	return 0.5 * norm2(vVec);
}

static __host__ __device__
var_t	calculate_potential_energy(var_t mu, const vec_t* rVec)
{
	return -mu / norm(rVec);
}

static __host__ __device__
var_t	calculate_energy(var_t mu, const vec_t* rVec, const vec_t* vVec)
{
	return calculate_kinetic_energy(vVec) + calculate_potential_energy(mu, rVec);
}

static __host__ __device__
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

static __host__ __device__
int_t	calculate_phase(var_t mu, const pp_disk::orbelem_t* oe, vec_t* rVec, vec_t* vVec)
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
static __host__ __device__
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
static __host__ __device__
int_t	calculate_orbelem(var_t mu, const vec_t* rVec, const vec_t* vVec, pp_disk::orbelem_t* oe)
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

static __host__ __device__
var_t	orbital_period(var_t mu, var_t sma)
{
	return TWOPI * sqrt(CUBE(sma)/mu);
}

static __host__ __device__
var_t	orbital_frequency(var_t mu, var_t sma) 
{
	return 1.0 / orbital_period(mu, sma);
}

//static __host__ __device__
//var_t	calculate_gamma_stokes(var_t stokes, var_t density, var_t radius)
//{
//	if (density == 0.0 || radius == 0.0)
//		return 0.0;
//	else
//		return (3.0/8.0)*stokes/(density*radius);
//}


// Calculate acceleration caused by particle j on particle i 
static __host__ __device__ 
vec_t calculate_grav_accel_pair(const vec_t ci, const vec_t cj, var_t mass, vec_t ai)
{
	vec_t dVec;
	
	dVec.x = cj.x - ci.x;
	dVec.y = cj.y - ci.y;
	dVec.z = cj.z - ci.z;

	dVec.w = SQR(dVec.x) + SQR(dVec.y) + SQR(dVec.z);	// = r2
	var_t r = sqrt(dVec.w);								// = r

	dVec.w = mass / (r*dVec.w);

	ai.x += dVec.w * dVec.x;
	ai.y += dVec.w * dVec.y;
	ai.z += dVec.w * dVec.z;

	return ai;
}

static __global__
void	calculate_grav_accel_kernel(interaction_bound iBound, const pp_disk::param_t* params, const vec_t* coor, vec_t* acce)
{
	int	bodyIdx = iBound.sink.x + blockIdx.x * blockDim.x + threadIdx.x;

	if (bodyIdx < iBound.sink.y) {
		acce[bodyIdx].x = 0.0;
		acce[bodyIdx].y = 0.0;
		acce[bodyIdx].z = 0.0;
		acce[bodyIdx].w = 0.0;
		for (int j = iBound.source.x; j < iBound.source.y; j++) 
		{
			if (j == bodyIdx) {
				continue;
			}
			acce[bodyIdx] = calculate_grav_accel_pair(coor[bodyIdx], coor[j], params[j].mass, acce[bodyIdx]);
		}
	}
	acce[bodyIdx].x *= K2;
	acce[bodyIdx].y *= K2;
	acce[bodyIdx].z *= K2;
}

static __global__
void	calculate_orbelem_kernel(
		int_t total, int_t refBodyId, 
		const pp_disk::param_t *params, const vec_t *coor, const vec_t *velo, 
		pp_disk::orbelem_t *orbelem)
{
	int	bodyIdx = blockIdx.x * blockDim.x + threadIdx.x;
	if (total > bodyIdx && refBodyId != bodyIdx) {
		var_t mu = K2 * (params[refBodyId].mass + params[bodyIdx].mass);
		vec_t rVec = vector_subtract((vec_t*)(&coor[bodyIdx]), (vec_t*)(&coor[refBodyId]));
		vec_t vVec = vector_subtract((vec_t*)(&velo[bodyIdx]), (vec_t*)(&velo[refBodyId]));
		calculate_orbelem(mu, &rVec, &vVec, (pp_disk::orbelem_t*)(&orbelem[bodyIdx]));
	}
}

pp_disk::pp_disk(number_of_bodies *nBodies) :
	ode(2),
	nBodies(nBodies)
{
	//round_up_n();
	allocate_vectors();
}

pp_disk::~pp_disk()
{
	delete nBodies;
}

void pp_disk::allocate_vectors()
{
	// Parameters
	int	npar = sizeof(param_t) / sizeof(var_t);
	h_p.resize(npar * nBodies->total);

	// Aliases to coordinates and velocities
	int ndim = sizeof(vec_t) / sizeof(var_t);
	h_y[0].resize(ndim * nBodies->total);
	h_y[1].resize(ndim * nBodies->total);

	//if (0 != gasDisc) {
	//	acceGasDrag.resize(ndim * (bodies.super_planetesimal + bodies.planetesimal));
	//	// TODO:
	//	// ask Laci, how to find out if there was an error during these 2 cuda function calls
	//	cudaMalloc((void**)&d_gasDisc, sizeof(gas_disc));
	//	cudaMemcpy(d_gasDisc, gasDisc, sizeof(gas_disc), cudaMemcpyHostToDevice );
	//}
}

cudaError_t pp_disk::call_calculate_accel_kernel(const param_t* params, const vec_t* coor, vec_t* acce)
{
	cudaError_t cudaStatus = cudaSuccess;
	
	interaction_bound iBound = nBodies->get_self_interacting();
	int		nBodyToCalculate = nBodies->n_self_interacting();
	int		nThread = std::min(THREADS_PER_BLOCK, nBodyToCalculate);
	int		nBlock = (nBodyToCalculate + nThread - 1)/nThread;
	dim3	grid(nBlock);
	dim3	block(nThread);

	calculate_grav_accel_kernel<<<grid, block>>>(iBound, params, coor, acce);
	cudaStatus = HANDLE_ERROR(cudaGetLastError());
	if (cudaSuccess != cudaStatus) {
		throw nbody_exception("calculate_grav_accel_kernel launch failed", cudaStatus);
	}

	iBound = nBodies->get_nonself_interacting();
	nBodyToCalculate = nBodies->super_planetesimal + nBodies->planetesimal;
	if (0 < nBodyToCalculate) {
		nThread		= std::min(THREADS_PER_BLOCK, nBodyToCalculate);
		nBlock		= (nBodyToCalculate + nThread - 1)/nThread;
		grid.x		= nBlock;
		block.x		= nThread;

		calculate_grav_accel_kernel<<<grid, block>>>(iBound, params, coor, acce);
		cudaStatus = HANDLE_ERROR(cudaGetLastError());
		if (cudaSuccess != cudaStatus) {
			throw nbody_exception("calculate_grav_accel_kernel launch failed", cudaStatus);
		}
	}

	iBound = nBodies->get_non_interacting();
	nBodyToCalculate = nBodies->test_particle;
	if (0 < nBodyToCalculate) {
		nThread		= std::min(THREADS_PER_BLOCK, nBodyToCalculate);
		nBlock		= (nBodyToCalculate + nThread - 1)/nThread;
		grid.x		= nBlock;
		block.x		= nThread;

		calculate_grav_accel_kernel<<<grid, block>>>(iBound, params, coor, acce);
		cudaStatus = HANDLE_ERROR(cudaGetLastError());
		if (cudaSuccess != cudaStatus) {
			throw nbody_exception("calculate_grav_accel_kernel launch failed", cudaStatus);
		}
	}

	return cudaStatus;
}

void pp_disk::calculate_dy(int i, int r, ttt_t t, const d_var_t& p, const std::vector<d_var_t>& y, d_var_t& dy)
{
	switch (i)
	{
	case 0:
		// Copy velocities from previous step
		thrust::copy(y[1].begin(), y[1].end(), dy.begin());
		break;
	case 1:
		// Calculate accelerations
		param_t	*params = (param_t*)p.data().get();
		vec_t	*coor = (vec_t*)y[0].data().get();
		vec_t	*velo = (vec_t*)y[1].data().get();
		vec_t	*acce = (vec_t*)dy.data().get();
		call_calculate_accel_kernel(params, coor, acce);

		break;
	}
}

void pp_disk::calculate_orbelem(int_t refBodyId)
{
	static const int noe = sizeof(orbelem_t)/sizeof(var_t);

	cudaError_t cudaStatus = cudaSuccess;

	// There are noe orbital elements
	d_orbelem.resize(noe * nBodies->total);

	// Calculate orbital elements of the bodies
	int		nThread = std::min(THREADS_PER_BLOCK, nBodies->total);
	int		nBlock = (nBodies->total + nThread - 1)/nThread;
	dim3	grid(nBlock);
	dim3	block(nThread);

	param_t	*params = (param_t*)d_p.data().get();
	vec_t	*coor = (vec_t*)d_y[0].data().get();
	vec_t	*velo = (vec_t*)d_y[1].data().get();

	calculate_orbelem_kernel<<<grid, block>>>(nBodies->total, refBodyId, params, coor, velo, d_orbelem.data().get());
	cudaStatus = HANDLE_ERROR(cudaGetLastError());
	if (cudaSuccess != cudaStatus) {
		throw nbody_exception("calculate_orbelem_kernel launch failed", cudaStatus);
	}
}

void pp_disk::load(string filename, int n)
{
	vec_t* h_coord = (vec_t*)h_y[0].data();
	vec_t* h_veloc = (vec_t*)h_y[1].data();
	param_t* h_param = (param_t*)h_p.data();

	ifstream input(filename.c_str());

	var_t dummy;
	if (input) {
        int		id;
		ttt_t	time;
        
        for (int i = 0; i < n; i++) { 
            input >> id;
			input >> time;

			if (nBodies->n_massive() <= i) {
				input >> dummy;			// advance reader and discard mass
				input >> dummy;			// advance reader and discard radius
				h_param[i].mass = 0.0;
				h_param[i].radius = 0.0;
			}
			else {
				input >> h_param[i].mass;
				input >> h_param[i].radius;
			}

			input >> h_coord[i].x;
			input >> h_coord[i].y;
			input >> h_coord[i].z;

			input >> h_veloc[i].x;
			input >> h_veloc[i].y;
			input >> h_veloc[i].z;
        }
        input.close();
	}
	else {
		throw nbody_exception("Cannot open file.");
	}
}

// Print body positions
int pp_disk::print_positions(ostream& sout)
{
	param_t* h_param = (param_t*)h_p.data();
	vec_t* h_coord = (vec_t*)h_y[0].data();
	vec_t* h_veloc = (vec_t*)h_y[1].data();
	
	for (int i = 0; i < nBodies->total; i ++)
	{
		sout << i << '\t';
		sout << t << '\t';
		sout << h_param[i].mass << '\t';
		sout << h_param[i].radius << '\t';
		sout << h_coord[i].x - h_coord[0].x << '\t';
		sout << h_coord[i].y - h_coord[0].y << '\t';
		sout << h_coord[i].z - h_coord[0].z << '\t';
		sout << h_veloc[i].x - h_veloc[0].x << '\t';
		sout << h_veloc[i].y - h_veloc[0].y << '\t';
		sout << h_veloc[i].z - h_veloc[0].z;

		//sout << h_coord[i].x << '\t';
		//sout << h_coord[i].y << '\t';
		//sout << h_coord[i].z << '\t';
		//sout << h_veloc[i].x << '\t';
		//sout << h_veloc[i].y << '\t';
		//sout << h_veloc[i].z;

		sout << endl;
	}

	return 0;
}

// Print body orbital elements
int pp_disk::print_orbelem(ostream& sout)
{
	param_t *h_param = (param_t*)h_p.data();
	orbelem_t *oe	 = (orbelem_t*)h_orbelem.data();
	
	for (int i = 0; i < nBodies->total; i ++)
	{
		sout << i << '\t';
		sout << t << '\t';
		sout << h_param[i].mass << '\t';
		sout << h_param[i].radius << '\t';
		sout << oe[i].sma << '\t';
		sout << oe[i].ecc << '\t';
		sout << oe[i].inc << '\t';
		sout << oe[i].peri << '\t';
		sout << oe[i].node << '\t';
		sout << oe[i].mean;

		sout << endl;
	}

	return 0;
}
