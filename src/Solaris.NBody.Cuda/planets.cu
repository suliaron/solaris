// includes system 
#include <iostream>
#include <iomanip>
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
#include "number_of_bodies.h"
#include "nbody_exception.h"
#include "planets.h"

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


__host__ __device__
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

__host__ __device__
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
__host__ __device__
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

__host__ __device__
var_t	calculate_kinetic_energy(const vec_t* vVec)
{
	return 0.5 * norm2(vVec);
}

__host__ __device__
var_t	calculate_potential_energy(var_t mu, const vec_t* rVec)
{
	return -mu / norm(rVec);
}

__host__ __device__
var_t	calculate_energy(var_t mu, const vec_t* rVec, const vec_t* vVec)
{
	return calculate_kinetic_energy(vVec) + calculate_potential_energy(mu, rVec);
}

__host__ __device__
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

__host__ __device__
int_t	calculate_phase(var_t mu, const planets::orbelem_t* oe, vec_t* rVec, vec_t* vVec)
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
__host__ __device__
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
__host__ __device__
int_t	calculate_orbelem(var_t mu, const vec_t* rVec, const vec_t* vVec, planets::orbelem_t* oe)
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

__host__ __device__
var_t	orbital_period(var_t mu, var_t sma)
{
	return TWOPI * sqrt(CUBE(sma)/mu);
}

__host__ __device__
var_t	orbital_frequency(var_t mu, var_t sma) 
{
	return 1.0 / orbital_period(mu, sma);
}

__host__ __device__ var_t	calculate_gamma_stokes(var_t stokes, var_t density, var_t radius)
{
    return (3.0/8.0)*stokes/(density*radius);
}



// Calculate acceleration caused by particle j on particle i 
__host__ __device__ 
vec_t calculate_grav_accel_pair(const vec_t ci, const vec_t cj, var_t mass, vec_t a)
{
	vec_t d;
	
	d.x = cj.x - ci.x;
	d.y = cj.y - ci.y;
	d.z = cj.z - ci.z;

	d.w = SQR(d.x) + SQR(d.y) + SQR(d.z);
	d.w = d.w * d.w * d.w;
	d.w = K2 * mass / sqrt(d.w);

	a.x += d.x * d.w;
	a.y += d.y * d.w;
	a.z += d.z * d.w;

	return a;
}

__global__
void	calculate_grav_accel_kernel(interaction_bound iBound, const planets::param_t* params, const vec_t* coor, vec_t* acce)
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
}

__global__
void calculate_drag_accel_kernel(interaction_bound iBound, var_t timeF, const gas_disc* gasDisc, const planets::param_t* params, const vec_t* coor, const vec_t* velo, vec_t* acce)
{
	int tid		= blockIdx.x * blockDim.x + threadIdx.x;
	int	bodyIdx = iBound.sink.x + tid;

	if (bodyIdx < iBound.sink.y) {
		// TODO: ask Laci, why coor[bodyIdx] does not work?
		//vec_t rVec  = coor[bodyIdx];
		vec_t vGas	= gas_velocity(gasDisc->eta, K2*params[0].mass, (vec_t*)&coor[bodyIdx]);
		var_t rhoGas= gas_density_at(gasDisc, (vec_t*)&coor[bodyIdx]) * timeF;
		//printf("rhoGas: %.15lf\n", rhoGas);
		//printf("vGas: %.15lf, %.15lf, %.15lf\n", vGas.x, vGas.y, vGas.z);
		//printf("velo[%d]: %.15lf, %.15lf, %.15lf\n", bodyIdx, velo[bodyIdx].x, velo[bodyIdx].y, velo[bodyIdx].z);

		vec_t u;
		u.x	= velo[bodyIdx].x - vGas.x;
		u.y	= velo[bodyIdx].y - vGas.y;
		u.z	= velo[bodyIdx].z - vGas.z;

		//printf("u: %.15lf, %.15lf, %.15lf\n", u.x, u.y, u.z);

		var_t C		= 0.0;
		// TODO: implement the different regimes according to the mean free path of the gas molecules
		// Epstein-regime:
		{

		}
		// Stokes-regime:
		{
			//var_t uLength = norm(&u);
			C = params[bodyIdx].gamma_stokes * norm(&u) * rhoGas;
			//printf("params[bodyIdx].gamma_stokes: %.15lf\n", params[bodyIdx].gamma_stokes);
			//printf("norm(&u): %.15lf\n", norm(&u));
			//printf("C: %.15lf\n", C);
		}
		// Transition regime:
		{

		}

		acce[tid].x = -C * u.x;
		acce[tid].y = -C * u.y;
		acce[tid].z = -C * u.z;
		acce[tid].w = 0.0;

		//printf("acce[%d]: %.15lf, %.15lf, %.15lf\n", tid, acce[tid].x, acce[tid].y, acce[tid].z);
	}
}

// a = a + b
__global__
void add_two_vector(int_t n, var_t *a, const var_t *b)
{
	int tid = blockIdx.x * blockDim.x + threadIdx.x;
	if (n > tid) {
		a[tid] += b[tid];
	}
}

cudaError_t	planets::call_calculate_grav_accel_kernel(const planets::param_t* params, const vec_t* coor, vec_t* acce)
{
	cudaError_t cudaStatus = cudaSuccess;
	
	interaction_bound iBound = bodies.get_self_interacting();

	int		nBodyToCalculate = bodies.n_self_interacting();
	int		nThread = std::min(THREADS_PER_BLOCK, nBodyToCalculate);
	int		nBlock = (nBodyToCalculate + nThread - 1)/nThread;
	dim3	grid(nBlock);
	dim3	block(nThread);

	calculate_grav_accel_kernel<<<grid, block>>>(iBound, params, coor, acce);
	cudaStatus = HANDLE_ERROR(cudaGetLastError());
	if (cudaSuccess != cudaStatus) {
		throw nbody_exception("calculate_grav_accel_kernel launch failed", cudaStatus);
	}

	iBound = bodies.get_nonself_interacting();
	nBodyToCalculate = bodies.super_planetesimal + bodies.planetesimal;
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

	iBound = bodies.get_non_interacting();
	nBodyToCalculate = bodies.test_particle;
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

cudaError_t planets::call_calculate_drag_accel_kernel(ttt_t time, const gas_disc* gasDisc, const planets::param_t* params, const vec_t* coor, const vec_t* velo, vec_t* acce)
{
	cudaError_t cudaStatus = cudaSuccess;

	// TODO: calculate it using the value of the time
	var_t timeF = 1.0;
	
	interaction_bound iBound = bodies.get_bodies_gasdrag();

	int		nBodyToCalculate = bodies.super_planetesimal + bodies.planetesimal;
	int		nThread = std::min(THREADS_PER_BLOCK, nBodyToCalculate);
	int		nBlock = (nBodyToCalculate + nThread - 1)/nThread;
	dim3	grid(nBlock);
	dim3	block(nThread);

	calculate_drag_accel_kernel<<<grid, block>>>(iBound, timeF, gasDisc, params, coor, velo, acce);
	cudaStatus = HANDLE_ERROR(cudaGetLastError());
	if (cudaSuccess != cudaStatus) {
		throw nbody_exception("calculate_grav_accel_kernel launch failed", cudaStatus);
	}

	return cudaStatus;
}

planets::planets(number_of_bodies bodies, gas_disc* gasDisc) :
	ode(2),
	bodies(bodies),
	gasDisc(gasDisc),
	d_gasDisc(0),
	acceGasDrag(d_vec_t())
{
	//round_up_n();
	allocate_vectors();

}

planets::~planets()
{
	cudaFree(d_gasDisc);
	free(gasDisc);
}

void planets::allocate_vectors()
{
	// Parameters
	int	npar = sizeof(planets::param_t) / sizeof(var_t);
	h_p.resize(npar * bodies.total);

	// Aliases to coordinates and velocities
	int ndim = sizeof(vec_t) / sizeof(var_t);
	h_y[0].resize(ndim * bodies.total);
	h_y[1].resize(ndim * bodies.total);

	if (0 != gasDisc) {
		acceGasDrag.resize(ndim * (bodies.super_planetesimal + bodies.planetesimal));
		// TODO:
		// ask Laci, how to find out if there was an error during these 2 cuda function calls
		cudaMalloc((void**)&d_gasDisc, sizeof(gas_disc));
		cudaMemcpy(d_gasDisc, gasDisc, sizeof(gas_disc), cudaMemcpyHostToDevice );
	}
}

void planets::round_up_n()
{
	// Round up n to the number of bodies per tile
	int m = ((bodies.total + NTILE - 1) / NTILE) * NTILE;

	if (bodies.total != m) {
		cerr << "Number of bodies rounded up to " << m << endl;
	}
	bodies.total_rounded_up = m;
}

void planets::calculate_dy(int i, int r, ttt_t t, const d_var_t& p, const std::vector<d_var_t>& y, d_var_t& dy)
{
	//h_var_t coor = y[0];
	//for (int i = 0; i < coor.size(); i++)
	//	cout << "coor[" << i << "] = " << setprecision(20) << coor[i] << endl;
	//
	//h_var_t velo;
	//h_var_t acce;
	switch (i)
	{
	case 0:
		// Copy velocities from previous step
		thrust::copy(y[1].begin(), y[1].end(), dy.begin());
		//velo = dy;
		//for (int i = 0; i < coor.size(); i++)
		//	cout << "velo[" << i << "] = " << setprecision(20) << velo[i] << endl;
		break;
	case 1:
		// Calculate accelerations originated from gravity
		call_calculate_grav_accel_kernel((param_t*)p.data().get(), (vec_t*)d_y[0].data().get(), (vec_t*)dy.data().get());
		//acce = dy;
		//for (int i = 0; i < coor.size(); i++)
		//	cout << "acce[" << i << "] = " << setprecision(20) << acce[i] << endl;

		if (0 != gasDisc) {
			if (0 == r) {
			// Calculate accelerations originated from gas drag
				call_calculate_drag_accel_kernel(t, d_gasDisc, (param_t*)p.data().get(), (vec_t*)d_y[0].data().get(), (vec_t*)d_y[1].data().get(), acceGasDrag.data().get());
			}
			// Add acceGasDrag to dy
			int_t offset = bodies.n_self_interacting() * 4;
			int nBodyToCalculate = bodies.super_planetesimal + bodies.planetesimal;
			int_t nData = (nBodyToCalculate) * 4;

			var_t* a = (var_t*)(dy.data().get() + offset);
			var_t* aGD = (var_t*)acceGasDrag.data().get();

			int		nThread = std::min(THREADS_PER_BLOCK, nBodyToCalculate);
			int		nBlock = (nBodyToCalculate + nThread - 1)/nThread;
			dim3	grid(nBlock);
			dim3	block(nThread);

			add_two_vector<<<grid, block>>>(nData, a, aGD);
			cudaError_t cudaStatus = cudaGetLastError();
			if (cudaStatus != cudaSuccess) {
				throw nbody_exception("add_two_vector kernel failed", cudaStatus);
			}
		}

		// Calculate accelerations originated from migration

		// Calculate accelerations originated from other forces

		break;
	}
}

void planets::load(string filename)
{
	if (bodies.total > this->bodies.total) {
		throw nbody_exception("Too many lines in file.");
	}

	vec_t* h_coord = (vec_t*)h_y[0].data();
	vec_t* h_veloc = (vec_t*)h_y[1].data();
	param_t* h_param = (param_t*)h_p.data();

	ifstream input(filename.c_str());

	if (input) {
        int		id;
		ttt_t	time;
        
		for (int i = 0; i < this->bodies.total; i++) { 
            input >> id;
			input >> time;

			input >> h_param[i].mass;
			input >> h_param[i].radius;

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
int planets::print_positions(ostream& sout)
{
	param_t* h_param = (param_t*)h_p.data();
	vec_t* h_coord = (vec_t*)h_y[0].data();
	vec_t* h_veloc = (vec_t*)h_y[1].data();
	
	for (int i = 0; i < bodies.total; i ++)
	{
		sout << i << '\t';
		sout << t << '\t';
		sout << h_param[i].mass << '\t';
		sout << h_param[i].radius << '\t';

		sout << h_coord[i].x << '\t';
		sout << h_coord[i].y << '\t';
		sout << h_coord[i].z << '\t';
		sout << h_veloc[i].x << '\t';
		sout << h_veloc[i].y << '\t';
		sout << h_veloc[i].z;

		sout << endl;
	}

	return bodies.total;
}
