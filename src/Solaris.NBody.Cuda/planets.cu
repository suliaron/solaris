#include "config.h"
#include "InteractionBound.h"
#include "NumberOfBodies.h"
#include "planets.h"

#include "thrust\device_vector.h"
#include "thrust\host_vector.h"
#include "thrust\generate.h"
#include "thrust\copy.h"

#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#define THREADS_PER_BLOCK	256

// Calculate acceleration caused by particle j on parrticle i 
__device__ 
vec_t calculate_grav_accel_pair(const vec_t ci, const vec_t cj, var_t mass, vec_t a)
{
	vec_t d;
	
	d.x = cj.x - ci.x;
	d.y = cj.y - ci.y;
	d.z = cj.z - ci.z;

	d.w = d.x * d.x + d.y * d.y + d.z * d.z;
	d.w = d.w * d.w * d.w;
	d.w = - K2 * mass / sqrt(d.w);

	a.x += d.x * d.w;
	a.y += d.y * d.w;
	a.z += d.z * d.w;

	return a;
}

__device__ 
vec_t circular_velocity(var_t mu, var_t r, var_t alpha)
{
	vec_t	result;

	var_t v		= sqrt(mu/r);
	result.x	=-v*sin(alpha);
	result.y	= v*cos(alpha);
	result.z	= 0.0;

	return result;
}

__device__
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

__device__
var_t	gas_density_at(const planets::gaspar_t* gaspar, var_t r, var_t z)
{
	var_t result = 0.0;

	var_t	h	= gaspar->sch.x * pow(r, gaspar->sch.y);
	var_t	arg	= SQR(z/h);
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

__device__
var_t	calculate_kinetic_energy(const vec_t* velo)
{
	return 0.5 * (SQR(velo->x) + SQR(velo->y) + SQR(velo->z));
}

__device__
var_t	calculate_potential_energy(var_t mu, const vec_t* coor)
{
	return -mu / sqrt((SQR(coor->x) + SQR(coor->y) + SQR(coor->z)));
}

__device__
var_t	calculate_energy(const var_t mu, const vec_t* coor, const vec_t* velo)
{
	return calculate_kinetic_energy(velo) + calculate_potential_energy(mu, coor);
}

__device__
vec_t	cross_product(const vec_t* v, const vec_t* u)
{
	vec_t result;

	result.x = v->y*u->z - v->z*u->y;
    result.y = v->z*u->x - v->x*u->z;
    result.z = v->x*u->y - v->y*u->x;

	return result;
}

__device__
var_t	norm(const vec_t* v)
{
	return sqrt(SQR(v->x) + SQR(v->y) + SQR(v->z));
}

__device__
var_t	norm2(const vec_t* v)
{
	return SQR(v->x) + SQR(v->y) + SQR(v->z);
}

#define	sq3	1.0e-14
__device__
int		calculate_orbelem(const var_t mu, const vec_t* coor, const vec_t* velo, var_t* sma, var_t* ecc)
{
	// Calculate energy, h
    var_t h = calculate_energy(mu, coor, velo);
    if (h >= 0.0) {
        return 1;
    }

    vec_t c = cross_product(coor, velo);
	var_t cNorm2 = norm2(&c);
    /*
    * Calculate eccentricity, e
    */
    var_t e2 = 1.0 + 2.0 * cNorm2 * h / SQR(mu);
	*ecc = abs(e2) < sq3 ? 0.0 : sqrt(e2); 
    /*
    * Calculate semi-major axis, a
    */
    *sma = -mu / (2.0 * h);

    return 0;
}
#undef	sq3

#define	sq2 1.0e-14
#define	sq3	1.0e-14
__device__
	int		calculate_orbelem(const var_t mu, const vec_t* coor, const vec_t* velo, planets::orbelem_t* orbelem)
{
	// Calculate energy, h
    var_t h = calculate_energy(mu, coor, velo);
    if (h >= 0.0) {
        return 1;
    }

    vec_t c = cross_product(coor, velo);
	var_t cNorm2 = norm2(&c);
    /*
    * Calculate eccentricity, e
    */
    var_t e2 = 1.0 + 2.0*cNorm2*h/SQR(mu);
	orbelem->ecc = abs(e2) < sq3 ? 0.0 : sqrt(e2); 
    /*
    * Calculate semi-major axis, a
    */
    orbelem->sma = -mu / (2.0 * h);

    /*
    * Calculate inclination, incl
    */
	cNorm2 = 2(&c);
    var_t cosi = c.z / cNorm;
    var_t sini = sqrt(c.x * c.x + c.y * c.y) / c.Length();
    var_t incl = acos(cosi);
    if (incl < sq2)
    {
        incl = 0.0;
    }
    /*
    * Calculate longitude of node, O
    */
    double node = 0.0;
    if (incl != 0.0)
    {
        double tmpx = -c.y / (c.Length() * sini);
        double tmpy = c.x / (c.Length() * sini);
		node = atan2(tmpy, tmpx);
		ShiftIntoRange(0.0, 2.0*Constants::Pi, node);
    }
    /*
    * Calculate argument of pericenter, w
    */
    double E = 0.0;
    double peri = 0.0;
    if (e2 != 0.0)
    {
        double tmpx = (l.x * cos(node) + l.y * sin(node)) / l.Length();
        double tmpy = (-l.x * sin(node) + l.y * cos(node)) / (l.Length() * cosi);
        peri = atan2(tmpy, tmpx);
        ShiftIntoRange(0.0, 2.0*Constants::Pi, peri);

        tmpx = 1.0 / e * (1.0 - r.Length() / a);
        tmpy = rv / (sqrt(mu * a) * e);
        E = atan2(tmpy, tmpx);
        ShiftIntoRange(0.0, 2.0*Constants::Pi, E);
    }
    else
    {
        peri = 0.0;
        E = atan2(r.y, r.x);
        ShiftIntoRange(0, 2.0*Constants::Pi, E);
    }
    /*
    * Calculate mean anomaly, M
    */
    double M = E - e * sin(E);
    ShiftIntoRange(0, 2.0*Constants::Pi, M);

	orbitalElement->semiMajorAxis			= a;
	orbitalElement->eccentricity			= e;
	orbitalElement->inclination				= incl;
	orbitalElement->argumentOfPericenter	= peri;
	orbitalElement->longitudeOfNode			= node;
	orbitalElement->meanAnomaly				= M;

	return 0;
}
#undef	sq2
#undef	sq3

__global__
void	calculate_grav_accel_kernel(InteractionBound iBound, const planets::param_t* params, const vec_t* coor, vec_t* acce)
{
	int	bodyIdx = iBound.sink.x + blockIdx.x * blockDim.x + threadIdx.x;

	if (bodyIdx < iBound.sink.y) {
		for (int j = iBound.source.x; j < iBound.source.y; j++) 
		{
			if (bodyIdx == j) {
				continue;
			}
			acce[bodyIdx] = calculate_grav_accel_pair(coor[bodyIdx], coor[j], params[j].mass, acce[bodyIdx]);

			//rVec.x = coor[j].x - coor[bodyIdx].x;
			//rVec.y = coor[j].y - coor[bodyIdx].y;
			//rVec.z = coor[j].z - coor[bodyIdx].z;
			//var_t r2 = SQR(rVec.x) + SQR(rVec.y) + SQR(rVec.z);
			//// TODO: how to find out to call the right function in compile time?
			//var_t r = sqrt(r2);
			//var_t invr3 = (var_t)1.0/(r2*r);
			//var_t s = params[j].mass * invr3;

			//acce[bodyIdx].x += s * rVec.x;
			//acce[bodyIdx].y += s * rVec.y;
			//acce[bodyIdx].z += s * rVec.z;
		}
	}
}

__global__
void calculate_drag_accel_kernel(InteractionBound iBound, var_t timeF, const planets::gaspar_t* gaspar, const planets::param_t* params, const vec_t* coor, const vec_t* velo, vec_t* acce)
{
	int	bodyIdx = iBound.sink.x + blockIdx.x * blockDim.x + threadIdx.x;

	if (bodyIdx < iBound.sink.y) {
		var_t r		= sqrt(SQR(coor[bodyIdx].x) + SQR(coor[bodyIdx].y) + SQR(coor[bodyIdx].z));
		vec_t vGas	= gas_velocity(gaspar->eta, K2*params[0].mass, r, atan2(coor[bodyIdx].y, coor[bodyIdx].x));
		var_t rhoGas= gas_density_at(gaspar, r, coor[bodyIdx].z) * timeF;

		vec_t u;
		u.x			= velo[bodyIdx].x -vGas.x;
		u.y			= velo[bodyIdx].y -vGas.y;
		u.z			= velo[bodyIdx].z -vGas.z;

		var_t C		= 0.0;
		// TODO: implement the different regimes according to the mean free path of the gas molecules
		// Epstein-regime:
		{

		}
		// Stokes-regime:
		{
			var_t uLength = sqrt(SQR(vGas.x) + SQR(vGas.y) + SQR(vGas.z));
			C = params[bodyIdx].gamma_stokes * uLength * rhoGas;
		}
		// Transition regime:
		{

		}

		acce[bodyIdx].x = -C * u.x;
		acce[bodyIdx].y = -C * u.y;
		acce[bodyIdx].z = -C * u.z;
	}
}

cudaError_t	planets::call_calculate_grav_accel_kernel(NumberOfBodies nBodies, const planets::param_t* params, const vec_t* coor, vec_t* acce)
{
	cudaError_t cudaStatus = cudaSuccess;
	
	InteractionBound iBound = nBodies.get_self_interacting();

	int		nBodyToCalculate = nBodies.n_self_interacting();
	int		nThread = std::min(THREADS_PER_BLOCK, nBodyToCalculate);
	int		nBlock = (nBodyToCalculate + nThread - 1)/nThread;
	dim3	grid(nBlock);
	dim3	block(nThread);

	calculate_grav_accel_kernel<<<grid, block>>>(iBound, params, coor, acce);
	cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess) {
		std::cerr << "calculate_grav_accel_kernel launch failed: " << cudaGetErrorString(cudaStatus) << std::endl;
		return cudaStatus;
	}

	iBound = nBodies.get_nonself_interacting();
	nBodyToCalculate = nBodies.superPlanetesimal + nBodies.planetesimal;
	nThread		= std::min(THREADS_PER_BLOCK, nBodyToCalculate);
	nBlock		= (nBodyToCalculate + nThread - 1)/nThread;
	grid.x		= nBlock;
	block.x		= nThread;

	calculate_grav_accel_kernel<<<grid, block>>>(iBound, params, coor, acce);
	cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess) {
		std::cerr << "calculate_grav_accel_kernel launch failed: " << cudaGetErrorString(cudaStatus) << std::endl;
		return cudaStatus;
	}

	iBound = nBodies.get_non_interacting();
	nBodyToCalculate = nBodies.testParticle;
	nThread		= std::min(THREADS_PER_BLOCK, nBodyToCalculate);
	nBlock		= (nBodyToCalculate + nThread - 1)/nThread;
	grid.x		= nBlock;
	block.x		= nThread;

	calculate_grav_accel_kernel<<<grid, block>>>(iBound, params, coor, acce);
	cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess) {
		std::cerr << "calculate_grav_accel_kernel launch failed: " << cudaGetErrorString(cudaStatus) << std::endl;
		return cudaStatus;
	}

	return cudaStatus;
}

cudaError_t planets::call_calculate_drag_accel_kernel(NumberOfBodies nBodies, ttt_t time, const planets::gaspar_t* gaspar, const planets::param_t* params, const vec_t* coor, const vec_t* velo, vec_t* acce)
{
	cudaError_t cudaStatus = cudaSuccess;

	// TODO: calculate it using the value of the time
	var_t timeF = 1.0;
	
	InteractionBound iBound = nBodies.get_bodies_gasdrag();

	int		nBodyToCalculate = nBodies.superPlanetesimal + nBodies.planetesimal;
	int		nThread = std::min(THREADS_PER_BLOCK, nBodyToCalculate);
	int		nBlock = (nBodyToCalculate + nThread - 1)/nThread;
	dim3	grid(nBlock);
	dim3	block(nThread);

	calculate_drag_accel_kernel<<<grid, block>>>(iBound, timeF, gaspar, params, coor, velo, acce);
	cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess) {
		std::cerr << "calculate_drag_accel_kernel launch failed: " << cudaGetErrorString(cudaStatus) << std::endl;
		return cudaStatus;
	}

	return cudaStatus;
}

