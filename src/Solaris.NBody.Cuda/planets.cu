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

// includes project's header
#include "config.h"
#include "interaction_bound.h"
#include "number_of_bodies.h"
#include "nbody_exception.h"
#include "planets.h"
#include "nbody_util.h"

#define THREADS_PER_BLOCK	256


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
void calculate_drag_accel_kernel(interaction_bound iBound, var_t timeF, const planets::gaspar_t* gaspar, const planets::param_t* params, const vec_t* coor, const vec_t* velo, vec_t* acce)
{
	int	bodyIdx = iBound.sink.x + blockIdx.x * blockDim.x + threadIdx.x;

	if (bodyIdx < iBound.sink.y) {
		var_t r		= sqrt(SQR(coor[bodyIdx].x) + SQR(coor[bodyIdx].y) + SQR(coor[bodyIdx].z));
		vec_t vGas	= gas_velocity(gaspar->eta, K2*params[0].mass, r, atan2(coor[bodyIdx].y, coor[bodyIdx].x));
		var_t rhoGas= gas_density_at(gaspar, r, coor[bodyIdx].z) * timeF;

		vec_t u;
		u.x			= velo[bodyIdx].x - vGas.x;
		u.y			= velo[bodyIdx].y - vGas.y;
		u.z			= velo[bodyIdx].z - vGas.z;

		var_t C		= 0.0;
		// TODO: implement the different regimes according to the mean free path of the gas molecules
		// Epstein-regime:
		{

		}
		// Stokes-regime:
		{
			var_t uLength = norm(&u);
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

cudaError_t	planets::call_calculate_grav_accel_kernel(number_of_bodies nBodies, const planets::param_t* params, const vec_t* coor, vec_t* acce)
{
	cudaError_t cudaStatus = cudaSuccess;
	
	interaction_bound iBound = nBodies.get_self_interacting();

	int		nBodyToCalculate = nBodies.n_self_interacting();
	int		nThread = std::min(THREADS_PER_BLOCK, nBodyToCalculate);
	int		nBlock = (nBodyToCalculate + nThread - 1)/nThread;
	dim3	grid(nBlock);
	dim3	block(nThread);

	calculate_grav_accel_kernel<<<grid, block>>>(iBound, params, coor, acce);
	if ((cudaStatus = cudaGetLastError()) != cudaSuccess) {
		throw nbody_exception("calculate_grav_accel_kernel launch failed", cudaStatus);
	}

	iBound = nBodies.get_nonself_interacting();
	nBodyToCalculate = nBodies.super_planetesimal + nBodies.planetesimal;
	if (nBodyToCalculate > 0) {
		nThread		= std::min(THREADS_PER_BLOCK, nBodyToCalculate);
		nBlock		= (nBodyToCalculate + nThread - 1)/nThread;
		grid.x		= nBlock;
		block.x		= nThread;

		calculate_grav_accel_kernel<<<grid, block>>>(iBound, params, coor, acce);
		cudaStatus = cudaGetLastError();
		if (cudaStatus != cudaSuccess) {
			throw nbody_exception("calculate_grav_accel_kernel launch failed", cudaStatus);
		}
	}

	iBound = nBodies.get_non_interacting();
	nBodyToCalculate = nBodies.test_particle;
	if (nBodyToCalculate > 0) {
		nThread		= std::min(THREADS_PER_BLOCK, nBodyToCalculate);
		nBlock		= (nBodyToCalculate + nThread - 1)/nThread;
		grid.x		= nBlock;
		block.x		= nThread;

		calculate_grav_accel_kernel<<<grid, block>>>(iBound, params, coor, acce);
		cudaStatus = cudaGetLastError();
		if (cudaStatus != cudaSuccess) {
			throw nbody_exception("calculate_grav_accel_kernel launch failed", cudaStatus);
		}
	}

	return cudaStatus;
}

cudaError_t planets::call_calculate_drag_accel_kernel(number_of_bodies nBodies, ttt_t time, const planets::gaspar_t* gaspar, const planets::param_t* params, const vec_t* coor, const vec_t* velo, vec_t* acce)
{
	cudaError_t cudaStatus = cudaSuccess;

	// TODO: calculate it using the value of the time
	var_t timeF = 1.0;
	
	interaction_bound iBound = nBodies.get_bodies_gasdrag();

	int		nBodyToCalculate = nBodies.super_planetesimal + nBodies.planetesimal;
	int		nThread = std::min(THREADS_PER_BLOCK, nBodyToCalculate);
	int		nBlock = (nBodyToCalculate + nThread - 1)/nThread;
	dim3	grid(nBlock);
	dim3	block(nThread);

	calculate_drag_accel_kernel<<<grid, block>>>(iBound, timeF, gaspar, params, coor, velo, acce);
	cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess) {
		throw nbody_exception("calculate_drag_accel_kernel launch failed", cudaStatus);
	}

	return cudaStatus;
}

planets::planets(number_of_bodies bodies) :
	ode(2),
	bodies(bodies)
{
	//round_up_n();
	allocate_vectors();
}

planets::~planets()
{
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
	switch (i)
	{
	case 0:
		// Copy velocities from previous step
		thrust::copy(y[1].begin(), y[1].end(), dy.begin());
		break;
	case 1:
		// Calculate accelerations originated from gravity
		call_calculate_grav_accel_kernel(bodies, (param_t*)p.data().get(), (vec_t*)d_y[0].data().get(), (vec_t*)dy.data().get());

		// Calculate accelerations originated from gas drag

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
		if (h_param[i].mass == 0)
		{
			continue;
		}

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
