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

// Calculate acceleration caused by one particle on another
static __device__ 
vec_t calculate_accel_pair(const vec_t c1, const vec_t c2, var_t m, vec_t a)
{
	vec_t d;
	
	d.x = c1.x - c2.x;
	d.y = c1.y - c2.y;
	d.z = c1.z - c2.z;

	d.w = d.x * d.x + d.y * d.y + d.z * d.z;
	if (d.w > 0)
	{
		d.w = d.w * d.w * d.w;
		d.w = - K2 * m / sqrt(d.w);

		a.x += d.x * d.w;
		a.y += d.y * d.w;
		a.z += d.z * d.w;
	}

	return a;
}

// Calculate and sum up accelerations
static __global__
void calculate_accel_kernel(const pp_disk::param_t* p, const vec_t* c, vec_t* a)
{
	// Index of this particle
	int i = blockDim.x * blockIdx.x + threadIdx.x;

	// Start index of other particle
	int j = NTILE * blockIdx.y;

	// Copy data into shared memory
	//__shared__ vec_t s_c1[NTILE];
	__shared__ vec_t s_c2[NTILE];
	__shared__ pp_disk::param_t s_p[NTILE];

	//s_c1[threadIdx.x] = c[i];
	s_c2[threadIdx.x] = c[j + threadIdx.x];
	s_p[threadIdx.x] = p[j + threadIdx.x];

	// Loop over other particles in the tile and sum up acceleration
	vec_t aa = { 0, 0, 0, 0};

	__syncthreads();

//#pragma unroll NTILE
	for (int k = 0; k < NTILE; k ++)
	{
		//aa = calculate_accel_pair(s_c1[threadIdx.x], s_c2[k], s_p[k].mass, aa);
		aa = calculate_accel_pair(c[i], s_c2[k], s_p[k].mass, aa);
	}

	a[gridDim.x * blockDim.x * blockIdx.y + i] = aa;
}

static __global__
void sum_accel_kernel(const vec_t* a, vec_t* suma)
{
	int i = blockDim.x * blockIdx.x + threadIdx.x;

	// Loop over blocks and sum up acceleration
	vec_t aa = { 0, 0, 0, 0};
	vec_t ab;

#pragma unroll
	for (int k = 0; k < gridDim.x; k++)
	{
		//ab = a[NBODY * k + i];
		ab = a[gridDim.x * blockDim.x * k + i];
		
		aa.x += ab.x;
		aa.y += ab.y;
		aa.z += ab.z;
	}

	suma[i] = aa;
}

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
		// Calculate accelerations in tiles

		call_calculate_accel_kernel((param_t*)p.data().get(), (vec_t*)y[0].data().get(), (vec_t*)dy.data().get());

		// Now d_atemp contains the accelerations
		// Sum up results from tiles into the output vector
		//call_sum_accel_kernel(d_accelerations.data().get(), (vec_t*)dy.data().get());

		break;
	}
}

void pp_disk::load(string filename, int n)
{
	vec_t* h_coord = (vec_t*)h_y[0].data();
	vec_t* h_veloc = (vec_t*)h_y[1].data();
	param_t* h_param = (param_t*)h_p.data();

	ifstream input(filename.c_str());

	if (input) {
        int		id;
		ttt_t	time;
        
        for (int i = 0; i < n; i++) { 
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
