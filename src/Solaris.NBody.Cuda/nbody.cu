#include "nbody.h"

#include <cuda_runtime.h>
#include "device_launch_parameters.h"

#include <iostream>
#include <fstream>

#include "config.h"
#include "nbody_exception.h"
#include "util.h"

using namespace std;

// Calculate acceleration caused by one particle on another
__inline__ __device__ 
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
__global__
	void calculate_accel_kernel(const nbody::param_t* p, const vec_t* c, vec_t* a)
{
	// Index of this particle
	int i = blockDim.x * blockIdx.x + threadIdx.x;

	// Start index of other particle
	int j = NTILE * blockIdx.y;

	// Copy data into shared memory
	//__shared__ vec_t s_c1[NTILE];
	__shared__ vec_t s_c2[NTILE];
	__shared__ nbody::param_t s_p[NTILE];

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

__global__
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

// Determines if two intervals are disjoint
__inline__ __device__ 
	bool detect_disjoint_interval(var_t a1, var_t a2, var_t b1, var_t b2, var_t buffer)
{
	return
		min(b1, b2) - max(a1, a2) > buffer |
		min(a1, a2) - max(b1, b2) > buffer;
}

// Determines if bounding boxes of trajectory elements are disjoint
__inline__ __device__ 
	bool detect_disjoint_pair(var_t r1, var_t r2, vec_t c1in, vec_t c1out, vec_t c2in, vec_t c2out, var_t buffer)
{
	return
		detect_disjoint_interval(c1in.x, c1out.x, c2in.x, c2out.x, (r1 + r2) * buffer) |
		detect_disjoint_interval(c1in.y, c1out.y, c2in.y, c2out.y, (r1 + r2) * buffer) |
		detect_disjoint_interval(c1in.z, c1out.z, c2in.z, c2out.z, (r1 + r2) * buffer);
}

// Detect potentially intersecting trajectories
// Ezt itt lehetne koordinátánként is csinálni!!!
__global__ 
	void detect_intersections_kernel(const nbody::param_t* p, const vec_t* cin, const vec_t* cout, var_t buffer, int2_t *interactions, int* interactions_end)
{
	// Index of this particle
	int_t i = blockDim.x * blockIdx.x + threadIdx.x;

	// Start index of other particle
	int_t j = NTILE * blockIdx.y;

	// Copy data into shared memory
	__shared__ nbody::param_t s_p1[NTILE];
	__shared__ nbody::param_t s_p2[NTILE];
	__shared__ vec_t s_c1in[NTILE];
	__shared__ vec_t s_c1out[NTILE];
	__shared__ vec_t s_c2in[NTILE];
	__shared__ vec_t s_c2out[NTILE];
	
	bool_t rr[NTILE];

	s_p1[threadIdx.x] = p[i];
	s_p2[threadIdx.x] = p[j + threadIdx.x];
	s_c1in[threadIdx.x] = cin[i];
	s_c1out[threadIdx.x] = cout[i];
	s_c2in[threadIdx.x] = cin[j + threadIdx.x];
	s_c2out[threadIdx.x] = cout[j + threadIdx.x];

	__syncthreads();

	// Loop over other particles in the tile to detect intersections
	// of bounding boxes along the trajectory

	int count = 0;

#pragma unroll NTILE
	for (int_t k = 0; k < NTILE; k ++)
	{
		bool_t intersecting =
			j + k > i &
			s_p1[threadIdx.x].mass != 0 & s_p2[k].mass != 0 &
			!detect_disjoint_pair(
				s_p1[threadIdx.x].radius, s_p2[k].radius,
				s_c1in[threadIdx.x], s_c1out[threadIdx.x],
				s_c2in[k], s_c2out[k],
				buffer);

		rr[k] = intersecting;

		count += intersecting ? 1 : 0;
	}

	if (count == 0)
	{
		return;
	}

	int end = atomicAdd(interactions_end, count) + count;

#pragma unroll NTILE
	for (int_t k = 0; k < NTILE; k ++)
	{
		if (rr[k])
		{
			interactions[end - count].x = i;
			interactions[end - count].y = j + k;
			count--;
		}
	}
}

// Detect colliding particle pair
__inline__ __device__ 
	bool detect_collision_pair(nbody::param_t p1, nbody::param_t p2, vec_t c1, vec_t c2, vec_t v1, vec_t v2, var_t buffer)
{
	// Assume all intersecting pairs colliding
	vec_t d;
	
	d.x = c1.x - c2.x;
	d.y = c1.y - c2.y;
	d.z = c1.z - c2.z;

	d.w = d.x * d.x + d.y * d.y + d.z * d.z;

	var_t r = (p1.radius + p2.radius) * buffer;
	r *= r;

	return d.w < r;
}

__inline__ __device__ 
	void update_colliding_pair(nbody::param_t* p,	vec_t* c, vec_t* v,	int2_t ii)
{
	var_t M = p[ii.x].mass + p[ii.y].mass;
		
	c[ii.x].x = (p[ii.x].mass * c[ii.x].x + p[ii.y].mass * c[ii.y].x) / M;
	c[ii.x].y = (p[ii.x].mass * c[ii.x].y + p[ii.y].mass * c[ii.y].y) / M;
	c[ii.x].z = (p[ii.x].mass * c[ii.x].z + p[ii.y].mass * c[ii.y].z) / M;

	v[ii.x].x = (p[ii.x].mass * v[ii.x].x + p[ii.y].mass * v[ii.y].x) / M;
	v[ii.x].y = (p[ii.x].mass * v[ii.x].y + p[ii.y].mass * v[ii.y].y) / M;
	v[ii.x].z = (p[ii.x].mass * v[ii.x].z + p[ii.y].mass * v[ii.y].z) / M;

	p[ii.x].mass = M;
	p[ii.y].mass = 0;

	var_t V =
		p[ii.x].radius * p[ii.x].radius * p[ii.x].radius +
		p[ii.y].radius * p[ii.y].radius * p[ii.y].radius;
	p[ii.x].radius = pow(V, 1.0/3.0);
}

// Detect colliding particles and handle collisins
// by sticking particles together
__global__ 
	void detect_collisions_kernel(nbody::param_t* p,	vec_t* c, vec_t* v,	ttt_t time,	var_t buffer, int2_t *interactions, int* interactions_end, nbody::collision_t *collisions, int* collisions_end)
{

	// Each thread will check the collision of one potentially colliding
	// pair of particles in the interactions list

	int_t i = blockDim.x * blockIdx.x + threadIdx.x;

	// Allocate shared storege for collision detections
	__shared__ bool s_collisions[NTILE];

	// Check collision
	int2_t ii;
	bool_t colliding;
	if (i < *interactions_end)
	{
		ii = interactions[i];
		colliding = detect_collision_pair(
			p[ii.x], p[ii.y], 
			c[ii.x], c[ii.y], 
			v[ii.x], v[ii.y], 
			buffer);
	}
	else
	{
		colliding = false;
	}

	s_collisions[threadIdx.x] = colliding;

	// Count colliding particles in the block
	int count = __syncthreads_count(colliding);

	// Expand collisions list to accommodate new collision data
	int end;
	if (threadIdx.x == 0)
	{
		end = atomicAdd(collisions_end, count) + count;
	}

	// Update colliding particles' parameters
	if (i < *interactions_end && colliding)
	{
		update_colliding_pair(p, c, v, ii);
	}

	__syncthreads();

	// Copy collision information into the list
	if (threadIdx.x == 0)
	{
		for (int k = 0; k < NTILE; k ++)
		{
			if (s_collisions[k])
			{
				collisions[end - count].time = time;
				collisions[end - count].x = ii.x;
				collisions[end - count].y = ii.y;

				count--;
			}
		}
	}
}


nbody::nbody(int n, ttt_t t0) :
	ode(2, t0),
	n(n),
	d_accelerations(d_vec_t()),
	d_interactions(d_int2_t()),
	d_interactions_end(d_int_t()),
	h_collisions(h_collision_t()),
	h_collisions_end(h_int_t()),
	d_collisions(d_collision_t()),
	d_collisions_end(d_int_t())
{
	round_up_n();
	allocate_vectors();
}

nbody::~nbody()
{
}

void nbody::allocate_vectors()
{
	// Allocate vector for acceleration intermediate results
	d_accelerations.resize(n * n / NTILE);	
	
	h_interactions_end.resize(1);
	d_interactions.resize(n * n);
	d_interactions_end.resize(1);

	h_collisions.resize(n * n / 2);
	h_collisions_end.resize(1);
	d_collisions.resize(n * n / 2);
	d_collisions_end.resize(1);

	// Parameters
	h_p.resize(NPAR * n);

	// Aliases to coordinates and velocities
	h_y[0].resize(NDIM * n);
	h_y[1].resize(NDIM * n);
}

void nbody::round_up_n()
{
	// Round up n to the number of bodies per tile
	int m = ((n + NTILE - 1) / NTILE) * NTILE;

	if (n != m) {
		cerr << "Number of bodies rounded up to " << m << endl;
	}
	n = m;
}

void nbody::call_calculate_accel_kernel(const param_t* p, const vec_t* c, vec_t* atemp)
{
	dim3 blocks(n / NTILE, n / NTILE);
	dim3 threads(NTILE);

	calculate_accel_kernel<<<blocks, threads>>>(p, c, atemp);
}

void nbody::call_sum_accel_kernel(const vec_t* atemp, vec_t* a)
{
	dim3 blocks(n / NTILE);
	dim3 threads(NTILE);

	sum_accel_kernel<<<blocks, threads>>>(atemp, a);
}

int nbody::call_detect_intersections_kernel(d_var_t& cin, d_var_t& cout)
{
	h_interactions_end[0] = 0;
	thrust::copy(h_interactions_end.begin(), h_interactions_end.end(), d_interactions_end.begin());

	dim3 blocks(n / NTILE, n / NTILE);
	dim3 threads(NTILE);

	detect_intersections_kernel<<<blocks, threads>>>
		(
			(param_t*)thrust::raw_pointer_cast(d_p.data()),
			(vec_t*)thrust::raw_pointer_cast(cin.data()),
			(vec_t*)thrust::raw_pointer_cast(cout.data()),
			buffer_radius,
			(int2_t*)thrust::raw_pointer_cast(d_interactions.data()),
			(int_t*)thrust::raw_pointer_cast(d_interactions_end.data())
		);

	thrust::copy(d_interactions_end.begin(), d_interactions_end.end(), h_interactions_end.begin());
	
	return h_interactions_end[0];
}

int nbody::call_detect_collisions_kernel()
{
	dim3 blocks((h_interactions_end[0] + (NTILE - 1)) / NTILE);
	dim3 threads(NTILE);

	detect_collisions_kernel <<< blocks, threads >>> (
		(param_t*)thrust::raw_pointer_cast(d_p.data()),
		(vec_t*)thrust::raw_pointer_cast(d_y[0].data()),
		(vec_t*)thrust::raw_pointer_cast(d_y[1].data()),
		t,
		buffer_radius,
		(int2_t*)thrust::raw_pointer_cast(d_interactions.data()),
		(int*)thrust::raw_pointer_cast(d_interactions_end.data()),
		(collision_t*)thrust::raw_pointer_cast(d_collisions.data()),
		(int*)thrust::raw_pointer_cast(d_collisions_end.data()));

	thrust::copy(d_collisions_end.begin(), d_collisions_end.end(), h_collisions_end.begin());
	
	// Ez még lehet, hogy máshova kerül majd
	if (h_collisions_end[0] > 0)
	{
		thrust::copy(d_collisions.begin(), d_collisions.begin() + h_collisions_end[0], h_collisions.begin());
	}

	return h_collisions_end[0];
}

void nbody::calculate_dy(int i, int r, ttt_t t, const d_var_t& p, const std::vector<d_var_t>& y, d_var_t& dy)
{
	switch (i)
	{
	case 0:
		// Copy velocities from previous step
		thrust::copy(y[1].begin(), y[1].end(), dy.begin());
		break;
	case 1:
		// Calculate accelerations in tiles

		call_calculate_accel_kernel((param_t*)p.data().get(), (vec_t*)y[0].data().get(), (vec_t*)d_accelerations.data().get());

		// Now d_atemp contains the accelerations
		// Sum up results from tiles into the output vector
		call_sum_accel_kernel(d_accelerations.data().get(), (vec_t*)dy.data().get());

		break;
	}
}

// itt valahogy vissza kellene kapni a d_ bla-bla-t!
int nbody::detect_collisions()
{
	int ic = call_detect_intersections_kernel(d_y[0], d_yout[0]);

	if (ic > 0)
	{
		 int cc = call_detect_collisions_kernel();
		 return cc;
	}
	else
	{
		return 0;
	}
}

void nbody::load(string filename, int n)
{
	if (n > this->n) {
		throw nbody_exception("Too many lines in file.");
	}

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
int nbody::print_positions(ostream& sout)
{
	param_t* h_param = (param_t*)h_p.data();
	vec_t* h_coord = (vec_t*)h_y[0].data();
	vec_t* h_veloc = (vec_t*)h_y[1].data();
	
	for (int i = 0; i < n; i ++)
	{
		if (h_param[i].mass == 0)
		{
			continue;
		}

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

	return n;
}

int nbody::print_collisions(ostream& sout, int start)
{
	for (int i = start; i < h_collisions_end[0]; i++)
	{
		sout << h_collisions[i].time << '\t';
		sout << h_collisions[i].x << '\t';
		sout << h_collisions[i].y;
		sout << endl;
	}

	return h_collisions_end[0];
}

#ifdef FALSE
// Generate random positions and velocities
void nbody::generate_random(int n)
{
	if (n > this->n)
	{
		throw nbody_exception("Too many objects to generate.");
	}

	vec_t* h_coord = (vec_t*)h_y[0].data();
	vec_t* h_veloc = (vec_t*)h_y[1].data();
	param_t* h_param = (param_t*)h_p.data();

    for (int i = 0; i < n; i++)
    { 
		h_coord[i].x = randf(-5, 5);
		h_coord[i].y = randf(-5, 5);
		h_coord[i].z = randf(-5, 5);

		h_veloc[i].x = randf((var_t)-0.001, (var_t)0.001);
		h_veloc[i].y = randf((var_t)-0.001, (var_t)0.001);
		h_veloc[i].z = randf((var_t)-0.001, (var_t)0.001);

		h_param[i].mass = 1;
    }

	// Fill up rest with zeros
	for (int i = n; i < this->n; i++)
    { 
		h_coord[i].x = 0;
		h_coord[i].y = 0;
		h_coord[i].z = 0;

		h_veloc[i].x = 0;
		h_veloc[i].y = 0;
		h_veloc[i].z = 0;

		h_param[i].mass = 0;
    }
}

// Generate some objects that are about to collide
void nbody::generate_colliding()
{
	vec_t* h_coord = (vec_t*)h_y[0].data();
	vec_t* h_veloc = (vec_t*)h_y[1].data();
	param_t* h_param = (param_t*)h_p.data();

	h_coord[0].x = 1;
	h_coord[0].y = 0.001;
	h_coord[0].z = 0;

	h_veloc[0].x = -0.001;
	h_veloc[0].y = 0;
	h_veloc[0].z = 0;

	h_param[0].mass = 0.01;
	h_param[0].radius = 0.01;

	h_coord[1].x = -1;
	h_coord[1].y = -0.001;
	h_coord[1].z = 0;

	h_veloc[1].x = 0.001;
	h_veloc[1].y = 0;
	h_veloc[1].z = 0;

	h_param[1].mass = 0.01;
	h_param[1].radius = 0.01;

	// Fill up rest with zeros
	for (int i = 2; i < this->n; i++)
    { 
		h_coord[i].x = 0;
		h_coord[i].y = 0;
		h_coord[i].z = 0;

		h_veloc[i].x = 0;
		h_veloc[i].y = 0;
		h_veloc[i].z = 0;

		h_param[i].mass = 0;
    }
}
#endif