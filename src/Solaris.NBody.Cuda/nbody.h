#pragma once

#include <string>

#include "ode.h"
#include "config.h"

using namespace std;

class nbody : public ode
{
public:
	// Type for parameters
	typedef struct param
			{
				var_t mass;		// particle mass
				var_t radius;	// particle radius
			} param_t;

	typedef struct collision
			{
				ttt_t time;
				int_t x;
				int_t y;
			} collision_t;

	typedef thrust::host_vector<param_t> host_param_t;
	typedef thrust::device_vector<param_t> device_param_t;
	
	typedef thrust::host_vector<collision_t> host_collision_t;
	typedef thrust::device_vector<collision_t> device_collision_t;

	var_t buffer_radius;

private:
	// Number of bodies
	int n;
	// Temporary storage for accelerations, required between kernel calls
	d_vec_t d_accelerations;
	
	// Dtorage for intersection detection
	host_int_t h_interactions_end;
	device_int2_t d_interactions;
	device_int_t d_interactions_end;
	
	// Storage for collisions
	host_collision_t h_collisions;
	host_int_t h_collisions_end;
	device_collision_t d_collisions;
	device_int_t d_collisions_end;

public:
	nbody(int n);
	~nbody();

private:
	void round_up_n();
	void allocate_vectors();

	void call_calculate_accel_kernel(const param_t* p, const vec_t* c, vec_t* atemp);
	void call_sum_accel_kernel(const vec_t* atemp, vec_t* a);
	int call_detect_intersections_kernel(d_var_t& cin, d_var_t& cout);
	int call_detect_collisions_kernel();

public:
	void calculate_dy(int i, int r, ttt_t t, const d_var_t& p, const std::vector<d_var_t>& y, d_var_t& dy);
	int detect_collisions();

	void load(string filename, int n);
	int print_positions(ostream& sout);
	int print_collisions(ostream& sout, int start);
};