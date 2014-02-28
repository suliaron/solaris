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

	typedef thrust::host_vector<param_t>		h_param_t;
	typedef thrust::device_vector<param_t>		d_param_t;
	
	typedef thrust::host_vector<collision_t>	h_collision_t;
	typedef thrust::device_vector<collision_t>	d_collision_t;

	var_t buffer_radius;

	nbody(int n, ttt_t t0);
	~nbody();

	void calculate_dy(int i, int r, ttt_t t, const d_var_t& p, const std::vector<d_var_t>& y, d_var_t& dy);

	int detect_collisions();

	void load(string filename, int n);
	int print_positions(ostream& sout);
	int print_collisions(ostream& sout, int start);

private:
	// Number of bodies
	int				n;
	// Temporary storage for accelerations, required between kernel calls
	d_vec_t			d_accelerations;
	
	// Storage for intersection detection
	h_int_t			h_interactions_end;
	d_int2_t		d_interactions;
	d_int_t			d_interactions_end;
	
	// Storage for collisions
	h_collision_t	h_collisions;
	h_int_t			h_collisions_end;
	d_collision_t	d_collisions;
	d_int_t			d_collisions_end;

	void round_up_n();
	void allocate_vectors();

	void call_calculate_accel_kernel(const param_t* p, const vec_t* c, vec_t* atemp);
	void call_sum_accel_kernel(const vec_t* atemp, vec_t* a);
	int call_detect_intersections_kernel(d_var_t& cin, d_var_t& cout);
	int call_detect_collisions_kernel();

};