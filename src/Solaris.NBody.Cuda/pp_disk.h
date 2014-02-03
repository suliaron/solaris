#pragma once

#include <string>

#include "ode.h"
#include "config.h"

class number_of_bodies;

using namespace std;

class pp_disk : public ode
{
public:
	// Type for parameters
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
	} param_t;

	typedef thrust::host_vector<param_t>		h_param_t;
	typedef thrust::device_vector<param_t>		d_param_t;
	
	pp_disk(number_of_bodies *nBodies);
	~pp_disk();

	void calculate_dy(int i, int r, ttt_t t, const d_var_t& p, const std::vector<d_var_t>& y, d_var_t& dy);

	void load(string filename, int n);
	int print_positions(ostream& sout);

private:
	number_of_bodies	*nBodies;

	void allocate_vectors();

	cudaError_t call_calculate_accel_kernel(const param_t* params, const vec_t* coor, vec_t* acce);

};