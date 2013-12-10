#pragma once

#include <string>

#include "ode.h"
#include "config.h"

using namespace std;

class planets : public ode
{
public:
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

	typedef thrust::host_vector<param_t>	h_param_t;
	typedef thrust::device_vector<param_t>	d_param_t;

	typedef struct gaspar
	{
		var_t c_rho;
		var_t p_rho;
		var_t c_h;
		var_t p_h;
		var_t c_eta;
		var_t p_eta;
		var_t c_tau;
		var_t p_tau;
	} gaspar_t;

	typedef enum migration
	{
		type1,
		type2
	} migration_t;

private:

	//! Calls the kernel that calculates the accelerations from gravitational
	/*  interactions. Done in tiles.
		\param p Vector of parameters of the bodies
		\param c Vector of coordinates of the bodies
		\param bounds Vector of indices limiting the interacting pairs
		\param atemp Will hold the accelerations for each body per each tile
	*/
	cudaError_t	call_calculate_grav_accel_kernel(NumberOfBodies nBodies, const planets::param_t* p, const vec_t* c, vec_t* atemp);
	// void call_calculate_grav_accel_kernel(const param_t* p, const vec_t* c, const int4_t bounds, vec_t* atemp);

	//! Calls the kernel that sums up accelerations by tiles
	/*
		\param atemp Vector of accelerations for each body for each tile
		\param a Will hold the summed acceleration for each body
	*/
	void call_sum_grav_accel_kernel(const vec_t* atemp, vec_t* a);

	//! Calls the kernel that calculates the acceleration due to drag force on bodies
	cudaError_t call_calculate_drag_accel_kernel(NumberOfBodies nBodies, ttt_t time, const planets::gaspar_t& gaspar, const planets::param_t* params, const vec_t* coor, const vec_t* velo, vec_t* acce);

	cudaError_t call_calculate_epheremis_kernel(const param_t* p, const vec_t* c, const vec_t* v, int2_t bounds);

	cudaError_t call_calculate_migration_accel_kernel(const param_t* p, const vec_t* c, const vec_t* v, int2_t bounds, gaspar_t gaspar, vec_t* atemp);

};