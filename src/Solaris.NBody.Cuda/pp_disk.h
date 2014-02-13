#pragma once

#include <string>

#include "ode.h"
#include "config.h"

class number_of_bodies;
class gas_disc;

using namespace std;

class pp_disk : public ode
{
public:
	typedef enum migration_type
	{
		NO,
		TYPE_I,
		TYPE_II
	} migration_type_t;

	// Type for parameters
	typedef struct param
	{
		//! Unique number to identify the object
		int_t	id;
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
		//! Type of the migration
		migration_type_t migType;
		//! The migration stop at this distance measured from the star
		var_t	migStopAt;
	} param_t;

	typedef struct orbelem
	{
		//! Semimajor-axis of the body
		var_t sma;
		//! Eccentricity of the body
		var_t ecc;
		//! Inclination of the body
		var_t inc;
		//! Argument of the pericenter
		var_t peri;
		//! Longitude of the ascending node
		var_t node;
		//! Mean anomaly
		var_t mean;
	} orbelem_t;

	typedef thrust::host_vector<param_t>		h_param_t;
	typedef thrust::device_vector<param_t>		d_param_t;

	typedef thrust::host_vector<orbelem_t>		h_orbelem_t;
	typedef thrust::device_vector<orbelem_t>	d_orbelem_t;

	d_orbelem_t			d_orbelem;
	h_orbelem_t			h_orbelem;
	
	pp_disk(number_of_bodies *nBodies, gas_disc *gasDisc);
	~pp_disk();

	void calculate_orbelem(int_t refBodyId);

	void calculate_dy(int i, int r, ttt_t t, const d_var_t& p, const std::vector<d_var_t>& y, d_var_t& dy);

	void load(string filename, int n);
	int print_positions(ostream& sout);
	int print_orbelem(ostream& sout);

private:
	number_of_bodies	*nBodies;
	gas_disc			*gasDisc;
	gas_disc			*d_gasDisc;

	d_var_t				acceGasDrag;
	d_var_t				acceMigrateI;
	d_var_t				acceMigrateII;

	void allocate_vectors();

	//! Calls the kernel that calculates the accelerations from gravitational
	/*  interactions.
		\param params Vector of parameters of the bodies
		\param coor Vector of coordinates of the bodies
		\param acce Will hold the accelerations for each body
	*/
	cudaError_t call_calculate_grav_accel_kernel(const param_t* params, const vec_t* coor, vec_t* acce);

	//! Calls the kernel that calculates the acceleration due to drag force.
	/*
		\param time The actual time of the simulation
		\param gasDisc The parameters describing the gas disk
		\param params Vector of parameters of the bodies
		\param coor Vector of coordinates of the bodies
		\param velo Vector of velocities of the bodies
		\param acce Will hold the accelerations for each body
	*/
	cudaError_t call_calculate_drag_accel_kernel(ttt_t time, const gas_disc* gasDisc, const param_t* params, const vec_t* coor, const vec_t* velo, vec_t* acce);

};

static __host__ __device__ void		shift_into_range(var_t lower, var_t upper, var_t* value);
static __host__ __device__ vec_t	cross_product(const vec_t* v, const vec_t* u);
static __host__ __device__ var_t	dot_product(const vec_t* v, const vec_t* u);
static __host__ __device__ var_t	norm2(const vec_t* v);
static __host__ __device__ var_t	norm(const vec_t* v);
static __host__ __device__ vec_t	circular_velocity(var_t mu, const vec_t* rVec);
static __host__ __device__ vec_t	gas_velocity(var2_t eta, var_t mu, const vec_t* rVec);
static __host__ __device__ var_t	gas_density_at(const gas_disc* gasDisc, const vec_t* rVec);
static __host__ __device__ var_t	calculate_kinetic_energy(const vec_t* vVec);
static __host__ __device__ var_t	calculate_potential_energy(var_t mu, const vec_t* rVec);
static __host__ __device__ var_t	calculate_energy(var_t mu, const vec_t* rVec, const vec_t* vVec);
static __host__ __device__ int_t	kepler_equation_solver(var_t ecc, var_t mean, var_t eps, var_t* E);
static __host__ __device__ int_t	calculate_phase(var_t mu, const pp_disk::orbelem_t* oe, vec_t* rVec, vec_t* vVec);
static __host__ __device__ int_t	calculate_sma_ecc(var_t mu, const vec_t* coor, const vec_t* velo, var_t* sma, var_t* ecc);
static __host__ __device__ int_t	calculate_orbelem(var_t mu, const vec_t* coor, const vec_t* velo, pp_disk::orbelem_t* orbelem);
static __host__ __device__ var_t	orbital_period(var_t mu, var_t sma);
static __host__ __device__ var_t	orbital_frequency(var_t mu, var_t sma);
static __host__ __device__ var_t	calculate_gamma_stokes(var_t cd, var_t density, var_t radius);
static __host__ __device__ var_t	calculate_gamma_epstein(var_t density, var_t radius);
