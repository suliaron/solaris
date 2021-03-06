#pragma once

#include <string>

#include "config.h"
#include "interaction_bound.h"
#include "ode.h"

class number_of_bodies;
class gas_disk;

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
		int_t id;
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
	
	pp_disk(number_of_bodies *nBodies, gas_disk *gasDisk, ttt_t t0);
	~pp_disk();

	h_orbelem_t calculate_orbelem(int_t refBodyId);

	void calculate_dy(int i, int r, ttt_t t, const d_var_t& p, const std::vector<d_var_t>& y, d_var_t& dy);

	void load(string filename, int n);
	int print_positions(ostream& sout);
	int print_orbelem(ostream& sout);

	//! Calls the cpu function that calculates the accelerations from gravitational
	/*  interactions.
		\param iBound containes the staring and ending indicies of the sink and source bodies
		\param params Vector of parameters of the bodies
		\param coor Vector of coordinates of the bodies
		\param acce Will hold the accelerations for each body
	*/
	void calculate_grav_accel(interaction_bound iBound, const param_t* params, const vec_t* coor, vec_t* acce);

private:
	dim3	grid;
	dim3	block;

	number_of_bodies	*nBodies;
	gas_disk			*h_gasDisk;
	gas_disk			*d_gasDisk;

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
		\param params Vector of parameters of the bodies
		\param coor Vector of coordinates of the bodies
		\param velo Vector of velocities of the bodies
		\param acce Will hold the accelerations for each body
	*/
	cudaError_t call_calculate_drag_accel_kernel(ttt_t time, const param_t* params, const vec_t* coor, const vec_t* velo, vec_t* acce);

	//! Calls the kernel that calculates the acceleration due to type I migration.
	/*
		\param time The actual time of the simulation
		\param params Vector of parameters of the bodies
		\param coor Vector of coordinates of the bodies
		\param velo Vector of velocities of the bodies
		\param acce Will hold the accelerations for each body
	*/
	cudaError_t call_calculate_migrateI_accel_kernel(ttt_t time, param_t* params, const vec_t* coor, const vec_t* velo, vec_t* acce);
};

static __host__ __device__ void		shift_into_range(var_t lower, var_t upper, var_t* value);
static __host__ __device__ vec_t	cross_product(const vec_t* v, const vec_t* u);
static __host__ __device__ var_t	dot_product(const vec_t* v, const vec_t* u);
static __host__ __device__ var_t	norm2(const vec_t* v);
static __host__ __device__ var_t	norm(const vec_t* v);
static __host__ __device__ vec_t	circular_velocity(var_t mu, const vec_t* rVec);
static __host__ __device__ vec_t	gas_velocity(var2_t eta, var_t mu, const vec_t* rVec);
static __host__ __device__ var_t	gas_density_at(const gas_disk* gasDisk, const vec_t* rVec);
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
static __host__ __device__ var_t	reduction_factor(const gas_disk* gasDisk, ttt_t t);
static __host__ __device__ var_t	midplane_density(const gas_disk* gasDisk, var_t r);
static __host__ __device__ var_t	typeI_migration_time(const gas_disk* gasDisk, var_t C, var_t O, var_t ar, var_t er, var_t h);
static __host__ __device__ var_t	typeI_eccentricity_damping_time(var_t C, var_t O, var_t ar, var_t er, var_t h);
