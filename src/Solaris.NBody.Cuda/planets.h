#pragma once

#include <string>

#include "ode.h"
#include "config.h"
#include "gas_disc.h"
#include "number_of_bodies.h"

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

	typedef thrust::host_vector<param_t>	h_param_t;
	typedef thrust::device_vector<param_t>	d_param_t;

	typedef enum migration
	{
		type1,
		type2
	} migration_t;

private:
	number_of_bodies	bodies;
	gas_disc			*gasDisc;
	gas_disc			*d_gasDisc;
	d_vec_t				acceGasDrag;

public:
	planets(number_of_bodies bodies, gas_disc* gasDisc);
	~planets();

	void round_up_n();
	void allocate_vectors();

	void transform_to_barycentric();
	void transform_to_objectcentric();

	void calculate_dy(int i, int r, ttt_t t, const d_var_t& p, const std::vector<d_var_t>& y, d_var_t& dy);

	void load(string filename);
	int print_positions(ostream& sout);
private:
	//! Calls the kernel that calculates the accelerations from gravitational
	/*  interactions. Done in tiles.
		\param p Vector of parameters of the bodies
		\param c Vector of coordinates of the bodies
		\param bounds Vector of indices limiting the interacting pairs
		\param atemp Will hold the accelerations for each body per each tile
	*/
	cudaError_t	call_calculate_grav_accel_kernel(const planets::param_t* p, const vec_t* c, vec_t* atemp);

	//! Calls the kernel that calculates the acceleration due to drag force on bodies
	cudaError_t call_calculate_drag_accel_kernel(ttt_t time, const gas_disc* gasDisc, const planets::param_t* params, const vec_t* coor, const vec_t* velo, vec_t* acce);

	//cudaError_t call_calculate_epheremis_kernel(const param_t* p, const vec_t* c, const vec_t* v, int2_t bounds);

	//cudaError_t call_calculate_migration_accel_kernel(const param_t* p, const vec_t* c, const vec_t* v, int2_t bounds, gas_disc* gasDisc, vec_t* atemp);

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
static __host__ __device__ int_t	calculate_phase(var_t mu, const planets::orbelem_t* oe, vec_t* rVec, vec_t* vVec);
static __host__ __device__ int_t	calculate_sma_ecc(var_t mu, const vec_t* coor, const vec_t* velo, var_t* sma, var_t* ecc);
static __host__ __device__ int_t	calculate_orbelem(var_t mu, const vec_t* coor, const vec_t* velo, planets::orbelem_t* orbelem);
static __host__ __device__ var_t	orbital_period(var_t mu, var_t sma);
static __host__ __device__ var_t	orbital_frequency(var_t mu, var_t sma);
//static __host__ __device__ var_t	calculate_gamma_stokes(var_t stokes, var_t density, var_t radius);
