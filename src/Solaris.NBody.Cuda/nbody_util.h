#pragma once

#include "config.h"
#include "planets.h"

__host__ __device__ void	shift_into_range(var_t lower, var_t upper, var_t* value);
__host__ __device__ vec_t	cross_product(const vec_t* v, const vec_t* u);
__host__ __device__ var_t	dot_product(const vec_t* v, const vec_t* u);
__host__ __device__ var_t	norm2(const vec_t* v);
__host__ __device__ var_t	norm(const vec_t* v);
__host__ __device__ vec_t	circular_velocity(var_t mu, vec_t* rVec);
__host__ __device__ vec_t	gas_velocity(var2_t eta, var_t mu, vec_t* rVec);
__host__ __device__ var_t	gas_density_at(const planets::gaspar_t* gaspar, const vec_t* rVec);
__host__ __device__ var_t	calculate_kinetic_energy(const vec_t* vVec);
__host__ __device__ var_t	calculate_potential_energy(var_t mu, const vec_t* rVec);
__host__ __device__ var_t	calculate_energy(var_t mu, const vec_t* rVec, const vec_t* vVec);
__host__ __device__ int_t	kepler_equation_solver(var_t ecc, var_t mean, var_t eps, var_t* E);
__host__ __device__ int_t	calculate_phase(var_t mu, const planets::orbelem_t* oe, vec_t* rVec, vec_t* vVec);
__host__ __device__ int_t	calculate_sma_ecc(var_t mu, const vec_t* coor, const vec_t* velo, var_t* sma, var_t* ecc);
__host__ __device__ int_t	calculate_orbelem(var_t mu, const vec_t* coor, const vec_t* velo, planets::orbelem_t* orbelem);
