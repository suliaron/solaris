#pragma once

#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "thrust\device_vector.h"
#include "thrust\host_vector.h"


// General settings for the integrator

//! Type of time variables
typedef double		ttt_t;
//! Type of variables
typedef double		var_t;
//! Type of tuple
typedef double2		var2_t;
//! Type of vectors
typedef double4		vec_t;
//! Type of boolean variables
typedef bool		bool_t;
//! Type of integer variables
typedef int			int_t;
//! Type of integer tuples variables
typedef int2		int2_t;

typedef thrust::host_vector<var_t>		h_var_t;
typedef thrust::device_vector<var_t>	d_var_t;

typedef thrust::host_vector<vec_t>		h_vec_t;
typedef thrust::device_vector<vec_t>	d_vec_t;

typedef thrust::host_vector<int_t>		h_int_t;
typedef thrust::device_vector<int_t>	d_int_t;

typedef thrust::device_vector<int2_t>	d_int2_t;

// NBody settings

#define NDIM		4		// Number of dimensions, 4 to coalesce memory copies
#define NTILE		256

#define	NVAR		2		// Number of vector variables per body (coordinate, velocity)
#define NPAR		2		// Number of parameters per body (mass, radius)

#define K			(var_t)0.01720209895
#define K2			(var_t)0.0002959122082855911025

#define	PI			(var_t)3.1415926535897932384626
#define	TWOPI		(var_t)6.2831853071795864769253
#define	TORAD		(var_t)0.0174532925199432957692
#define TODEG		(var_t)57.295779513082320876798

// These macro functions must be enclosed in parentheses in order to give
// correct results in the case of a division i.e. 1/SQR(x) -> 1/((x)*(x))
#define	SQR(x)		((x)*(x))
#define	CUBE(x)		((x)*(x)*(x))
