#pragma once

#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "thrust\device_vector.h"
#include "thrust\host_vector.h"

//
//#define NVAR 2			// Number of variables per body (coord, veloc)
//#define NPAR 2			// Number of parameters per body
#define NDIM 4				// Number of dimensions, 4 to coalesce memory copies
//
#define NTILE 256
//
//#define K2 (var_t)0.0002959122082855911025
//#define MASS_SUN (var_t)1.9891E+30
//#define DT (var_t)0.1

typedef double		ttt_t;		// Type of time variables
typedef double		var_t;		// Type of variables
typedef double4		vec_t;		// Type of vectors
typedef bool		bool_t;		// Type of boolean variables
typedef int			int_t;		// Type of integer variables
typedef int2		int2_t;		// Type of integer tuples variables

typedef thrust::host_vector<var_t>		h_var_t;
typedef thrust::device_vector<var_t>	d_var_t;

typedef thrust::host_vector<vec_t>		h_vec_t;
typedef thrust::device_vector<vec_t>	d_vec_t;

typedef thrust::host_vector<int_t>		host_int_t;
typedef thrust::device_vector<int_t>	device_int_t;

typedef thrust::device_vector<int2_t>	device_int2_t;

#define	NVAR		2			// Number of vector variables per body (coordinate, velocity)
#define NPAR		2			// Number of parameters per body (mass, radius)

#define K			(var_t)0.01720209895
#define K2			(var_t)0.0002959122082855911025

#define	SQR(x)		(x)*(x)
