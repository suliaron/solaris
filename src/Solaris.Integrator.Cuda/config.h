//#pragma once
//
//#include "cuda_runtime.h"
//#include <thrust/host_vector.h>
//#include <thrust/device_vector.h>
//
//typedef double	ttt_t;		// Type of time variables
//typedef double	var_t;		// Type of variables
//typedef double4 vec_t;		// Type of vectors, use 4dim for coalesced memory access
//typedef bool	bool_t;		// Type of booleans
//typedef int		int_t;
//typedef int2	int2_t;
//typedef int4	int4_t;
//
//typedef thrust::host_vector<var_t>		host_var_t;
//typedef thrust::device_vector<var_t>	device_var_t;
//
//typedef thrust::host_vector<vec_t>		host_vec_t;
//typedef thrust::device_vector<vec_t>	device_vec_t;
//
//typedef thrust::host_vector<int_t>		host_int_t;
//typedef thrust::device_vector<int_t>	device_int_t;
//
//typedef thrust::device_vector<int2_t>	device_int2_t;

#pragma once

#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "thrust\device_vector.h"
#include "thrust\host_vector.h"

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
