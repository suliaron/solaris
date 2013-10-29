#pragma once

#include "cuda_runtime.h"
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>

typedef double ttt_t;	// Type of time variables
typedef double var_t;	// Type of variables
typedef double4 vec_t;   // Type of vectors, use 4dim for coalesced memory access
typedef bool bool_t;	// Type of booleans
typedef int int_t;
typedef int2 int2_t;

typedef thrust::host_vector<var_t> host_var_t;
typedef thrust::device_vector<var_t> device_var_t;

typedef thrust::host_vector<vec_t> host_vec_t;
typedef thrust::device_vector<vec_t> device_vec_t;

typedef thrust::host_vector<int_t> host_int_t;
typedef thrust::device_vector<int_t> device_int_t;

typedef thrust::device_vector<int2_t> device_int2_t;