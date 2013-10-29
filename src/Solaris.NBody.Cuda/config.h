#pragma once

#include "cuda_runtime.h"
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>

#define NVAR 2			// Number of variables per body (coord, veloc)
#define NPAR 2			// Number of parameters per body
#define NDIM 4			// Number of dimensions, 4 to coalesce memory copies

#define NTILE 256

#define K2 (var_t)0.0002959122082855911025
#define MASS_SUN (var_t)1.9891E+30
#define DT (var_t)0.1