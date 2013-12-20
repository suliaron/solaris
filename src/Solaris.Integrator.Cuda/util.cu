#include "util.h"
#include <cstdlib>

#include "config.h"

var_t sum_vec_functor::operator()(var_t x, var_t y)
{
	return x + a * y;
}

var_t absdiff_vec_functor::operator()(var_t x, var_t y)
{
	return a * fabs(x - y);
}

var_t randf(var_t min, var_t max)
{
	return min + (max - min) * (var_t)rand() / RAND_MAX;
}

void copy_vec(d_var_t& xout, const d_var_t& xin)
{
	thrust::copy(xin.begin(), xin.end(), xout.begin());
}

// Compute xout = x + a * y
void sum_vec(d_var_t& xout, const d_var_t& xin, const d_var_t& yin, var_t a)
{
	thrust::transform(
			xin.begin(), xin.end(),
			yin.begin(),
			xout.begin(),
			sum_vec_functor(a));
}


// Compute xout = a * |x - y|
void absdiff_vec(d_var_t& xout, const d_var_t& xin, const d_var_t& yin, var_t a)
{
	thrust::transform(
			xin.begin(), xin.end(),
			yin.begin(),
			xout.begin(),
			absdiff_vec_functor(a));
}

var_t max_vec(const d_var_t& x)
{
	return *thrust::max_element(x.begin(), x.end());
}