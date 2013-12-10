#pragma once 

#include "config.h"

struct sum_vec_functor
{
	var_t a;

	sum_vec_functor(var_t a) : a(a) {}
	__host__ __device__ var_t operator()(var_t x, var_t y);
};

struct absdiff_vec_functor
{
	var_t a;

	absdiff_vec_functor(var_t a) : a(a) {}
	__host__ __device__ var_t operator()(var_t x, var_t y);
};


var_t randf(var_t min, var_t max);
void copy_vec(d_var_t& xout, const d_var_t& xin);
void sum_vec(d_var_t& xout, const d_var_t& xin, const d_var_t& yin, var_t a);
void absdiff_vec(d_var_t& xout, const d_var_t& xin, const d_var_t& yin, var_t a);
var_t max_vec(const d_var_t& x);
