#pragma once

#include "config.h"
#include "integrator.h"
#include "ode.h"

class opt_midpoint_method : public integrator
{
public:
	static var_t a[];
	static var_t b[];
	static ttt_t c[];

private:
	int		RKOrder;
	bool	adaptive;
	var_t	tolerance;

	std::vector<std::vector<d_var_t>>	d_k;
	std::vector<d_var_t>				d_ytemp;

	dim3	grid;
	dim3	block;

	void calculate_grid(int nData, int threads_per_block);
	void call_calc_2nd_arg_of_k2_kernel();
	void call_calc_final_dependent_variables_kernel();

public:
	opt_midpoint_method(ode& f, ttt_t dt, bool adaptive, var_t tolerance);
	~opt_midpoint_method();

	ttt_t	step();
};
