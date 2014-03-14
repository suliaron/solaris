#pragma once

#include "config.h"
#include "integrator.h"
#include "ode.h"

class rkn76 : public integrator
{
public:
	//! Constants representing the seventh order RKN formulae (Celestial Mechanics, Vol. 18(1978), 223-232.)
	static var_t c[];
	static var_t a[];
	static var_t bh[];
	static var_t b[];

private:
	//! The order of the embedded RK formulae
	int		RKOrder;
	//! The maximum number of the force calculation
	int		r_max;
	//! True if the method estimates the error and accordingly adjusts the step-size
	bool	adaptive;
	//! The maximum of the allowed local truncation error
	var_t	tolerance;

	//! Holds the derivatives for the differential equations
	std::vector<std::vector<d_var_t>>	d_k;
	//! Holds the temporary solution approximation along the step
	std::vector<d_var_t>				d_ytemp;
	//! Holds the leading local truncation error for each variable
	std::vector<d_var_t>				d_err;

	dim3	grid;
	dim3	block;

	void calculate_grid(int nData, int threads_per_block);
	void call_calc_ytemp_for_kr_kernel(int r);
	void call_calc_y_kernel();
	void call_calc_k8_sub_k9_kernel();

public:
	rkn76(ode& f, ttt_t, bool adaptive, var_t tolerance);
	~rkn76();

	ttt_t step();

};
