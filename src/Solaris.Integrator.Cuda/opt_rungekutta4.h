#pragma once

#include "config.h"
#include "integrator.h"
#include "ode.h"

class opt_rungekutta4 : public integrator
{
public:
	static var_t a[];
	static var_t b[];
	static var_t bh[];
	static ttt_t c[];

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
	std::vector<std::vector <d_var_t> >	d_k;
	//! Holds the temporary solution approximation along the step
	std::vector<d_var_t>				d_ytemp;
	//! Holds the leading local truncation error for each variable
	std::vector<d_var_t>				d_err;

	dim3	grid;
	dim3	block;

	void calculate_grid(int nData, int threads_per_block);
	void call_calc_ytemp_for_kr_kernel(int r);
	void call_calc_yHat_kernel();
	void call_calc_y_kernel();
	void call_calc_k4_sub_k5_kernel();

public:
	opt_rungekutta4(ode& f, ttt_t, bool adaptive, var_t tolerance);
	~opt_rungekutta4();

	ttt_t step();

};
