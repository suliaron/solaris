#pragma once

#include "config.h"

#include "integrator.h"
#include "ode.h"

template <int RKOrder>
class rungekuttanystrom : public integrator
{
public:
	static var_t a[];
	static var_t bh[];
	static var_t b[];
	static var_t bd[];
	static ttt_t c[];
	static var_t err;

private:
	bool adaptive;
	var_t tolerance;
	
	std::vector<d_var_t> d_k;		// Differentials on the device
	std::vector<d_var_t> d_ytemp;
	d_var_t d_yhtemp;	
	d_var_t d_dytemp;
	d_var_t d_err;

public:
	rungekuttanystrom(ode& f, ttt_t dt, bool adaptive, var_t tolerance);
	~rungekuttanystrom();

	ttt_t step();
};