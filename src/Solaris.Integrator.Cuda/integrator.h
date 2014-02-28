#pragma once

#include "config.h"
#include "ode.h"

class integrator
{
protected:
	ode& f;
	ttt_t dt;
	ttt_t dt_try;
	ttt_t dt_did;

	int_t n_failed_step;
	int_t n_step;

public:
	integrator(ode& f, ttt_t dt);
	~integrator();

	int_t get_n_failed_step();
	int_t get_n_step();

	virtual ttt_t step() = 0;	
};