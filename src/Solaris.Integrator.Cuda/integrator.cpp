#include "integrator.h"

integrator::integrator(ode& f, ttt_t dt)
	: f(f),	dt(dt)
{
	n_failed_step	= 0;
	n_step			= 0;

	dt_try			= 0.0;
	dt_did			= 0.0;
}

integrator::~integrator()
{
}

int_t integrator::get_n_failed_step()
{
	return n_failed_step;
}

int_t integrator::get_n_step()
{
	return n_step;
}
