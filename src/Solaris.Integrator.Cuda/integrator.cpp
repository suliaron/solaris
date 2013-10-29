#include "integrator.h"

integrator::integrator(ode& f, ttt_t dt)
	: f(f),	dt(dt)
{
}

integrator::~integrator()
{
}
