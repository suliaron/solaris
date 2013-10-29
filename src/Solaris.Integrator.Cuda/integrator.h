#pragma once

#include "config.h"
#include "ode.h"

class integrator
{
protected:
	ode& f;
	ttt_t dt;

public:
	integrator(ode& f, ttt_t dt);
	~integrator();

	virtual ttt_t step() = 0;	
};