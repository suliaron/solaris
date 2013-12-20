#pragma once

#include "config.h"
#include "integrator.h"
#include "ode.h"


class euler : public integrator
{
private:
	std::vector<d_var_t> d_dy;		// Differentials on the device

public:
	euler(ode& f, ttt_t dt);
	~euler();

	ttt_t step();
};