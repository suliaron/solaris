#pragma once

#include <thrust/device_vector.h>
#include "config.h"
#include "ode.h"
#include "integrator.h"


class euler : public integrator
{
private:
	std::vector<device_var_t> d_dy;		// Differentials on the device

public:
	euler(ode& f, ttt_t dt);
	~euler();

	ttt_t step();
};