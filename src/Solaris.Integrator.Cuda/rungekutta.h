#pragma once

#include "config.h"

#include "integrator.h"
#include "ode.h"

template <int RKOrder>
class rungekutta : public integrator
{
public:
	static var_t a[];
	static var_t b[];
	static ttt_t c[];

private:
	bool adaptive;
	var_t tolerance;
	
	std::vector<device_var_t> d_ytemp;				// Values on the device
	std::vector< std::vector< device_var_t> > d_k;		// Differentials on the device

public:
	rungekutta(ode& f, ttt_t dt, bool adaptive, var_t tolerance);
	~rungekutta();

	ttt_t step();
};
