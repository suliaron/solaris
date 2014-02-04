#pragma once

#include "config.h"

class gas_disc
{
public:
	gas_disc();
	gas_disc(var2_t rho, var2_t sch, var2_t eta, var2_t tau);

	var2_t	rho;
	var2_t	sch;
	var2_t	eta;
	var2_t	tau;

	// Input/Output streams
	friend std::ostream& operator<<( std::ostream&, gas_disc);
};
