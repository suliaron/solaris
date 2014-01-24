#include "gas_disc.h"

gas_disc::gas_disc() 
{
	rho.x = rho.y = 0.0;
	sch.x = sch.y = 0.0;
	eta.x = eta.y = 0.0;
	tau.x = tau.y = 0.0;
}

gas_disc::gas_disc(var2_t rho, var2_t sch, var2_t eta, var2_t tau) :
	rho(rho),
	sch(sch),
	eta(eta),
	tau(tau)
{}

std::ostream& operator<<(std::ostream& output, gas_disc gasDisc)
{
	output << "eta: " << gasDisc.eta.x << ", " << gasDisc.eta.y << std::endl;
	output << "rho: " << gasDisc.rho.x << ", " << gasDisc.rho.y << std::endl;
	output << "sch: " << gasDisc.sch.x << ", " << gasDisc.sch.y << std::endl;
	output << "tau: " << gasDisc.tau.x << ", " << gasDisc.tau.y << std::endl;

	return output;
}
