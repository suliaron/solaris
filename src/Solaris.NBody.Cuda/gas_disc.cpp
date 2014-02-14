#include "gas_disc.h"

gas_disc::gas_disc() 
{

	rho.x = rho.y = 0.0;
	sch.x = sch.y = 0.0;
	eta.x = eta.y = 0.0;
	tau.x = tau.y = 0.0;

	gas_decrease = CONSTANT;
	t0 = t1 = timeScale = 0.0;
}

gas_disc::gas_disc(var2_t rho, var2_t sch, var2_t eta, var2_t tau, gas_decrease_t gas_decrease, ttt_t t0, ttt_t t1, ttt_t timeScale) :
	rho(rho),
	sch(sch),
	eta(eta),
	tau(tau),
	gas_decrease(gas_decrease),
	t0(t0),
	t1(t1),
	timeScale(timeScale)
{}

__host__ __device__
var_t	gas_disc::reduction_factor(ttt_t t)
{
	switch (gas_decrease) 
	{
	case CONSTANT:
		return 1.0;
		break;
	case LINEAR:
		if (t <= t0) {
			return 1.0;
		}
		else if (t > t0 && t <= t1) {
			return 1.0 - (t - t0)/(t1 - t0);
		}
		else {
			return 0.0;
		}
		break;
	case EXPONENTIAL:
		return exp(-t/timeScale);
		break;
	default:
		return 1.0;
		break;
	}
}

std::ostream& operator<<(std::ostream& output, gas_disc gasDisc)
{
	output << "eta: " << gasDisc.eta.x << ", " << gasDisc.eta.y << std::endl;
	output << "rho: " << gasDisc.rho.x << ", " << gasDisc.rho.y << std::endl;
	output << "sch: " << gasDisc.sch.x << ", " << gasDisc.sch.y << std::endl;
	output << "tau: " << gasDisc.tau.x << ", " << gasDisc.tau.y << std::endl;
	output << "gas_decrease: " << gasDisc.gas_decrease << std::endl;
	output << "          t0: " << gasDisc.t0 << " [d]" << std::endl;
	output << "          t1: " << gasDisc.t1 << " [d]" << std::endl;
	output << "   timeScale: " << gasDisc.timeScale << " [d]" << std::endl;

	return output;
}
