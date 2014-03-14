#include "gas_disk.h"

gas_disk::gas_disk() 
{

	rho.x = rho.y = 0.0;
	sch.x = sch.y = 0.0;
	eta.x = eta.y = 0.0;
	tau.x = tau.y = 0.0;

	gas_decrease = CONSTANT;
	t0 = t1 = timeScale = 0.0;
}

gas_disk::gas_disk(var2_t rho, var2_t sch, var2_t eta, var2_t tau, gas_decrease_t gas_decrease, ttt_t t0, ttt_t t1, ttt_t timeScale) :
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
var_t	gas_disk::reduction_factor(ttt_t t)
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

std::ostream& operator<<(std::ostream& output, gas_disk gasDisk)
{
	output << "eta: " << gasDisk.eta.x << ", " << gasDisk.eta.y << std::endl;
	output << "rho: " << gasDisk.rho.x << ", " << gasDisk.rho.y << std::endl;
	output << "sch: " << gasDisk.sch.x << ", " << gasDisk.sch.y << std::endl;
	output << "tau: " << gasDisk.tau.x << ", " << gasDisk.tau.y << std::endl;
	output << "gas_decrease: " << gasDisk.gas_decrease << std::endl;
	output << "          t0: " << gasDisk.t0 << " [d]" << std::endl;
	output << "          t1: " << gasDisk.t1 << " [d]" << std::endl;
	output << "   timeScale: " << gasDisk.timeScale << " [d]" << std::endl;

	return output;
}
