#include "Integrator.h"

Integrator::Integrator()
{
	name		= "rungekutta78";
	accuracy = -10.0;
	epsilon	 = pow(10, accuracy);
}
