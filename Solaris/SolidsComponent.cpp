#include <cstdio>
#include "SolidsComponent.h"

SolidsComponent::SolidsComponent()
{
	_iceCondensationFactor = 3.0;
	// TODO: Convert 7 g/cm2 to Solar/AU2
	_solidsDensityFunction = PowerLaw(7, -1.5);
}
