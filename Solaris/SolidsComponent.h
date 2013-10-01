#ifndef SOLIDSCOMPONENT_H_
#define SOLIDSCOMPONENT_H_

#include "PowerLaw.h"

class SolidsComponent
{
public:
	SolidsComponent();

	//TODO remove inline functions
	inline double GetIceCondensationFactor() { return _iceCondensationFactor; }
	inline void SetIceCondensationFactor(double iceCondensationFactor) { _iceCondensationFactor = iceCondensationFactor; }

	inline PowerLaw* GetSolidsDensityFunction() { return &_solidsDensityFunction; }
	inline void SetSolidsDensityFunction(PowerLaw solidsDensityFunction) {
		_solidsDensityFunction = solidsDensityFunction;
	}

private:
	double		_iceCondensationFactor;
	PowerLaw	_solidsDensityFunction;
};

#endif