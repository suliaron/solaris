#include <math.h>
#include "PowerLaw.h"

PowerLaw::PowerLaw()
{
	c = 0.0;
	index = 0.0;
}

PowerLaw::PowerLaw(double co, double idx)
{
	c = co;
	index = idx;
}

// TODO: make efficency/speed test comparing the active code with the inactive one
double PowerLaw::Evaluate(double x)
{
    return this->c*pow(x, this->index);

	//if      (this->index == 0)
	//	return this->c;
	//else if (this->index == 0.5)
	//	return this->c * sqrt(x);
	//else if (this->index == 1)
	//	return this->c * x;
	//else if (this->index == 2)
	//	return this->c * x * x;
	//else if (this->index == -0.5)
	//	return this->c / sqrt(x);
	//else if (this->index == -1)
	//	return this->c / x;
	//else if (this->index == -2)
	//	return this->c / (x * x);
	//else return this->c*pow(x, this->index);
}
