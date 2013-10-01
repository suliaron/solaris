#ifndef POWERLAW_H_
#define POWERLAW_H_

class PowerLaw
{
public:
	PowerLaw();
	PowerLaw(double c, double index);

	double	Evaluate(double x);

	double	c;
	double	index;
};

#endif
