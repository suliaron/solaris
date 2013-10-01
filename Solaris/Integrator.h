#ifndef INTEGRATOR_H_
#define INTEGRATOR_H_

#include <string>

class Integrator {
public:
	Integrator();

	std::string	name;
	std::string	xsitype;
	double		accuracy;
	double		epsilon;
};

#endif
