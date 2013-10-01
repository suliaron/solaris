#ifndef CHARACTERISTICS_H_
#define CHARACTERISTICS_H_

#include <list>
#include "Component.h"

class Characteristics
{
public:

	Characteristics();
	Characteristics(double mass);
	Characteristics(const Characteristics &characteristics);

	double CalculateDensity();
	double CalculateMass();
	double CalculateRadius();
	double CalculateVolume();
	double GammaStokes();
	double GammaEpstein();

	static double CalculateMass(double density, double radius);

	std::list<Component> componentList;

	double	mass;
	double	radius;
	double	density;
	double	stokes;
	double	absVisMag;

};

#endif