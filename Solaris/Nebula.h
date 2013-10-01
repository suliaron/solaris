#ifndef NEBULA_H_
#define NEBULA_H_

#include <string>
#include "SolidsComponent.h"
#include "GasComponent.h"

class Nebula {
public:
	Nebula();

	std::string name;
	std::string description;
	std::string file;
	std::string path;

	double	massFactor;
	double	gasToDustRatio;
	double	snowLine;

	SolidsComponent	solidsComponent;
	GasComponent	gasComponent;
};

#endif