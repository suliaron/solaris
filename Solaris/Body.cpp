#include <iostream>
#include <algorithm>

#include "Body.h"
#include "Constants.h"
#include "Ephemeris.h"

int Body::_bodyId = 0;

Body::Body() 
{
	_id				= Body::_bodyId++;

	type			= UndefinedBodyType;
	ln				= UndefinedLn;
	mPCOrbitType	= UndefinedMPCOrbitType;

	migrationType	= No;
	migrationStopAt = 0.0;

	phase			= 0;
	orbitalElement	= 0;
	characteristics	= 0;
}

Body::Body(int id) 
{
	_id = id;

	type = UndefinedBodyType;
	ln = UndefinedLn;
	mPCOrbitType = UndefinedMPCOrbitType;

	migrationType = No;
	migrationStopAt = 0.0;

	phase = 0;
	orbitalElement = 0;
	characteristics = 0;
}

Body::Body(BodyType type)
{
	_id = Body::_bodyId++;

	this->type = type;
	ln = UndefinedLn;
	mPCOrbitType = UndefinedMPCOrbitType;

	migrationType = No;
	migrationStopAt = 0.0;

	phase = 0;
	orbitalElement = 0;
	characteristics = 0;
}

Body::Body(BodyType type, std::string name)
{
	_id = Body::_bodyId++;

	this->name = name;

	this->type = type;
	ln = UndefinedLn;
	mPCOrbitType = UndefinedMPCOrbitType;

	migrationType = No;
	migrationStopAt = 0.0;

	phase = 0;
	orbitalElement = 0;
	characteristics = 0;
}

double Body::GetGm() 
{ 
	return Constants::Gauss2*(characteristics->mass); 
}

/// <summary>
/// Computes the orbital period of the body with respect to the central body.
/// It uses the Phase of the object during the computation and the Energy is
/// computed using the one-centrum formalism.
/// If the Energy is positive, zero is returned.
/// </summary>
/// <param name="mu">(Gauss constant)^2*(m1 + m2)</param>
/// <returns>The orbital period of the body.</returns>
double Body::CalculateOrbitalPeriod(double mu)
{
	double h = Ephemeris::CalculateEnergy(mu, *(this->phase));
    if (h > 0.0)
    {
        return -1;
    }
    double a = -mu / (2.0 * h);
    return 2.0*Constants::Pi*sqrt((a*a*a)/mu);
}

void Body::Print()
{
	std::cout << "body.Id: " << this->_id << std::endl;
	if (this->name.length() > 0 )
		std::cout << "body.name: " << this->name << std::endl;
	if (this->designation.length() > 0 )
		std::cout << "body.designation: " << this->designation << std::endl;
	if (this->provisionalDesignation.length() > 0 )
		std::cout << "body.provisionalDesignation: " << this->provisionalDesignation << std::endl;
	if (this->reference.length() > 0 )
		std::cout << "body.reference: " << this->reference << std::endl;
	if (this->opposition.length() > 0 )
		std::cout << "body.opposition: " << this->opposition << std::endl;
	if (this->guid.length() > 0 )
		std::cout << "body.guid: " << this->guid << std::endl;

	std::cout << "body.type: " << this->type << std::endl;
	std::cout << "body.ln: " << this->ln << std::endl;
	std::cout << "body.migration: " << this->migrationType << std::endl;
	std::cout << "body.migrationStopAt : " << this->migrationStopAt << " AU" << std::endl;

	if (phase != 0) {
		std::cout << "body.phase: " << phase << std::endl;
	}
	if (orbitalElement != 0) {
		std::cout << orbitalElement << std::endl;
	}
	if (characteristics != 0) {
		std::cout << "mass:    " << characteristics->mass << std::endl;
		std::cout << "radius:  " << characteristics->radius << std::endl;
		std::cout << "density: " << characteristics->density << std::endl;
		std::cout << "stokes : " << characteristics->stokes << std::endl;
		if (!characteristics->componentList.empty()) {
			std::list<Component>::iterator it;
			for (it = characteristics->componentList.begin(); it != characteristics->componentList.end(); it++) {
					std::cout << it->name << ": " << it->ratio << std::endl;
			}
		}
	}
}
