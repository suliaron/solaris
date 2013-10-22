#include <math.h>

#include "Characteristics.h"
#include "Constants.h"

Characteristics::Characteristics()
{
	stokes		= 0.0;
	absVisMag	= 0.0;

	mass		= 0.0;
	radius		= 0.0;
	density		= 0.0;
}

Characteristics::Characteristics(double mass)
{
	this->mass		= mass;

	this->stokes	= 0.0;
	this->absVisMag	= 0.0;
	this->radius	= 0.0;
	this->density	= 0.0;
}

Characteristics::Characteristics(const Characteristics &characteristics)
{
	mass		= characteristics.mass;
	radius		= characteristics.radius;
	density		= characteristics.density;
	stokes		= characteristics.stokes;
	absVisMag	= characteristics.absVisMag;
	// TODO: check
	componentList = characteristics.componentList;
}

/// <summary>
/// Computes the density of the body assuming a spherical form.
/// The radius must be defined.
/// </summary>
/// <returns>The volume of the body</returns>
double Characteristics::CalculateDensity()
{
	return mass / CalculateVolume();
}

double Characteristics::CalculateMass()
{
    return density*(4.0 / 3.0 * Constants::Pi * (radius*radius*radius));
}

/// <summary>
/// Computes the radius of the body assumed that the body is spherical in shape.
/// The mass and the volume density must be defined.
/// </summary>
/// <returns>The radius of the body</returns>
double Characteristics::CalculateRadius()
{
	return pow(3.0/(4.0*Constants::Pi)*mass/density, 1.0/3.0);
}

/// <summary>
/// Computes the volume of the body assuming a spherical form.
/// The radius must be defined.
/// </summary>
/// <returns>The volume of the body</returns>
double Characteristics::CalculateVolume()
{
	return 4.0 / 3.0 * Constants::Pi * pow(radius, 3);
}

/// <summary>
/// Coefficient for the Epstein-drag force (gamma = m/rho*R = 4/3*Pi*R^2 )
/// </summary>
double Characteristics::GammaEpstein()
{
    return 1.0/(density*radius);
}

/// <summary>
/// Coefficient for the Stokes-drag force (= 3/8 * Cd /(rho*R) )
/// </summary>
//[XmlIgnore]
double Characteristics::GammaStokes()
{
    return (3.0/8.0)*stokes/(density*radius);
}

double Characteristics::CalculateMass(double density, double radius)
{
	double volume = 4.0 / 3.0 * Constants::Pi * (radius*radius*radius);
    return density*volume;
}
