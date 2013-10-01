#ifndef FARGOPARAMETERS_H_
#define FARGOPARAMETERS_H_

#include <string>

class FargoParameters
{
public:
	FargoParameters();

	int			ReadConfigFile(std::string& path);
	int			ParseConfig(bool verbose);

	double		aspectRatio;         // Thickness over Radius in the disc
	double		sigma0;              // Surface Density at r=1
	double		alphaViscosity;      // Uniform kinematic viscosity
	double		sigmaSlope;          // Slope of surface density profile.
	double		flaringIndex;
	std::string	excludeHill;

//Planet parameters

	std::string	planetConfig;
	double		thicknessSmoothing;  // Smoothing parameters in disk thickness

// Numerical method parameters

	std::string	transport;
	std::string	innerBoundary;		 // choose : OPEN or RIGID or NONREFLECTING
	std::string	disk;
	double		omegaFrame;
	std::string	frame;
	std::string	indirectTerm;

// Mesh parameters

	int			nRad;                // Radial number of zones
	int			nSec;                // Azimuthal number of zones (sectors)
	double		rMin;                // Inner boundary radius
	double		rMax;                // Outer boundary radius
	std::string	radialSpacing;       // Zone interfaces evenly spaced

// Output control parameters

	int			nTot;                // Total number of time steps
	int			nInterm;             // Time steps between outputs
	double		dT;                  // Time step length. 2PI = 1 orbit
	std::string	outputDir;

// Viscosity damping due to a dead zone

	double		viscModR1;           // Inner radius of dead zone
	double		viscModDeltaR1;      // Width of viscosity transition at inner radius
	double		viscModR2;           // Outer radius of dead zone
	double		viscModDeltaR2;      // Width of viscosity transition at outer radius
	double		viscMod;             // Viscosity damp

private:
	int			SetParameter(std::string& key, std::string& value, const bool verbose);

	std::string config;
};

#endif