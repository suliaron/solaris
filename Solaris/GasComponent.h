#ifndef GASCOMPONENT_H_
#define GASCOMPONENT_H_

#include "GasDecreaseType.h"
#include "PowerLaw.h"
#include "Vector.h"

class GasComponent
{
public:
	GasComponent();

	double		alpha;
	double      meanMolecularWeight;
    double      particleDiameter;

	GasDecreaseType type;
	double		timeScale;
	double		t0;
	double		t1;
	double		innerEdge;
	double		a;

	PowerLaw	eta;
	PowerLaw	tau;
	PowerLaw	scaleHeight;
	PowerLaw	density;
    PowerLaw    meanFreePath;
    PowerLaw    temperature;

	double	ReductionFactor(const double t);
	double	MidplaneDensity(const double r);
	double	GasDensityAt(double r, double z);
	Vector	GasVelocity(double mu, double r, double alpha);
	Vector	CircularVelocity(double mu, double r, double alpha);

    double  MeanFreePath_SI(const double rho);
    double  MeanFreePath_CMU(const double rho);
    
	double  SoundSpeed_SI(const double temperature);
    double  SoundSpeed_CMU(const double temperature);

	double  Temperature_SI(const double mC, const double r);
    double  Temperature_CMU(const double mC, const double r);

	double  MeanThermalSpeed_SI(const double mC, const double r);
	double  MeanThermalSpeed_CMU(const double mC, const double r);
};

#endif
