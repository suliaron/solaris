#ifndef GASCOMPONENT_H_
#define GASCOMPONENT_H_

#include "GasDecreaseType.h"
#include "PowerLaw.h"

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

	PowerLaw	eta;
	PowerLaw	scaleHeight;
	PowerLaw	tau;
	PowerLaw	gasDensityFunction;
    PowerLaw    meanFreePath;
    PowerLaw    temperature;

	double	FactorAt(const double t);
	double	MidplaneDensity(const double r);

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
