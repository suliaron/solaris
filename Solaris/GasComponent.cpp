#include <cstdio>
#include <cmath>

#include "Constants.h"
#include "GasComponent.h"

#define SQR(a)		((a)*(a))

GasComponent::GasComponent()
{
	// TODO: create an xml attribute for innerEdge in the xml input file.
	innerEdge	        = 10.0 * Constants::SolarRadiusToAu;
	alpha		        = 2.0e-3;
	
	// TODO: create an xml attribute for meanMolecularWeight and particleRadius
    meanMolecularWeight = 2.3;
    particleDiameter    = 3.0e-10;  // m

	timeScale 	        = 0.0;
	t0			        = 0.0;
	t1			        = 0.0;
	type		        = Constant;

	eta					= PowerLaw(0.0019, 0.5); 
	tau					= PowerLaw(2.0/3.0, 2.0);
	scaleHeight			= PowerLaw(0.02, 1.25);
	density				= PowerLaw(1.4e-9 * Constants::GramPerCm3ToSolarPerAu3, -2.75);
	a					= density.Evaluate(innerEdge) / SQR(SQR(innerEdge));
    
    double Clambda      = meanMolecularWeight * Constants::ProtonMass_CMU / (sqrt(2.0) * Constants::Pi * SQR(particleDiameter * Constants::MeterToAu) * density.c);
    double plambda      =-density.index;
    meanFreePath        = PowerLaw(Clambda, plambda);

}

double GasComponent::FactorAt(const double t)
{
	switch (type) 
	{
	case Constant:
		return 1.0;
		break;
	case Linear:
		if (t <= t0) {
			return 1.0;
		}
		else if (t > t0 && t <= t1) {
			return 1.0 - (t - t0)/(t1 - t0);
		}
		else {
			return 0.0;
		}
		break;
	case Exponential:
		return exp(-t/timeScale);
		break;
	default:
		return 1.0;
		break;
	}
}

double GasComponent::MidplaneDensity(const double r)
{
	double a1 = this->density.Evaluate(r);
	double a2 = this->scaleHeight.Evaluate(r);
	double a3 = a1*a2*Constants::SqrtTwoPi;
	return a3;
}

double GasComponent::GasDensityAt(double r, double z)
{
//	static double a = density.Evaluate(innerEdge) / SQR(SQR(innerEdge));

	double h = scaleHeight.Evaluate(r);
	double arg = SQR(z/h);
	if (r > innerEdge) {
		return density.Evaluate(r) * exp(-arg);
	}
	else {
		return a * SQR(SQR(r)) * exp(-arg);
	}
}

Vector GasComponent::GasVelocity(double mu, double r, double alpha)
{
	double f = sqrt(1.0 - 2.0*eta.Evaluate(r));
	return f * CircularVelocity(mu, r, alpha);
	//Vector v = f * CircularVelocity(mu, r, alpha);
	//v.x *= f;
	//v.y *= f;
	//v.z *= f;

	//return v;
}

Vector GasComponent::CircularVelocity(double mu, double r, double alpha)
{
	double v = sqrt( mu/r );
	return Vector(-v*sin(alpha), v*cos(alpha), 0.0);
}

double  GasComponent::MeanFreePath_SI(const double rho)
{
    static const double sqrtTwoPi = 4.4428829381583662470158809900607;
    static double a1 = meanMolecularWeight * Constants::ProtonMass_SI;
    static double d2 = SQR(particleDiameter);

    double result = a1 / sqrtTwoPi;
    result /= (d2*rho);

    return result;
}

double  GasComponent::MeanFreePath_CMU(const double rho)
{
    // http://www.terrapub.co.jp/journals/EPS/pdf/2003/5505/55050263.pdf
    double sqrtTwoPi = 4.4428829381583662470158809900607;
    double a1 = meanMolecularWeight * Constants::ProtonMass_CMU;
    double d2 = SQR(particleDiameter * Constants::MeterToAu);

    double result = a1 / sqrtTwoPi;
    result /= (d2*rho);

    return result;
}

double  GasComponent::SoundSpeed_SI(const double temperature)
{
    double result = Constants::Boltzman_SI * temperature;
    result /= meanMolecularWeight * Constants::ProtonMass_SI;

    return sqrt(result);
}

double  GasComponent::SoundSpeed_CMU(const double temperature)
{
    double result = Constants::Boltzman_CMU * temperature;
    result /= meanMolecularWeight * Constants::ProtonMass_CMU;

    return sqrt(result);
}

double  GasComponent::Temperature_SI(const double mC, const double r)
{
    /*
     *   cT = 1,034.15396587001420313006 / 1.3806488 * 10 ^ (-4 - 11 + 30 - 27 + 23) = 7.4903477688896278556144038947486e13
     */
    static double cTp = Constants::NewtonG * Constants::ProtonMassBoltzman_SI;
    static double pT = 2.0 * scaleHeight.index - 3.0;
    double cT = SQR(scaleHeight.c) * mC * meanMolecularWeight * cTp;
    cT *= pow(Constants::AuToMeter, -1.0 - pT);
    double result = cT * pow(r, pT);

    return result;
}

double  GasComponent::Temperature_CMU(const double mC, const double r)
{
    static double cTp = Constants::Gauss2 * Constants::ProtonMassBoltzman_CMU;
    static double pT = 2.0 * scaleHeight.index - 3.0;
    double cT = SQR(scaleHeight.c) * mC * meanMolecularWeight * cTp;
    double result = cT * pow(r, pT);

    return result;
}

double  GasComponent::MeanThermalSpeed_SI(const double mC, const double r)
{
	static double Cvth = sqrt((8.0 * Constants::Boltzman_SI)/(Constants::Pi * meanMolecularWeight * Constants::ProtonMass_SI));
	double result = Cvth * sqrt(Temperature_SI(r, mC));
	return result;
}

double  GasComponent::MeanThermalSpeed_CMU(const double mC, const double r)
{
	static double Cvth = sqrt((8.0 * Constants::Boltzman_CMU)/(Constants::Pi * meanMolecularWeight * Constants::ProtonMass_CMU));
	double result = Cvth * sqrt(Temperature_CMU(r, mC));
	return result;
}
