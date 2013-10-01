#ifndef EPHEMERIS_H_
#define EPHEMERIS_H_

#include "Phase.h"
#include "OrbitalElement.h"
#include "DateTime.h"

enum EpochFormat
{
	UndefinedFormat,
    JulianDate,
    Date,
    PackedDate
};

class Ephemeris {

public:
	static int CalculateOrbitalElement(const double mu, const Phase *phase, OrbitalElement *orbitalElement);
	static int CalculateOrbitalElement(const double mu, const Phase *phase, double *a, double *e);
	static int CalculatePhase(const double mu, const OrbitalElement *oe, Phase *phase);

	static double CalculateEnergy(const double mu, Phase phase);
	static double CalculateKineticEnergy(Vector v);
    static double CalculatePotentialEnergy(const double mu, Vector r);
	static int KeplerEquationSolver(const double e, const double m, const double eps, double &E);

	static void ShiftIntoRange(const double lower, const double upper, double &value);

	static int GetJulianDate(std::string epoch, double &julianDate);
	static EpochFormat GetFormat(std::string epoch);
	static void DateFormat(std::string epoch, DateTime &dateTime);
	static void PackedDateFormat(std::string epoch, DateTime &dateTime);

	static void ToJulianDate(DateTime epoch, double& julianDate);
	static void ToDateTime(double julianDate, DateTime& dateTime);
	static int MyInt(double x);
	static double MyFrac(double x);
};

#endif
