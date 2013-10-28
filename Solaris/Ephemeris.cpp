#include <algorithm>
#include <string>
#include <cctype>

#include "Ephemeris.h"
#include "Error.h"
#include "DateTime.h"
#include "Constants.h"

int Ephemeris::CalculateOrbitalElement(const double mu, const Phase *phase, double *a, double *e)
{
    const double sq3 = 1.0e-14;

	// Calculate energy, h
    double h = CalculateEnergy(mu, *phase);
    if (h >= 0.0)
    {
    	Error::_errMsg = "The energy is positive or zero!";
    	Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
        return 1;
    }
    Vector r = phase->position;
    Vector v = phase->velocity;

    Vector c = Vector::CrossProduct(r, v);
    /*
    * Calculate eccentricity, e
    */
    double e2 = 1.0 + 2.0 * c.LengthSQR() * h / (mu * mu);
    if (abs(e2) < sq3)
    {
        e2 = 0.0;
    }
    *e = sqrt(e2);
    /*
    * Calculate semi-major axis, a
    */
    *a = -mu / (2.0 * h);

    return 0;
}

int Ephemeris::CalculateOrbitalElement(const double mu, const Phase *phase, OrbitalElement *orbitalElement)
{
    const double sq2 = 1.0e-14;
    const double sq3 = 1.0e-14;

	// Calculate energy, h
    double h = CalculateEnergy(mu, *phase);
    if (h >= 0.0)
    {
    	Error::_errMsg = "The energy is positive or zero";
    	Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
        return 1;
    }
    Vector r = phase->position;
    Vector v = phase->velocity;

    Vector c = Vector::CrossProduct(r, v);
	Vector l = (-mu / r.Length()) * (r) + Vector::CrossProduct(v, c);
	double rv = Vector::DotProduct(r, v);
    /*
    * Calculate eccentricity, e
    */
    double e2 = 1.0 + 2.0 * c.LengthSQR() * h / (mu * mu);
    if (abs(e2) < sq3)
    {
        e2 = 0.0;
    }
    double e = sqrt(e2);
    /*
    * Calculate semi-major axis, a
    */
    double a = -mu / (2.0 * h);
    /*
    * Calculate inclination, incl
    */
    double cosi = c.z / c.Length();
    double sini = sqrt(c.x * c.x + c.y * c.y) / c.Length();
    double incl = acos(cosi);
    if (incl < sq2)
    {
        incl = 0.0;
    }
    /*
    * Calculate longitude of node, O
    */
    double node = 0.0;
    if (incl != 0.0)
    {
        double tmpx = -c.y / (c.Length() * sini);
        double tmpy = c.x / (c.Length() * sini);
		node = atan2(tmpy, tmpx);
		ShiftIntoRange(0.0, 2.0*Constants::Pi, node);
    }
    /*
    * Calculate argument of pericenter, w
    */
    double E = 0.0;
    double peri = 0.0;
    if (e2 != 0.0)
    {
        double tmpx = (l.x * cos(node) + l.y * sin(node)) / l.Length();
        double tmpy = (-l.x * sin(node) + l.y * cos(node)) / (l.Length() * cosi);
        peri = atan2(tmpy, tmpx);
        ShiftIntoRange(0.0, 2.0*Constants::Pi, peri);

        tmpx = 1.0 / e * (1.0 - r.Length() / a);
        tmpy = rv / (sqrt(mu * a) * e);
        E = atan2(tmpy, tmpx);
        ShiftIntoRange(0.0, 2.0*Constants::Pi, E);
    }
    else
    {
        peri = 0.0;
        E = atan2(r.y, r.x);
        ShiftIntoRange(0, 2.0*Constants::Pi, E);
    }
    /*
    * Calculate mean anomaly, M
    */
    double M = E - e * sin(E);
    ShiftIntoRange(0, 2.0*Constants::Pi, M);

	orbitalElement->semiMajorAxis			= a;
	orbitalElement->eccentricity			= e;
	orbitalElement->inclination				= incl;
	orbitalElement->argumentOfPericenter	= peri;
	orbitalElement->longitudeOfNode			= node;
	orbitalElement->meanAnomaly				= M;

	return 0;
}

/// <summary>
/// Computes the phase (position and velocity) of the body from the keplerian orbital elements.
/// </summary>
/// <param name="mu">(Gauss constant)^2*(m1 + m2)</param>
/// <param name="oe">The Keplerian orbital element</param>
/// <returns>The phase of the body</returns>
int Ephemeris::CalculatePhase(const double mu, const OrbitalElement* oe, Phase* phase)
{
    double e = oe->eccentricity;
	double E = 0;
	if (KeplerEquationSolver(e, oe->meanAnomaly, 1.0e-14, E) == 1) {
		Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
		return 1;
	}
    double v = 2.0 * atan(sqrt((1.0 + e) / (1.0 - e)) * tan(E / 2.0));

	double a = oe->semiMajorAxis;
    double p = a * (1.0 - e * e);
    double r = p / (1.0 + e * cos(v));
    double kszi = r * cos(v);
    double eta = r * sin(v);
    double vKszi = -sqrt(mu / p) * sin(v);
    double vEta = sqrt(mu / p) * (e + cos(v));

    double peri = oe->argumentOfPericenter;
    double node = oe->longitudeOfNode;
    double incl = oe->inclination;
    double cw = cos(peri);
    double sw = sin(peri);
    double cO = cos(node);
    double sO = sin(node);
    double ci = cos(incl);
    double si = sin(incl);

    Vector P = Vector(cw * cO - sw * sO * ci, cw * sO + sw * cO * ci, sw * si);
    Vector Q = Vector(-sw * cO - cw * sO * ci, -sw * sO + cw * cO * ci, cw * si);

	phase->position = Vector(kszi * P + eta * Q);
	phase->velocity = Vector(vKszi * P + vEta * Q);

	return 0;
}

/// <summary>
/// Solves the Kepler-equation using the standard procedure from the book of
/// B. Erdi page 53. If the error does not reach eps in 25 iterations,
/// the procedure stops and returns 1.
/// </summary>
/// <param name="e">The eccentricity</param>
/// <param name="m">The mean anomaly</param>
/// <param name="eps">The required precision of the solution</param>
/// <returns>Returns the eccentric anomaly</returns>
int Ephemeris::KeplerEquationSolver(const double e, const double m, const double eps, double &E)
{
	if (e == 0.0 || m == 0.0 || m == Constants::Pi)
    {
        E = m;
		return 0;
    }
    E = m + e * (sin(m)) / (1.0 - sin(m + e) + sin(m));
    double E1 = 0.0;
    double error;
    int step = 0;
    do {
        E1 = E - (E - e * sin(E) - m) / (1.0 - e * cos(E));
        error = abs(E1 - E);
        E = E1;
        step++;
    } while (error > eps && step <= 25);
	if (step > 25 ) {
		Error::_errMsg = "Could not compute the excentric anomaly E!";
		Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
		return 1;
	}

	return 0;
}

inline double Ephemeris::CalculateKineticEnergy(Vector v)
{
	return v.LengthSQR() / 2.0;
}

inline double Ephemeris::CalculatePotentialEnergy(const double mu, Vector r)
{
    return -mu / r.Length();
}

inline double Ephemeris::CalculateEnergy(const double mu, Phase phase)
{
	return CalculateKineticEnergy(phase.velocity) + CalculatePotentialEnergy(mu, phase.position);
}

void Ephemeris::ShiftIntoRange(const double lower, const double upper, double &value)
{
    double range = upper - lower;
    while (value >= upper)
    {
        value -= range;
    }
    while (value < lower)
    {
        value += range;
    }
}

/// <summary>
/// Calculates the Julian date from the epoch parameter. The format of
/// the epoch can be JulianDate, Date or Packed date format.
/// </summary>
/// <param name="epoch">The time instance in JulianDate, Date or Packed date format</param>
/// <returns>The time instance in double</returns>
int Ephemeris::GetJulianDate(std::string epoch, double &julianDate)
{
	DateTime dateTime;

	std::string data = epoch.substr(2, epoch.length() - 2);
	EpochFormat epochFormat = GetFormat(epoch);
	switch (epochFormat)
	{
		case JulianDate:
			julianDate = atof(data.c_str());
			break;
		case Date:
			DateFormat(epoch, dateTime);
			ToJulianDate(dateTime, julianDate);
			break;
		case PackedDate:
			PackedDateFormat(epoch, dateTime);
			ToJulianDate(dateTime, julianDate);
			break;
		case UndefinedFormat:
			Error::_errMsg = "Undefined format for epoch!";
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			return 1;
		default:
			Error::_errMsg = "Invalid value for epoch!";
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			return 1;
	}

	return 0;
}

/// <summary>
/// Determines the format of the epoch parameter. If the format is bad it
/// prints out an error message and returns with an UndefinedFormat.
/// </summary>
/// <param name="epoch">The time instance</param>
/// <returns>The format of the epoch</returns>
EpochFormat Ephemeris::GetFormat(std::string epoch)
{
    // Is epoch julian date (begins with JD)
	std::string prefix = epoch.substr(0, 2);
	std::transform(prefix.begin(), prefix.end(), prefix.begin(), ::tolower);
    if (prefix == "jd")
    {
		for (size_t i = 2; i < epoch.length(); i++) {
			char c = epoch[i];
			if (isdigit(c) || c == '.')
            {
                continue;
            }
			else {
				std::cout << "Bad JulianDate format: " << epoch << std::endl;
				return UndefinedFormat;
			}
		}
        return JulianDate;
    }

    // Is epoch date (begins with a number)
    if (isdigit(epoch[0]))
    {
        if (epoch.length() > 8)
        {
			std::cout << "Bad Date format: " << epoch << std::endl;
			return UndefinedFormat;
        }
		for (size_t i = 0; i < epoch.length(); i++) {
			char c = epoch[i];
			if (isdigit(c))
            {
                continue;
            }
			else {
				std::cout << "Bad Date format: " << epoch << std::endl;
				return UndefinedFormat;
			}
		}
        return Date;
    }

    // Is epoch packed date (begins with a Letter)
    if (epoch.length() > 5)
    {
		std::cout << "Bad PackedDate format: " << epoch << ". The length of the epoch in packed format must be 5 character long!" << std::endl;
		return UndefinedFormat;
    }
    char c = epoch[0];
    if (c < 'I' || c > 'K')
    {
		std::cout << "Bad PackedDate format: " << epoch << ". The first character of the epoch in packed format must be I, J or K!" << std::endl;
		return UndefinedFormat;
    }            
    if (!isdigit(epoch[1]) || !isdigit(epoch[2]))
    {
		std::cout << "Bad PackedDate format: " << epoch << ". The second and third character of the epoch in packed format must be a digit!" << std::endl;
		return UndefinedFormat;
    }
    c = epoch[3];
    bool b = isdigit(c) || (c >= 'A' && c <= 'C');
    if (!b)
    {
		std::cout << "Bad PackedDate format: " << epoch << ". The fourth character of the epoch in packed format must be a digit or A, B or C!" << std::endl;
		return UndefinedFormat;
    }
	c = epoch[4];
	b = isdigit(c) || (c >= 'A' && c <= 'V');
    if (!b)
    {
		std::cout << "Bad PackedDate format: " << epoch << ". The last character of the epoch in packed format must be a digit or a letter from A to V!" << std::endl;
		return UndefinedFormat;
    }
    return PackedDate;
}

void Ephemeris::DateFormat(std::string epoch, DateTime &dateTime)
{
    std::string data = epoch.substr(0, 4);
	dateTime.Year =  atoi(data.c_str());
    data = epoch.substr(4, 2);
	dateTime.Month = atoi(data.c_str());
    data = epoch.substr(6, 2);
	dateTime.Day = atoi(data.c_str());
}

void Ephemeris::PackedDateFormat(std::string epoch, DateTime &dateTime)
{
	// Get the first two digits of the year:
	char c = epoch[0];
	int year = (c - 'I' + 18) * 100;
	std::string s = epoch.substr(1, 2);
	dateTime.Year = year + atoi(s.c_str());

	// Get the month of the epoch
	c = epoch[3];
	dateTime.Month = isdigit(c) ? c - '0' : c - 'A' + 10;

	// Get the day of the epoch
	c = epoch[4];
	dateTime.Day = isdigit(c) ? c - '0' : c - 'A' + 10;
}

void Ephemeris::ToJulianDate(DateTime epoch, double& julianDate)
{
	int y = epoch.Month <= 2 ? epoch.Year - 1 : epoch.Year;
	int m = epoch.Month <= 2 ? epoch.Month + 12 : epoch.Month;
	/* 
	* The reason for the peculiar definition of B are to account for the "lost days" in October 1582
	* when the Gregorian Calendar replaced the Julian Calendar in Europe, and to deal with the
	* introduction of a leap day in the Gregorian calendar.
	*/

	int B;
	if (epoch.Year < 1582)
	{
		B = -2;
	}
	else if (epoch.Year == 1582 && epoch.Month < 10)
	{
		B = -2;
	}
	else if (epoch.Year == 1582 && epoch.Month == 10 && epoch.Day <= 4)
	{
		B = -2;
	}
	else
	{
		B = MyInt(y / 400.0) - MyInt(y / 100.0);
	}

	double universalTime = epoch.Millisec / 3600000.0 + epoch.Second / 3600.0 + epoch.Minute / 60.0 + epoch.Hour;
	julianDate = MyInt(365.25 * y) + MyInt(30.6001 * (m + 1)) + B + 1720996.5 + epoch.Day + universalTime / 24.0;
}

void Ephemeris::ToDateTime(double julianDate, DateTime& dateTime)
{
    long a = MyInt(julianDate + 0.5);
    long b = 0;
    long c = 0;
    if (a < 2299161)
    {
        b = 0;
        c = a + 1524;
    }
    if (a >= 2299161)
    {
        b = MyInt((a - 1867216.25) / 36524.25);
        c = a + b - MyInt(b / 4.0) + 1525;
    }
    long d = MyInt((c - 122.1) / 365.25);
    long e = MyInt(365.25 * d);
    long f = MyInt((c - e) / 30.6001);

    double day = c - e - MyInt(30.6001 * f) + MyFrac(julianDate + 0.5);
    int month = (f - 1 - 12 * MyInt(f / 14.0));
    int year = (d - 4715 - MyInt((7 + month) / 10.0));

    double universalTime = MyFrac(day) * 24.0;
    int hour = MyInt(universalTime);
    double r = MyFrac(universalTime) * 60.0;
    int minute = MyInt(r);
    r = MyFrac(r) * 60;
    int second = MyInt(r);

	dateTime.Year = year;
	dateTime.Month = month;
	dateTime.Day = MyInt(day);
	dateTime.Hour = hour;
	dateTime.Minute = minute;
	dateTime.Second = second;
}

int Ephemeris::MyInt(double x)
{
    return (int)floor(x);
}

double Ephemeris::MyFrac(double x)
{
    double fraction;
    if (x > 0)
    {
        fraction = x - floor(x);
    }
    else
    {
        fraction = x - ceil(x);
    }
    return fraction;
}

