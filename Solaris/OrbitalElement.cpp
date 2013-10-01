#include "OrbitalElement.h"
#include "Constants.h"
#include "Ephemeris.h"

OrbitalElement::OrbitalElement()
{
	semiMajorAxis = 0.0;
	eccentricity = 0.0;
	inclination = 0.0;
	argumentOfPericenter = 0.0;
	longitudeOfNode = 0.0;
	meanAnomaly = 0.0;
	meanMotion = 0.0;
}

OrbitalElement::OrbitalElement(double sma, double ecc, double inc, double peri, double node, double mean)
{
    semiMajorAxis = sma;
    eccentricity = ecc;
    inclination = inc;
    argumentOfPericenter = peri;
    longitudeOfNode = node;
    meanAnomaly = mean;
	meanMotion = 0.0;
}

OrbitalElement::OrbitalElement(const OrbitalElement &orbitalElement) {
	semiMajorAxis = orbitalElement.semiMajorAxis;
	eccentricity = orbitalElement.eccentricity;
	inclination = orbitalElement.inclination;
	argumentOfPericenter = orbitalElement.argumentOfPericenter;
	longitudeOfNode = orbitalElement.longitudeOfNode;
	meanAnomaly = orbitalElement.meanAnomaly;
	meanMotion = orbitalElement.meanMotion;
}

double OrbitalElement::CalculateLongitudeOfPericenter()
{
    double longitudeOfPericenter = argumentOfPericenter + longitudeOfNode;
	Ephemeris::ShiftIntoRange(0, Constants::Pi, longitudeOfPericenter);
    return longitudeOfPericenter;
}

// Read input OrbitalElement format: "(%f, %f, %f, %f, %f, %f)"
std::istream& operator>>(std::istream& input, OrbitalElement& o)
{
	char c;
    // skip '('
	input >> c;

	input >> o.semiMajorAxis;
	input >> c;		// skip ','

	input >> o.eccentricity;
	input >> c;		// skip ','

	input >> o.inclination;
	input >> c;		// skip ','

	input >> o.argumentOfPericenter;
	input >> c;		// skip ','

	input >> o.longitudeOfNode;
	input >> c;		// skip ','

	input >> o.meanAnomaly;
	input >> c;		// skip ','

	// skip ')'
	input >> c;

	return input;
}

// Write output OrbitalElement in format: "(%f, %f, %f, %f, %f, %f)"
std::ostream& operator<<(std::ostream& output, OrbitalElement o)
{
	output << "(" << o.semiMajorAxis << ", " << o.eccentricity << ", " << o.inclination << ", " << o.argumentOfPericenter << ", " << o.longitudeOfNode << ", " << o.meanAnomaly << ")";
	return output;
}

