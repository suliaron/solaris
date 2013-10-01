#ifndef ORBITAL_ELEMENT_H_
#define ORBITAL_ELEMENT_H_

#include <iostream>

class OrbitalElement {

public:

    double semiMajorAxis;
    double eccentricity;
    double inclination;
    double argumentOfPericenter;
    double longitudeOfNode;
    double meanAnomaly;
    double meanMotion;

	OrbitalElement();
    OrbitalElement(double sma, double ecc, double inc, double peri, double node, double mean);
	OrbitalElement(const OrbitalElement &orbitalElement);

	double CalculateLongitudeOfPericenter();

	// Input/Output streams
	friend std::istream& operator>>( std::istream&, OrbitalElement&);
	friend std::ostream& operator<<( std::ostream&, OrbitalElement);

private:
	std::string _errMsg;
};

#endif