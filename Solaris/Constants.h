#ifndef CONSTANTS_H_
#define CONSTANTS_H_

#include <string>

#define SQR(a)		((a)*(a))
#define CUBE(a)		((a)*(a)*(a))

namespace Constants
{
	const std::string CodeName		      = "Solaris";
	const std::string Version		      = "1.0";

	const std::string Usage				  = "Usage is -i|-c <infile>\n";

	const int	 CheckForSM			      = 100;
	const double SmallestNumber		      = 1.0e-50;

	const double Pi					      = 3.14159265358979323846;
	const double SqrtTwoPi			      = 2.50662827463100024161;

	const double Boltzman_SI		      = 1.3806488e-23;			        // J/K = (kg m^2) / (s^2 K)
	const double ProtonMass_SI		      = 1.672621777e-27;		        // kg
    const double BoltzmanProtonMass_SI    = Boltzman_SI / ProtonMass_SI;    // m^2 / (s^2 K) 8.25439928491377e3
    const double ProtonMassBoltzman_SI    = 1.0 / BoltzmanProtonMass_SI;    // (s^2 K) / m^2

    const double NewtonG			      = 6.67384e-11;		        	// m^3 / (kg s^2) = N m^2 / kg^2
	const double Gauss				      = 1.720209895e-2;                 // Au^3/2 / (Solar^1/2 day)
	const double Gauss2				      = 2.959122082855911025e-4;        // Au^3 / (Solar day^2)

	const double DegreeToRadian		      = 1.745329251994329509e-2;
	const double RadianToDegree		      = 1.0 / DegreeToRadian;

	const double SolarToMercury		      = 6.023600e6;
	const double SolarToVenus		      = 4.0852371e5;
	const double SolarToEarth		      = 3.3294605e5;
	const double SolarToEarthMoon	      = 3.2890056e5;
	const double SolarToMars		      = 3.098708e6;
	const double SolarToJupiter		      = 1.0473486e3;
	const double SolarToSaturn	 	      = 3.497898e3;
	const double SolarToUranus		      = 2.290298e4;
	const double SolarToNeptune		      = 1.941224e4;
	const double EarthToMoon		      = 8.130059e1;

	const double SolarToKilogram	      = 1.98911e30;
	const double SolarToGram		      = 1.0e3 * SolarToKilogram;

	const double MercuryToSolar		      = 1.0 / SolarToMercury;
	const double VenusToSolar		      = 1.0 / SolarToVenus;
	const double EarthToSolar		      = 1.0 / SolarToEarth;
	const double EarthMoonToSolar	      = 1.0 / SolarToEarthMoon; 
	const double MarsToSolar		      = 1.0 / SolarToMars;
	const double JupiterToSolar		      = 1.0 / SolarToJupiter;
	const double SaturnToSolar		      = 1.0 / SolarToSaturn;
	const double UranusToSolar		      = 1.0 / SolarToUranus;
	const double NeptuneToSolar		      = 1.0 / SolarToNeptune;
	const double MoonToEarth		      = 1.0 / EarthToMoon;

	const double KilogramToSolar	      = 1.0 / SolarToKilogram;
	const double GramToSolar		      = 1.0 / SolarToGram;

	const double AuToMeter			      = 1.495978707e11;
	const double AuToCentimeter		      = AuToMeter * 1.0e2;
	const double AuToKilometer		      = AuToMeter / 1.0e3;
	const double AuToSolarRadius	      = 215.094;

	const double MeterToAu			      = 1.0 / AuToMeter;
	const double KilometerToAu		      = 1.0 / AuToKilometer;
	const double SolarRadiusToAu	      = 1.0 / AuToSolarRadius;

	const double YearToDay			      = 365.25;
	const double DayToSecond		      = 86400.0;

	const double DayToYear			      = 1.0 / YearToDay;
	const double SecondToDay		      = 1.0 / DayToSecond;

	const double GramPerCm2ToEarthPerAu2  = GramToSolar*SolarToEarth / (SQR(1.0e-2*MeterToAu));
	const double EarthPerAu2ToGramPerCm2  = 1.0 / GramPerCm2ToEarthPerAu2;

	const double GramPerCm2ToSolarPerAu2  = GramToSolar / (SQR(1.0e-2*MeterToAu));
	const double SolarPerAu2ToGramPerCm2  = 1.0 / GramPerCm2ToSolarPerAu2;
	const double EarthPerAu2ToSolarPerAu2 = EarthToSolar;

	const double GramPerCm3ToSolarPerAu3  = GramToSolar / (CUBE(1.0e-2*MeterToAu));
	const double SolarPerAu3ToGramPerCm3  = 1.0 / GramPerCm3ToSolarPerAu3;

	const double Boltzman_CMU		      = Boltzman_SI * (KilogramToSolar * SQR(MeterToAu)) / (SQR(SecondToDay));	         // (Solar AU^2)/(day^2 K) 2.315266990160467e-66
	const double ProtonMass_CMU		      = ProtonMass_SI * KilogramToSolar; // Solar       8.408895320017495e-58 
    const double BoltzmanProtonMass_CMU   = Boltzman_CMU / ProtonMass_CMU;   // Au^2 / (d^2 K) 2.753354515722108e-9
    const double ProtonMassBoltzman_CMU   = 1.0 / BoltzmanProtonMass_CMU;    // (d^2 K) / Au^2  

}

#endif
