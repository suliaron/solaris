#include <algorithm>
#include <string>

#include "Units.h"
#include "Constants.h"
#include "Error.h"

int UnitTool::MassToSolar(std::string unit, double& value)
{
	std::transform(unit.begin(), unit.end(), unit.begin(), ::tolower);	
	if (		unit == "solar")
		value *= 1;
	else if (	unit == "mercury" ) {
		value *= Constants::MercuryToSolar;
	}
	else if (	unit == "venus" ) {
		value *= Constants::VenusToSolar;
	}
	else if (	unit == "earth" ) {
		value *= Constants::EarthToSolar;
	}
	else if (	unit == "earthmoon" ) {
		value *= Constants::EarthMoonToSolar;
	}
	else if (	unit == "mars" ) {
		value *= Constants::MarsToSolar;
	}
	else if (	unit == "jupiter" ) {
		value *= Constants::JupiterToSolar;
	}
	else if (	unit == "saturn" ) {
		value *= Constants::SaturnToSolar;
	}
	else if (	unit == "uranus" ) {
		value *= Constants::UranusToSolar;
	}
	else if (	unit == "neptune" ) {
		value *= Constants::NeptuneToSolar;
	}
	else if (	unit == "kilogram" || unit == "kg") {
		value *= Constants::KilogramToSolar;
	}
	else if (	unit == "gram" || unit == "g") {
		value *= Constants::GramToSolar;
	}
	else {
		Error::_errMsg = "Undefined mass unit!";
		Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
		return 1;
	}

	return 0;
}

int UnitTool::TimeToDay(std::string unit, double& value)
{
	std::transform(unit.begin(), unit.end(), unit.begin(), ::tolower);	
	if (		unit == "day") {
		value *= 1;
	}
	else if (	unit == "year" ) {
		value *= Constants::YearToDay;
	}
	else if (	unit == "second" ) {
		value *= Constants::SecondToDay;
	}
	else {
		Error::_errMsg = "Undefined time unit!";
		Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
		return 1;
	}

	return 0;
}

int UnitTool::DistanceToAu(std::string unit, double& value)
{
	std::transform(unit.begin(), unit.end(), unit.begin(), ::tolower);	
	if (		unit == "au") {
		value *= 1;
	}
	else if (	unit == "m" || unit == "meter" ) {
		value *= Constants::MeterToAu;
	}
	else if (	unit == "km" || unit == "kilometer") {
		value *= Constants::KilometerToAu;
	}
	else if (	unit == "solarradius") {
		value *= Constants::SolarRadiusToAu;
	}
	else {
		Error::_errMsg = "Undefined distance unit!";
		Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
		return 1;
	}

	return 0;
}

int UnitTool::AngleToRadian(std::string unit, double& value)
{
	std::transform(unit.begin(), unit.end(), unit.begin(), ::tolower);	
	if (		unit == "radian") {
		value *= 1;
	}
	else if (	unit == "degree" ) {
		value *= Constants::DegreeToRadian;
	}
	else {
		Error::_errMsg = "Undefined angle unit!";
		Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
		return 1;
	}

	return 0;
}

int UnitTool::VelocityToCM(std::string unit, double& value)
{
	std::transform(unit.begin(), unit.end(), unit.begin(), ::tolower);	
	if (		unit == "auday") {
		value *= 1;
	}
	else {
		Error::_errMsg = "Undefined velocity unit!";
		Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
		return 1;
	}

	return 0;
}

int UnitTool::SurfaceDensityToCM(std::string unit, double& value)
{
	std::transform(unit.begin(), unit.end(), unit.begin(), ::tolower);	
	if (		unit == "solarau2") {
		value *= 1;
	}
	else if (	unit == "earthau2" ) {
		value *= Constants::EarthPerAu2ToSolarPerAu2;
	}
	else if (	unit == "gcm2" ) {
		value *= Constants::GramPerCm2ToSolarPerAu2;
	}
	else {
		Error::_errMsg = "Undefined surface density unit!";
		Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
		return 1;
	}

	return 0;
}

int UnitTool::VolumeDensityToCM(std::string unit, double& value)
{
	std::transform(unit.begin(), unit.end(), unit.begin(), ::tolower);	
	if (		unit == "solarau3") {
		value *= 1;
	}
	else if (	unit == "gcm3" ) {
		value *= Constants::GramPerCm3ToSolarPerAu3;
	}
	else {
		Error::_errMsg = "Undefined volume density unit!";
		Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
		return 1;
	}

	return 0;
}
