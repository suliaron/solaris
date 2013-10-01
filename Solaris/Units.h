#ifndef UNITS_H_
#define UNITS_H_

enum DistanceUnit
{
    Au,
    m,
    km,
    SolarRadius
};

enum AngleUnit
{
    Degree,
    Radian
};

enum MassUnit
 {
    Solar,
    Mercury,
    Venus,
    Earth,
    Earthmoon,
    Mars,
    Jupiter,
    Saturn,
    Uranus,
    Neptune,
    Kilogram,
	Gram
};

enum TimeUnit
{
    Day,
    Year,
    Second
};

enum VelocityUnit
{
    AuDay
};

enum SurfaceDensityUnit
{
    SolarAu2,
    EarthAu2,
    GCm2
};

enum VolumeDensityUnit
{
    SolarAu3,
    GCm3
};


class UnitTool {
public:
	static int MassToSolar(std::string unit, double& value);
	static int TimeToDay(std::string unit, double& value);
	static int DistanceToAu(std::string unit, double& value);
	static int AngleToRadian(std::string unit, double& value);
	static int VelocityToCM(std::string unit, double& value);
	static int SurfaceDensityToCM(std::string unit, double& value);
	static int VolumeDensityToCM(std::string unit, double& value);
};

#endif