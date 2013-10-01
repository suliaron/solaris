#ifndef RUNGEKUTTA4_H_
#define RUNGEKUTTA4_H_

#include <string>

class Acceleration;
class BodyData;
class TimeLine;

class RungeKutta4
{
public:
	RungeKutta4();

	int			Driver(BodyData *bodyData, Acceleration *acceleration, TimeLine *timeLine);
	int 		Step(  BodyData *bodyData, Acceleration *acceleration);

	std::string	name;
	std::string	reference;
	double		accuracy;
	double		epsilon;
};

#endif