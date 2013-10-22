#ifndef SETTINGS_H_
#define SETTINGS_H_

#include "Integrator.h"
#include "Output.h"

class EventCondition;
class TimeLine;

class Settings {
public:

	Settings();

	bool			enableDistinctStartTimes;
	bool			baryCentric;
	Integrator		integrator;
	Output			output;
	TimeLine		*timeLine;

	double			ejection;
	double			hitCentrum;
	EventCondition	*closeEncounter;
	EventCondition	*collision;
	EventCondition	*weakCapture;
};

#endif