#ifndef SETTINGS_H_
#define SETTINGS_H_

class EventCondition;
class Integrator;
class Output;
class TimeLine;

class Settings {
public:

	Settings();

	bool			enableDistinctStartTimes;
	bool			baryCentric;
	Output			*output;
	Integrator		*integrator;
	TimeLine		*timeLine;

	double			ejection;
	double			hitCentrum;
	EventCondition	*closeEncounter;
	EventCondition	*collision;
	EventCondition	*weakCapture;
};

#endif