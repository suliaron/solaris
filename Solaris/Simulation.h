#ifndef SIMULATION_H_
#define SIMULATION_H_

#include <string>
#include "BodyData.h"
#include "BodyGroupList.h"
#include "Settings.h"

class	BinaryFileAdapter;
class	Body;
class	FargoParameters;
class	Nebula;

class Simulation {
public:
	
	Simulation();
	Simulation(std::string runType);

	int Initialize();
	int SetStartTimeOfMainPhase();
	int SetTimeAndSaveOfTimeLine();
	int InitializeTimeLineAndBodyGroups();
	int InitializePhases();

	void SortBodyGroupsByStartTime();

	int CheckBodyGroupList();
	int CheckStartTimes();

	//int CalculateStartTime(double &start);

	Body* FindBy(int bodyId);
	Body* FindCentralBody();

	BodyGroupList		bodyGroupList;
	std::list<Body *>	bodyList;
	BodyData			bodyData;

	Settings			settings;
	Nebula				*nebula;
	BinaryFileAdapter	*binary;
	FargoParameters		*fargoParameters;

	std::string			name;
	std::string			description;
	std::string			referenceFrame;
    std::string         runType;
};

#endif
