#ifndef SIMULATOR_H_
#define SIMULATOR_H_

#include <list>

#include "Counter.h"
#include "IntegratorType.h"
#include "BinaryFileAdapter.h"
#include "BodyData.h"
#include "BodyGroup.h"
#include "Event.h"

class Acceleration;
class RungeKutta4;
class RungeKuttaFehlberg78;
class DormandPrince;
class TimeLine;
class Simulation;

class Simulator
{
public:

	Simulator(Simulation* simulation);

    int     Continue();
	int		Run();

	IntegratorType			integratorType;
	RungeKuttaFehlberg78	*rungeKuttaFehlberg78;
	RungeKutta4				*rungeKutta4;
	DormandPrince			*dormandPrince;

	BodyData				bodyData;
	Counter					counter;
	BinaryFileAdapter::OutputType outputType;

private:
	int 	Synchronization();
	int		PreIntegration();
	int		MainIntegration();
	int		Integrate(TimeLine* timeLine);
	int		DecisionMaking(TimeLine* timeLine, bool& stop);
	int		DecisionMaking(const long int stepCounter, TimeLine* timeLine, double* hSum, bool& stop);

	int		Insert(double time, std::list<BodyGroup>::iterator &bgIt);
	int		Insert(double time, std::list<BodyGroup *>::iterator &bgIt);
	int		SetIteratorAfter(BodyType type, std::list<Body *>::iterator& it);
	int		PopulateBodyList(double time);

	double	ShortestPeriod();

	int 	BodyListToBodyData();
	void	UpdateBodyListAfterIntegration();

	int 	CheckEvent(double timeOfEvent);
	int		RemoveBody(int bodyId);
	int		HandleCollision(int idx1, int idx2, int& survivIdx, int &mergerIdx, int& survivId, int& mergerId);
	int		CalculatePhaseAfterCollision(int survivIdx, int mergerIdx);
	int		CalculateCharacteristicsAfterCollision(int survivId, int mergerId, int survivIdx, int mergerIdx);
	int		CalculateComponentListUnion(int survivId, int mergerId, std::list<Component> **survivList, std::list<Component> **mergerList, std::list<Component> *result);

	Event			_ejectionEvent;
	Event			_hitCentrumEvent;
	Event			_closeEncounterEvent;
	Event			_collisionEvent;
	Event			_weakCaptureEvent;

	time_t			_startTime;
	Simulation*		_simulation;
	Acceleration*	_acceleration;
};

#endif
