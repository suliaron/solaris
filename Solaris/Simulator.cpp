#include <cstdio>
#include <cstring>
#include <ctime>
#include <sstream>

#include "Acceleration.h"
#include "BinaryFileAdapter.h"
#include "Body.h"
#include "BodyGroupList.h"
#include "Calculate.h"
#include "Constants.h"
#include "DormandPrince.h"
#include "Error.h"
#include "EventCondition.h"
#include "Integrator.h"
#include "IntegratorType.h"
#include "RungeKutta4.h"
#include "RungeKuttaFehlberg78.h"
#include "Settings.h"
#include "Simulation.h"
#include "Simulator.h"
#include "TimeLine.h"
#include "Tools.h"

#define SQR(a)		((a)*(a))
#define CUBE(a)		((a)*(a)*(a))

 /**
  * Used during the decision making to decide if it is time to save the data or stop the integration
  */
enum StepSums {
	ALL_STEP,		/**< the sum of all the steps since the start of the integration */
	LAST_SAVE,		/**< the sum of the steps since the last save phases */
	LAST_NSTEP		/**< the sum of the last NSTEP steps, only for debugging purposes */
};

Simulator::Simulator(Simulation *simulation)
{
	_simulation			= simulation;

	_startTime			= (time_t)0;
	_acceleration		= 0;
	rungeKuttaFehlberg78= 0;
	rungeKutta4			= 0;
	dormandPrince		= 0;

	integratorType  = RUNGE_KUTTA_FEHLBERG78;
	nAllStep		= 0;
	nFailedStep		= 0;
	nSuccededStep	= 0;
}

int Simulator::Continue()
{


    return 0;
}

int Simulator::Run()
{
	if (     _simulation->settings.integrator.name == "rungekutta78" || _simulation->settings.integrator.name == "rungekuttafehlberg78") {
		integratorType = RUNGE_KUTTA_FEHLBERG78;
		rungeKuttaFehlberg78 =  new RungeKuttaFehlberg78();
	}
	else if (_simulation->settings.integrator.name == "rungekutta4") {
		integratorType = RUNGE_KUTTA4;
		rungeKutta4 = new RungeKutta4();
	}
	else if (_simulation->settings.integrator.name == "rungekutta56") {
		integratorType = RUNGE_KUTTA56;
	}
	else if (_simulation->settings.integrator.name == "dormandprince") {
		integratorType = DORMAND_PRINCE;
		dormandPrince = new DormandPrince();
	}
	else {
		Error::_errMsg = "Unknown integrator type!";
		Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
		return 1;
	}

	_acceleration = new Acceleration(integratorType, _simulation->settings.baryCentric, &bodyData, _simulation->nebula);

	if (_simulation->bodyGroupList.nOfDistinctStartTimes > 1) {
		_simulation->binary->Log("The synchronization phase of the simulation begins", false);
		_startTime = time(0);
		if (Synchronization() == 1) {
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			return 1;
		}
		_simulation->binary->LogTimeSpan("The synchronization phase of the simulation took ", _startTime);
	}

	// The pre-integration phase is only needed if the user has defined the start attribute
	// in the TimeLine Tag. Otherwise the start time will be the last or the first epoch
	// of the BodyGroupList, and therefore at the end of the synchronization process
	// the time equals to the start time.
	if (_simulation->settings.timeLine->startTimeDefined) {
		_simulation->binary->Log("The pre-integration phase of the simulation begins", false);
		_startTime = time(0);
		if (PreIntegration() == 1) {
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			return 1;
		}
		_simulation->binary->LogTimeSpan("The pre-integration phase of the simulation took ", _startTime);
	}

	if (_simulation->settings.timeLine->length != 0.0) {
		_simulation->binary->Log("The main-integration phase of the simulation begins", false);
		_startTime = time(0);
		if (MainIntegration() == 1) {
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			return 1;
		}
		_simulation->binary->LogTimeSpan("The main-integration phase of the simulation took ", _startTime);
	}

	return 0;
}

int Simulator::Integrate(TimeLine* timeLine)
{
	long int		stepCounter = 0;
	double			hSum[3] = { 0.0, 0.0, 0.0 };

	if (BodyListToBodyData() == 1) {
		Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
		return 1;
	}
	_simulation->binary->SavePhases(timeLine->time, bodyData.nBodies.total, bodyData.y0, bodyData.id);
	Calculate::Integrals(&bodyData);
	_simulation->binary->SaveIntegrals(timeLine->time, 16, bodyData.integrals);

	bool stop = false;
	switch (integratorType) {
		case RUNGE_KUTTA_FEHLBERG78:
			for ( ; ; ) {
				if (rungeKuttaFehlberg78->Driver(&bodyData, _acceleration, timeLine) == 1) {
					Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
					return 1;
				}
				stepCounter++;
				if (DecisionMaking(stepCounter, timeLine, hSum, stop) == 1) {
					Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
					return 1;
				}
				if (stepCounter % Constants::CheckForSM == 0) {
					int	n = 6*bodyData.nBodies.total;
					Tools::CheckAgainstSmallestNumber(n, bodyData.y);
					Tools::CheckAgainstSmallestNumber(n, bodyData.y0);
                    if (_simulation->binary->FileExists("Info"))
                    {
                        std::cout << "Time: " << timeLine->time * Constants::DayToYear << std::endl;
                    }
				}
				if (stop)
					break;
			}
			break;
		case RUNGE_KUTTA4:
			for ( ; ; ) {
				if (rungeKutta4->Driver(&bodyData, _acceleration, timeLine) == 1) {
					Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
					return 1;
				}
				stepCounter++;
				if (DecisionMaking(stepCounter, timeLine, hSum, stop) == 1) {
					Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
					return 1;
				}
				if (stepCounter % Constants::CheckForSM == 0) {
					int	n = 6*bodyData.nBodies.total;
					Tools::CheckAgainstSmallestNumber(n, bodyData.y);
					Tools::CheckAgainstSmallestNumber(n, bodyData.y0);
				}
				if (stop)
					break;
			}
			break;
		case DORMAND_PRINCE:
			for ( ; ; ) {
				if (dormandPrince->Driver(&bodyData, _acceleration, timeLine) == 1) {
					Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
					return 1;
				}
				stepCounter++;
				if (DecisionMaking(stepCounter, timeLine, hSum, stop) == 1) {
					Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
					return 1;
				}
				if (stepCounter % Constants::CheckForSM == 0) {
					int	n = 6*bodyData.nBodies.total;
					Tools::CheckAgainstSmallestNumber(n, bodyData.y);
					Tools::CheckAgainstSmallestNumber(n, bodyData.y0);
				}
				if (stop)
					break;
			}
			break;
		default:
			Error::_errMsg = "The integrator type is currently not supported!";
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			return 1;
			break;
	}
	_simulation->binary->SavePhases(timeLine->time, bodyData.nBodies.total, bodyData.y0, bodyData.id);
	Calculate::Integrals(&bodyData);
	_simulation->binary->SaveIntegrals(timeLine->time, 16, bodyData.integrals);

	return 0;
}

#define NSTEP 500
int	Simulator::DecisionMaking(const long int stepCounter, TimeLine* timeLine, double* hSum, bool& stop)
{
	hSum[ALL_STEP]	+= timeLine->hDid;
	hSum[LAST_SAVE] += timeLine->hDid;
	hSum[LAST_NSTEP]+= timeLine->hDid;

#ifdef _DEBUG
	if ((stepCounter) % NSTEP == 0) {
		fprintf(stderr, "stepCounter: %9ld t: %14.10g [hDid: %10.6g hNext: %10.6g] hAvg: {%16.10e, %16.10e}, %10.6g%% ready\n",
			stepCounter, timeLine->time, timeLine->hDid, timeLine->hNext, hSum[LAST_NSTEP]/NSTEP, hSum[ALL_STEP]/stepCounter, fabs(hSum[ALL_STEP]/timeLine->length)*100.0);
		hSum[LAST_NSTEP] = 0.0;
	}
#endif

// NOTE: az alábbi lehetöségek sorrendje fontos. Elöször ez eseményeket ellenörzöm, aztán
// pedig, hogy elértük-e az integrálás végét, az adatokat el kell-e menteni.
	if (CheckEvent(timeLine->time) == 1) {
		Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
		return 1;
	}

	// TODO: Check if the number of removed bodies exceeds a threshold, than use realloc.
	// The Acceleration::rm3
	// must be set to NULL in order to force a new allocation
	// in the subsequent call to Acceleration::Compute().
//	if (bodyData.nBodies.removed > 0 && bodyData.nBodies.removed % bodyData.nBodies.threshold == 0) {
//#ifdef _DEBUG
//		fprintf(stderr, "File: %40s, Function: %40s, Line: %10d (removed == threshold)\n", __FILE__, __FUNCTION__, __LINE__);
//#endif
//		delete[] acceleration.rm3;
//		delete[] acceleration.accelGasDrag;
//		delete[] acceleration.accelMigrationTypeI;
//		delete[] acceleration.accelMigrationTypeII;
//	}

	if (fabs(hSum[ALL_STEP]) >= fabs(timeLine->length)) {
		UpdateBodyListAfterIntegration();
		stop = true;
		return 0;
	}

	if (bodyData.nBodies.total <= 1) {
		_simulation->binary->Log("The total number of bodies reduced to 1!", true);
		UpdateBodyListAfterIntegration();
		stop = true;
		return 0;
	}

	if (fabs(hSum[ALL_STEP] + timeLine->hNext) > fabs(timeLine->length)) {
		timeLine->hNext = timeLine->length - hSum[ALL_STEP];
	}

	if (fabs(hSum[LAST_SAVE]) >= fabs(timeLine->output)) {

		_simulation->binary->SavePhases(timeLine->time, bodyData.nBodies.total, bodyData.y0, bodyData.id);
		Calculate::Integrals(&bodyData);
		_simulation->binary->SaveIntegrals(timeLine->time, 16, bodyData.integrals);

		timeLine->save += timeLine->output;
		hSum[LAST_SAVE] = 0.0;
	}

	if (fabs(hSum[LAST_SAVE] + timeLine->hNext) > fabs(timeLine->output)) {
		timeLine->hNext = timeLine->output - hSum[LAST_SAVE];
	}

	return 0;
}
#undef NSTEP

int Simulator::Synchronization()
{
	// TODO: check this assignment
	TimeLine syncTimeLine = *_simulation->settings.timeLine;

	std::list<double> startTimes;
	_simulation->bodyGroupList.DistinctStartTimes(startTimes, _simulation->settings.timeLine->Forward());

	// The first start time will be that of the massive bodies, since
	// the time field of syncTimeLine equals to that if enableDistinctStartTimes
	// is not enabled.
	if (!_simulation->settings.enableDistinctStartTimes) {
		startTimes.push_front(syncTimeLine.time);
	}

	std::list<BodyGroup *> listOfBG;
	for (std::list<double>::iterator it = startTimes.begin(); it != startTimes.end(); it++) {
		_simulation->bodyGroupList.FindBy((*it), listOfBG);
		for (std::list<BodyGroup *>::iterator bgIt = listOfBG.begin(); bgIt != listOfBG.end(); bgIt++) {
			if ((*bgIt)->inserted) continue;
			if (Insert(syncTimeLine.time, bgIt) == 1) {
				Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
				return 1;
			}
			(*bgIt)->inserted = true;
			_simulation->binary->Log("The bodygroup (" + (*bgIt)->description + ") with epoch: " + (*bgIt)->epoch + " was included", false);
		}
		if (++it != startTimes.end()) {
			syncTimeLine.length = *it - syncTimeLine.time;
			syncTimeLine.output = syncTimeLine.Forward() ? _simulation->settings.timeLine->output :
														  -_simulation->settings.timeLine->output;
			--it;
		}
		else
			break;

		if (syncTimeLine.length != 0.0 && _simulation->bodyList.size() > 1) {

			syncTimeLine.hDid = 0.0;
			syncTimeLine.hNext = syncTimeLine.Forward() ? ShortestPeriod() / 50.0 : -ShortestPeriod() / 50.0;
			if (fabs(syncTimeLine.hNext) > fabs(syncTimeLine.length)) {
				syncTimeLine.hNext = syncTimeLine.length;
			}

			// TODO: check whether the SOLARIS works with this comented lines of code
			// This was commented out since the Integrate() method calls the BodyListToBodyData() at the beginning.
			//if (BodyListToBodyData() == 1) {
			//	Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			//	return 1;
			//}

			if (Integrate(&syncTimeLine) == 1) {
				Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
				return 1;
			}

		}
		// TODO: the line below is only necessary for Test purposes. In real life the time
		// field of the syncTimeLine is updated by the Integrate() function.
		//syncTimeLine.time += syncTimeLine.length;

		listOfBG.clear();
	}

	_simulation->settings.timeLine->time = syncTimeLine.time;
	_simulation->settings.timeLine->save = syncTimeLine.time;

	return 0;
}

int Simulator::PreIntegration()
{
	// TODO: check this assignment
	TimeLine preTimeLine = *_simulation->settings.timeLine;

	preTimeLine.length = preTimeLine.start - preTimeLine.time;
	preTimeLine.output = preTimeLine.Forward() ? _simulation->settings.timeLine->output : -_simulation->settings.timeLine->output;

	preTimeLine.hDid = 0.0;
	preTimeLine.hNext  = preTimeLine.Forward() ? ShortestPeriod() / 50.0 : -ShortestPeriod() / 50.0;
	if (fabs(preTimeLine.hNext) > fabs(preTimeLine.length)) {
		preTimeLine.hNext = preTimeLine.length;
	}

	// If no synchronization was performed, than the bodyList is empty,
	// so it must be populated before the call to Integrate()
	if (_simulation->bodyList.size() == 0) {
		if (PopulateBodyList(preTimeLine.time) == 1) {
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			return 1;
		}
	}

	if (Integrate(&preTimeLine) == 1) {
		Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
		return 1;
	}
	
	return 0;
}

int Simulator::MainIntegration()
{
	// If neither synchronization nor pre-integration was performed, than the bodyList is empty,
	// so it must be populated before the call to Integrate()
	if (_simulation->bodyList.size() == 0) {
		if (PopulateBodyList(_simulation->settings.timeLine->time) == 1) {
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			return 1;
		}
	}

	_simulation->settings.timeLine->hDid = 0.0;
	_simulation->settings.timeLine->hNext = _simulation->settings.timeLine->Forward() ? ShortestPeriod() / 50000.0 : -ShortestPeriod() / 50000.0;
	if (fabs(_simulation->settings.timeLine->hNext) > fabs(_simulation->settings.timeLine->length)) {
		_simulation->settings.timeLine->hNext = _simulation->settings.timeLine->length;
	}

	if (Integrate(_simulation->settings.timeLine) == 1) {
		Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
		return 1;
	}

	return 0;
}

// Insert the bodies in the bgIt into the bodyList object sorted by BodyType
int Simulator::Insert(double time, std::list<BodyGroup>::iterator &bgIt)
{
	std::list<Body *>::iterator bodyListIt = _simulation->bodyList.begin();
	for (int i=0; i<NOfBodyType; i++ ) {
		if ((BodyType)i == UndefinedBodyType) continue;

		std::list<Body *> bodyListByType;
		bgIt->FindBy((BodyType)i, bodyListByType);
		if (bodyListByType.size() > 0) {
			_simulation->binary->SaveBodyProperties(time, bodyListByType);
			SetIteratorAfter((BodyType)i, bodyListIt);
			_simulation->bodyList.insert(bodyListIt, bodyListByType.begin(), bodyListByType.end());
		}
	}

	return 0;
}

int Simulator::Insert(double time, std::list<BodyGroup *>::iterator &bgIt)
{
	std::list<Body *>::iterator bodyListIt = _simulation->bodyList.begin();
	for (int i=0; i<NOfBodyType; i++ ) {
		if ((BodyType)i == UndefinedBodyType) continue;

		std::list<Body *> bodyListByType;
		(*bgIt)->FindBy((BodyType)i, bodyListByType);
		if (bodyListByType.size() > 0) {
			_simulation->binary->SaveBodyProperties(time, bodyListByType);
			SetIteratorAfter((BodyType)i, bodyListIt);
			_simulation->bodyList.insert(bodyListIt, bodyListByType.begin(), bodyListByType.end());
		}
	}

	return 0;
}

int Simulator::SetIteratorAfter(BodyType type, std::list<Body *>::iterator& it)
{
	while (it != _simulation->bodyList.end() && (*it)->type != type) {
		it++;
	}
	while (it != _simulation->bodyList.end() && (*it)->type == type) {
		it++;
	}

	return 0;
}

int Simulator::PopulateBodyList(double time)
{
	std::list<BodyGroup>::iterator bgIt = _simulation->bodyGroupList.items.begin();
	for ( ; bgIt != _simulation->bodyGroupList.items.end(); bgIt++) {
		if (Insert(time, bgIt) == 1) {
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			return 1;
		}
	}

	return 0;
}

double Simulator::ShortestPeriod()
{
	double centralGm = _simulation->bodyList.front()->GetGm();
	double period = 1.0e10;
	for (std::list<Body *>::iterator it = _simulation->bodyList.begin(); it != _simulation->bodyList.end(); it++) {
		if ((*it)->type == CentralBody)
			continue;
		double mu = centralGm + ((*it)->type == TestParticle ? 0.0 : (*it)->GetGm());
		double p = (*it)->CalculateOrbitalPeriod(mu);
		if (p > 0 && p < period)
			period = p;
	}
	return period;
}

int Simulator::BodyListToBodyData()
{
	bodyData.Free();

	bodyData.nBodies.Count(_simulation->bodyList);
	if (bodyData.Allocate() == 1) {
		Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
		return 1;
	}

	int i = 0;
	for (std::list<Body *>::iterator it = _simulation->bodyList.begin(); it != _simulation->bodyList.end(); it++) {
		bodyData.id[i]          = (*it)->GetId();
		bodyData.type[i]        = (*it)->type;
		bodyData.migStopAt[i]   = (*it)->migrationStopAt;
		bodyData.migType[i]     = (*it)->migrationType;

		if ((*it)->type != TestParticle) {
			bodyData.mass[i]    = (*it)->characteristics->mass;
			bodyData.radius[i]  = (*it)->characteristics->radius;
			bodyData.density[i] = (*it)->characteristics->density;
            bodyData.cD[i]      = (*it)->characteristics->stokes;

			if (bodyData.radius[i] > 0 ) { // && (*it)->characteristics->stokes > 0) {
				bodyData.gammaEpstein[i] = (*it)->characteristics->GammaEpstein();
                if ((*it)->characteristics->stokes > 0) {
                    bodyData.gammaStokes[i] = (*it)->characteristics->GammaStokes();
                } else {
    				bodyData.gammaStokes[i] = 0.0;
                }
			}
			else {
				bodyData.gammaEpstein[i]    = 0.0;
    			bodyData.gammaStokes[i]     = 0.0;
			}
		}
		else {
			bodyData.mass[i]        = 0.0;
			bodyData.radius[i]      = 0.0;
			bodyData.density[i]     = 0.0;
			bodyData.cD[i]          = 0.0;
			bodyData.gammaStokes[i] = 0.0;
			bodyData.gammaEpstein[i]= 0.0;
		}

		int i0 = 6*i;
		bodyData.y0[i0 + 0] = (*it)->phase->position.x;
		bodyData.y0[i0 + 1] = (*it)->phase->position.y;
		bodyData.y0[i0 + 2] = (*it)->phase->position.z;
		bodyData.y0[i0 + 3] = (*it)->phase->velocity.x;
		bodyData.y0[i0 + 4] = (*it)->phase->velocity.y;
		bodyData.y0[i0 + 5] = (*it)->phase->velocity.z;
		i++;
	}

	if (_simulation->settings.baryCentric) {
		Calculate::PhaseOfBC(&bodyData, bodyData.bc);
		Tools::ToPhase(bodyData.bc, &(bodyData.phaseOfBC));

		Calculate::PhaseWithRespectToBC(&bodyData, bodyData.bc);

		Calculate::AngularMomentum(&bodyData, &(bodyData.angularMomentum[0]));
		Calculate::KineticEnergy(&bodyData, bodyData.kineticEnergy[0]);
		Calculate::PotentialEnergy(&bodyData, bodyData.potentialEnergy[0]);
		bodyData.totalEnergy[0] = bodyData.kineticEnergy[0] - bodyData.potentialEnergy[0];
	}

	return 0;
}

void Simulator::UpdateBodyListAfterIntegration()
{
#ifdef _DEBUG
//	fprintf(stderr, "File: %40s, Function: %40s, Line: %10d\n", __FILE__, __FUNCTION__, __LINE__);
#endif
	int i = 0;
	for (std::list<Body *>::iterator it = _simulation->bodyList.begin(); it != _simulation->bodyList.end(); it++ ) {
		if ((*it)->type != TestParticle) {
			(*it)->characteristics->mass = bodyData.mass[i];
			(*it)->characteristics->radius = bodyData.radius[i];
			(*it)->characteristics->density = bodyData.density[i];
            (*it)->characteristics->stokes = bodyData.cD[i];
		}

		int i0 = 6*i;
		(*it)->phase->position.x = bodyData.y0[i0 + 0];
		(*it)->phase->position.y = bodyData.y0[i0 + 1];
		(*it)->phase->position.z = bodyData.y0[i0 + 2];
		(*it)->phase->velocity.x = bodyData.y0[i0 + 3];
		(*it)->phase->velocity.y = bodyData.y0[i0 + 4];
		(*it)->phase->velocity.z = bodyData.y0[i0 + 5];

		i++;
	}
}

int Simulator::CheckEvent(double timeOfEvent)
{
#ifdef _DEBUG
//	fprintf(stderr, "File: %40s, Function: %40s, Line: %10d\n", __FILE__, __FUNCTION__, __LINE__);
#endif
	static double ejection = _simulation->settings.ejection;
	static double hitCentrum = _simulation->settings.hitCentrum;
	static double e3 = ejection > 0 ? 1.0/(ejection*ejection*ejection) : 0.0;
	static double h3 = hitCentrum > 0 ? 1.0/(hitCentrum*hitCentrum*hitCentrum) : 0.0;

	for (int i=1; i<bodyData.nBodies.total; i++) {
		// If Ejection was set check if distance of the body is larger than it
		if (ejection > 0 && _acceleration->rm3[i] < e3) {
			int i0 = 6*i;
			TwoBodyAffair affair(Ejection, timeOfEvent, 0, i, bodyData.id[0], bodyData.id[i], bodyData.y0, &(bodyData.y0[i0]));
			_ejectionEvent.items.push_back(affair);
			_ejectionEvent.N++;
		}
		// If HitCentrum was set check if distance of the body is smaller than it
		if (hitCentrum > 0 && _acceleration->rm3[i] > h3) {
			int i0 = 6*i;
			TwoBodyAffair affair(HitCentrum, timeOfEvent, 0, i, bodyData.id[0], bodyData.id[i], bodyData.y0, &(bodyData.y0[i0]));
			_hitCentrumEvent.items.push_back(affair);
			_hitCentrumEvent.N++;
		}
	}

	if (_ejectionEvent.items.size() > 0) {
		_simulation->binary->SaveTwoBodyAffairs(_ejectionEvent.items);
		// Remove the ejected bodies from the simulation
		for (std::list<TwoBodyAffair>::iterator it = _ejectionEvent.items.begin(); it != _ejectionEvent.items.end(); it++) {
			std::ostringstream stream;
			stream << *it;
			_simulation->binary->Log(stream.str(), true);
			RemoveBody(it->body2Id);
		}
		_ejectionEvent.items.clear();
	}

	if (_hitCentrumEvent.items.size() > 0) {
		_simulation->binary->SaveTwoBodyAffairs(_hitCentrumEvent.items);
		// Remove the hit centrum bodies from the simulation
		for (std::list<TwoBodyAffair>::iterator it = _hitCentrumEvent.items.begin(); it != _hitCentrumEvent.items.end(); it++) {

			int survivIdx = -1;
			int mergerIdx = -1;
			int survivId = -1;
			int mergerId = -1;

			if (HandleCollision(it->idx1, it->idx2, survivIdx, mergerIdx, survivId, mergerId) == 1)	{
				Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
				return 1;
			}

			Body *body = _simulation->FindBy(survivId);
			body->characteristics->mass    = bodyData.mass[   survivIdx];
			body->characteristics->radius  = bodyData.radius[ survivIdx];
			body->characteristics->density = bodyData.density[survivIdx];
			body->characteristics->stokes  = bodyData.cD[     survivIdx];
			_simulation->binary->SaveVariableProperty(body, it->time);

			std::ostringstream stream;
			stream << *it;
			_simulation->binary->Log(stream.str(), true);
			RemoveBody(it->body2Id);
		}
		_hitCentrumEvent.items.clear();
	}

	if (_simulation->settings.collision != 0) {
		std::list<TwoBodyAffair> collisions;
		double factor = _simulation->settings.collision->factor;
		for (register int i=0; i < bodyData.nBodies.total; i++) {
			int j = bodyData.indexOfNN[i];
			if (j >= 0 && factor*(bodyData.radius[i] + bodyData.radius[j]) > bodyData.distanceOfNN[i]) {

				int survivIdx = -1;
				int mergerIdx = -1;
				int survivId = -1;
				int mergerId = -1;

				if (HandleCollision(i, j, survivIdx, mergerIdx, survivId, mergerId) == 1) {
					Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
					return 1;
				}

				double time = _simulation->settings.timeLine->time;
				// TODO: Check whether the affair class meeds the survivIdx, mergerIdx!
				TwoBodyAffair affair(Collision, time, survivIdx, mergerIdx, survivId, mergerId, &bodyData.y[6*survivIdx], &bodyData.y[6*mergerIdx]);
				collisions.push_back(affair);

				Body *body = _simulation->FindBy(survivId);
				body->characteristics->mass    = bodyData.mass[   survivIdx];
				body->characteristics->radius  = bodyData.radius[ survivIdx];
				body->characteristics->density = bodyData.density[survivIdx];
                body->characteristics->stokes  = bodyData.cD[     survivIdx];
				_simulation->binary->SaveVariableProperty(body, time);

				std::ostringstream stream;
				stream << "Collision: At " << time*Constants::DayToYear << " [yr] between body with id: " << survivId << " and id: " << mergerId;
				_simulation->binary->Log(stream.str(), true);
				// Remove the merger body
				RemoveBody(mergerId);

				// Since the nearest neighbor is symmetrical, i.e. if i is the neighbor of j then vice versa.
				// Therefore if the collision between i and j was handled, than the corresponding data in the
				// indexOfNN array must be set to -1 in order to avoid the artificial collision between j and i.
				bodyData.indexOfNN[mergerIdx] = -1;
			}
		}
		_simulation->binary->SaveTwoBodyAffairs(collisions);
	}

	return 0;
}

int Simulator::RemoveBody(int bodyId)
{
#ifdef _DEBUG
//	fprintf(stderr, "File: %40s, Function: %40s, Line: %10d\n", __FILE__, __FUNCTION__, __LINE__);
#endif

	int index = 0;
	for ( ; index < bodyData.nBodies.total; index++) {
		if (bodyData.id[index] == bodyId) {
			break;
		}
	}

	if (bodyData.nBodies.UpdateAfterRemove((BodyType)(bodyData.type[index])) == 1) {
		Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
		return 1;
	}

	for (int i=index; i<bodyData.nBodies.total; i++) {
		bodyData.id[i] = bodyData.id[i + 1];
		bodyData.type[i] = bodyData.type[i + 1];
		bodyData.migType[i] = bodyData.migType[i + 1];

		bodyData.mass[i] = bodyData.mass[i + 1];
		bodyData.radius[i] = bodyData.radius[i + 1];
		bodyData.density[i] = bodyData.density[i + 1];
		bodyData.gammaStokes[i] = bodyData.gammaStokes[i + 1];
		bodyData.gammaEpstein[i] = bodyData.gammaEpstein[i + 1];

		int i0 = 6*i;
		memcpy(&bodyData.y0[i0], &bodyData.y0[i0 + 6], 6*sizeof(double));
	}

	return 0;
}

int Simulator::HandleCollision(int idx1, int idx2, int& survivIdx, int &mergerIdx, int& survivId, int& mergerId)
{
	if (this->bodyData.mass[idx1] == this->bodyData.mass[idx2]) {
		survivIdx = idx1 < idx2 ? idx1 : idx2;
		mergerIdx = idx1 > idx2 ? idx1 : idx2;
	} else {
		survivIdx = this->bodyData.mass[idx1] > this->bodyData.mass[idx2] ? idx1 : idx2;
		mergerIdx = this->bodyData.mass[idx1] < this->bodyData.mass[idx2] ? idx1 : idx2;
	}

	if (bodyData.type[mergerIdx] == TestParticle)
		return 0;

	survivId = bodyData.id[survivIdx];
	mergerId = bodyData.id[mergerIdx];

	if (CalculatePhaseAfterCollision(survivIdx, mergerIdx) == 1) {
		Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
		return 1;
	}
	if (CalculateCharacteristicsAfterCollision(survivId, mergerId, survivIdx, mergerIdx) == 1) {
		Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
		return 1;
	}

	return 0;
}

int Simulator::CalculatePhaseAfterCollision(int survivIdx, int mergerIdx)
{
#ifdef _DEBUG
//	fprintf(stderr, "File: %40s, Function: %40s, Line: %10d\n", __FILE__, __FUNCTION__, __LINE__);
#endif
    // TODO: revisit this after the change of the equations of motions, i.e. when the ref. frame is inertial
	if (bodyData.type[survivIdx] != CentralBody) {
		double survivMass = bodyData.mass[survivIdx];
		double mergerMass = 0.0;
		if (bodyData.type[mergerIdx] == SuperPlanetesimal) {
			mergerMass = Characteristics::CalculateMass(bodyData.density[mergerIdx], bodyData.radius[mergerIdx]);
		}
		else {
			mergerMass = bodyData.mass[mergerIdx];
		}
		double M = survivMass + mergerMass;
		int i0 = 6*survivIdx;
		int j0 = 6*mergerIdx;
		for (int i=0; i<6; i++) {
			bodyData.y0[i0 + i] = 1.0/M*(survivMass * bodyData.y0[i0 + i] + mergerMass * bodyData.y0[j0 + i]);
		}
	}

	return 0;
}

int Simulator::CalculateCharacteristicsAfterCollision(int survivId, int mergerId, int survivIdx, int mergerIdx)
{
#ifdef _DEBUG
//	fprintf(stderr, "File: %40s, Function: %40s, Line: %10d\n", __FILE__, __FUNCTION__, __LINE__);
#endif
	double mergerMass = 0.0;
	double survivMass = bodyData.mass[survivIdx];
	if (bodyData.type[mergerIdx] == SuperPlanetesimal) {
		mergerMass = Characteristics::CalculateMass(bodyData.density[mergerIdx], bodyData.radius[mergerIdx]);
	}
	else {
		mergerMass = bodyData.mass[mergerIdx];
	}
	double totalMass = survivMass + mergerMass;

	// The original mass of the survivor body is needed to compute the composition list of the new body
	std::list<Component> *survivList = 0;
	std::list<Component> *mergerList = 0;
	std::list<Component> componentList;
	if (CalculateComponentListUnion(survivId, mergerId, &survivList, &mergerList, &componentList) == 1) {
		Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
		return 1;
	}
	if (componentList.size() > 0) {
		for (std::list<Component>::iterator it = componentList.begin(); it != componentList.end(); it++) {
			Component *cs = Tools::FindByName(survivList, it->name);
			Component *cm = Tools::FindByName(mergerList, it->name);

			if (cs != 0 && cm != 0)
			{
				it->ratio = 1.0 / totalMass * (survivMass*cs->ratio + mergerMass*cm->ratio);
			}
			else if (cs == 0 && cm != 0)
			{
				it->ratio = 1.0 / totalMass * mergerMass*cm->ratio;
			}
			else
			{
				it->ratio = 1.0 / totalMass * survivMass*cs->ratio;
			}
		}
		// Update the component list of the survivor body:
		Body *surviv = _simulation->FindBy(survivId);
		surviv->characteristics->componentList.clear();
		surviv->characteristics->componentList = componentList;
	}
	// Now it is possible to update the mass of the survivor body
	bodyData.mass[survivIdx] = totalMass;

    // Since the mass is obligatory (if it is not defined the XmlFileAdapter returns 1) than either the
	// radius or the density is defined the other is computed by the Integrator::ToBodyData() function.
	// Therefore it is enough to check only one of the characteristics.
	// NOTE: even if the radius and density of the merger is known, it will be not used.
	if (bodyData.radius[survivIdx] == 0) {
		return 0;
	}

	// If neither the radius nor the density of the merger is defined than use the density of the survivor,
	// and compute its radius
	if (bodyData.radius[mergerIdx] == 0) {
		bodyData.density[mergerIdx] = bodyData.density[survivIdx];
		bodyData.radius[mergerIdx] = pow(3.0/(4.0*Constants::Pi)*bodyData.mass[mergerIdx]/bodyData.density[mergerIdx], 1.0/3.0);
	}

	// Compute the volume of the two colliding bodies:
	double totalVolume = 4.0/3.0*Constants::Pi*(CUBE(bodyData.radius[survivIdx]) + CUBE(bodyData.radius[mergerIdx]));
	// Compute the radius of the new body assuming a spherical form:
	bodyData.radius[survivIdx] = pow(3.0/(4.0*Constants::Pi)*totalVolume, 1.0/3.0);
	// Compute the density of the new body:
	bodyData.density[survivIdx] = totalMass / totalVolume;

	return 0;
}

int Simulator::CalculateComponentListUnion(int survivId, int mergerId, std::list<Component> **survivList, std::list<Component> **mergerList, std::list<Component> *result)
{
#ifdef _DEBUG
//	fprintf(stderr, "File: %40s, Function: %40s, Line: %10d\n", __FILE__, __FUNCTION__, __LINE__);
#endif
	Body *surviv = _simulation->FindBy(survivId);
	if (surviv == 0) {
		Error::_errMsg = "The body was not found by FindBy(survivId)!";
		Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
		return 1;
	}
	Body *merger = _simulation->FindBy(mergerId);
	if (merger == 0) {
		Error::_errMsg = "The body was not found by FindBy(mergerId)!";
		Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
		return 1;
	}

	*survivList = &(surviv->characteristics->componentList);
	*mergerList = &(merger->characteristics->componentList);
	if ((*survivList)->size() == 0 && (*mergerList)->size() == 0)
		return 0;

	if ((*survivList)->size() == 0) {
		Component c("Unknown", 100);
		surviv->characteristics->componentList.push_back(c);
	}
	if ((*mergerList)->size() == 0) {
		Component c("Unknown", 100);
		merger->characteristics->componentList.push_back(c);
	}

	Tools::MergeComponentList(*survivList, *mergerList, result);

	return 0;
}

