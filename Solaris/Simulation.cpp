#include <stdio.h> // TODO: remove this after debug
#include <sstream>

#include "BinaryFileAdapter.h"
#include "Ephemeris.h"
#include "Error.h"
#include "Integrator.h"
#include "Simulation.h"
#include "Settings.h"
#include "TimeLine.h"

Simulation::Simulation()
{
	nebula		        = 0;
	binary		        = 0;
    fargoParameters     = 0;
}

Simulation::Simulation(std::string runType)
{
    this->runType       = runType;

	nebula		        = 0;
	binary		        = 0;
    fargoParameters     = 0;
}

int Simulation::Initialize()
{
	if (InitializeTimeLineAndBodyGroups() == 1) {
		Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
		return 1;
	}

	if (InitializePhases() == 1) {
		Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
		return 1;
	}

	return 0;
}

int Simulation::InitializeTimeLineAndBodyGroups()
{
    // Compute the start time of the main integration phase
    // If the start attribute in the TimeLine tag was not defined in the xml, than it will be computed from the epochs of the BodyGroups
	if (!settings.timeLine->startTimeDefined)
    {
		SetStartTimeOfMainPhase();
    }

	if (bodyGroupList.SetStartTime(settings.timeLine->start) == 1) {
		Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
		return 1;
	}

	if (CheckStartTimes() == 1) {
		Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
		return 1;
	}

	// Sort the Body Groups in the BodyGroupList based on their start time into increasing order.
	this->bodyGroupList.items.sort(BodyGroupList::CompareStartTime);
	if (!settings.timeLine->Forward()) {
		bodyGroupList.items.reverse();
	}

	if (SetTimeAndSaveOfTimeLine() == 1) {
		Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
		return 1;
	}

	/*
	* If the enableDistinctStartTimes is false, than the Time and Save
	* property of the TimeLine class must be set equal to the
	* startTime of the BodyGroup containing the massive bodies.
	*/
	//double time = 0.0;
	//if (!settings.enableDistinctStartTimes)
	//{
	//	std::list<BodyGroup>::iterator it;
	//	if (!bodyGroupList.GetBodyGroupWithMassiveBodies(it)) {
	//		binary->Log("None of the BodyGroups contains massive bodies!", true);
	//		Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
	//		return 1;
	//	}
	//	time = it->startTime;
	//}
	//// Because of the sort operation above the time equals to the startTime field of the first BodyGroup
	//else {
	//	time = bodyGroupList.items.front().startTime;
	//}
	//settings.timeLine->time = time;
	//settings.timeLine->save = time;

	return 0;
}

/// Sets the start time of the main integration phase based on the epochs of the body groups.
/// If the integration is a forward integration (i.e. length > 0) than the latest epoch will
/// be the start time, if it is a backward integration than the earliset epoch will be.
/// If none of the body groups contain epoch the start time is set to zero.
int Simulation::SetStartTimeOfMainPhase()
{
	double start = 0.0;
    // If the length attribute in the TimeLine tag is positive (forward integration) than
    // the start time is the last epoch, 
    int result = bodyGroupList.GetEpoch(start, settings.timeLine->Forward() ? Last : First);

	// No epochs were defined for the BodyGroups
	if (     result == -1) {
		settings.timeLine->start = 0.0;
	}
	// One or more epochs were found. If the direction of
	// the main integration is Forward in time, than the start time will be the latest epoch,
	// if the direction is Backward than it will be the first epoch.
	else if (result == 0) {
		settings.timeLine->start = start;
	}
	// An error occurred during the calculation.
	else {
		Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
		return 1;
	}

	return 0;
}

/// If enableDistinctStartTimes is false than the start time of the TimeLine will be
/// the start time of the Body Group that contains the massive bodies. Otherwise the
/// start time of the TimeLine will be the start time of the earliest (forward integration)
/// or latest (backward) start time.
int Simulation::SetTimeAndSaveOfTimeLine()
{
	double time = 0.0;
	if (settings.enableDistinctStartTimes == false)
	{
		std::list<BodyGroup>::iterator it;
		if (bodyGroupList.GetBodyGroupWithMassiveBodies(it) == false) {
			binary->Log("None of the BodyGroups contains massive bodies!", true);
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			return 1;
		}
		time = it->startTime;
	}
	// Because of the sort operation above the time equals to the startTime field of the first BodyGroup
	else {
		time = bodyGroupList.items.front().startTime;
	}
	settings.timeLine->time = time;
	settings.timeLine->save = time;

	return 0;
}

int Simulation::InitializePhases()
{
	Body *centralBody = FindCentralBody();
	for (std::list<BodyGroup>::iterator bgIt = this->bodyGroupList.items.begin(); bgIt != this->bodyGroupList.items.end(); bgIt++) {
		for (std::list<Body>::iterator bIt = bgIt->items.begin(); bIt != bgIt->items.end(); bIt++) {

			BodyType type = bIt->type;
			if (type != CentralBody && bIt->phase == 0) {
				if (centralBody == 0) {
					Error::_errMsg = "The central body is not defined, the phase cannot be computed!";
					Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
					return 1;
				}
				// Compute phase
				double mu = centralBody->GetGm() + (type != TestParticle ? bIt->GetGm() : 0.0);
				Phase phase(bIt->GetId());
				if (Ephemeris::CalculatePhase(mu, bIt->orbitalElement, &phase) == 1) {
					Error::_errMsg = "The phase could not be computed for body with Id: ";
					Error::_errMsg += static_cast<std::ostringstream *>(&(std::ostringstream() << bIt->GetId()))->str() + "!";
					Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
					return 1;
				}
				bIt->phase = new Phase(phase);
			}

			if (type == TestParticle || type == SuperPlanetesimal)
				continue;
			if (bIt->characteristics->radius == 0.0 && bIt->characteristics->density == 0.0)
				continue;
			if (bIt->characteristics->radius == 0.0) {
				bIt->characteristics->radius = bIt->characteristics->CalculateRadius();
			}
			else {
				bIt->characteristics->density = bIt->characteristics->CalculateDensity();
			}
		}
	}

	return 0;
}

// TODO: no calls to this function
/**
 * Sort the BodyGroups in the BodyGroupList into increasing order
 */
void Simulation::SortBodyGroupsByStartTime()
{
	this->bodyGroupList.items.sort(BodyGroupList::CompareStartTime);
}

int Simulation::CheckBodyGroupList()
{
	if (this->bodyGroupList.items.size() == 0) {
		Error::_errMsg = "There is no BodyGroup defined in the Simulation!";
		Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
		return 1;
	}
	if (this->bodyGroupList.items.size() == 1 && this->bodyGroupList.items.front().items.size() <= 1) {
		Error::_errMsg = "There is only one body defined in the Simulation!";
		Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
		return 1;
	}

	// Checks whether the BodyGroupList contains one and only one central body.
    // If zero or more than 1 are defined, an ApplicationException is thrown.
	int nCB = this->bodyGroupList.CountBy(CentralBody);
    if (nCB > 1)
    {
		Error::_errMsg = "There are more than 1 central body defined!";
		Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
		return 1;
    }
    if (nCB == 0)
    {
		Error::_errMsg = "The central body is not defined!";
		Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
		return 1;
    }

	std::list<std::string> distinctReferenceFrame;
	this->bodyGroupList.DistinctReferenceFrame(distinctReferenceFrame);
	if (distinctReferenceFrame.size() > 1)
    {
		Error::_errMsg = "The BodyGroups have different reference frames!";
		Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
		return 1;
    }
	if (distinctReferenceFrame.size() == 1)
    {
		this->referenceFrame = distinctReferenceFrame.front();
    }

	return 0;
}

/// If the enableDistinctStartTimes attribute is false and more than one BodyGroup contains
/// massive bodies with different start times, it returns with an error and logs it to the
/// log file.
int Simulation::CheckStartTimes()
{
	// If distinct start times for massive bodies was not enabled than massive bodies must have the same start time
	int count = 0;
	if (bodyGroupList.DistinctStartTimesOfMassiveBodies(count) == 1) {
		Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
		return 1;
	}
	if (!settings.enableDistinctStartTimes && count > 1)
    {
		binary->Log("More than 1 BodyGroup contain massive bodies with different start times!", true);
		Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
		return 1;
    }

	return 0;
}

// Calculate the start time of the simulation using the epochs of the BodyGroups.
//int Simulation::CalculateStartTime(double &start)
//{
//	return this->settings->timeLine->Forward() ? this->bodyGroupList.LastEpoch(start) : this->bodyGroupList.FirstEpoch(start);
//}

Body* Simulation::FindBy(int id)
{
	for (std::list<Body *>::iterator it = bodyList.begin(); it != bodyList.end(); it++) {
		if ((*it)->GetId() == id) {
			return *it;
		}
	}
	return 0;
}

Body* Simulation::FindCentralBody()
{
	std::list<Body *> centralBodyList;
	bodyGroupList.FindBy(CentralBody, centralBodyList);
	if (centralBodyList.size() == 0) {
		return 0;
	}
	else {
		return centralBodyList.front();
	}
}
