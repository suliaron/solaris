#include <cstdio>

#include "BodyGroup.h"
#include "Ephemeris.h"
#include "Error.h"

int BodyGroup::_bodyGroupId = 0;

BodyGroup::BodyGroup()
{
	this->_id = BodyGroup::_bodyGroupId++;

	this->offset = 0.0;
	this->startTime = 0.0;

	inserted = false;
}

BodyGroup::BodyGroup(std::string guid)
{
	this->_id = BodyGroup::_bodyGroupId++;

	this->guid = guid;

	this->offset = 0.0;
	this->startTime = 0.0;

	inserted = false;
}

BodyGroup::BodyGroup(std::string guid, std::string description)
{
	this->_id = BodyGroup::_bodyGroupId++;

	this->guid = guid;
    this->description = description;

	this->offset = 0.0;
	this->startTime = 0.0;

	inserted = false;
}

BodyGroup::BodyGroup(std::string guid, std::string description, std::string epoch)
{
	this->_id = BodyGroup::_bodyGroupId++;

	this->guid = guid;
    this->description = description;
    this->epoch = epoch;

	this->offset = 0.0;
	this->startTime = 0.0;

	inserted = false;
}

BodyGroup::BodyGroup(std::string guid, std::string description, std::string epoch, double offest)
{
	this->_id = BodyGroup::_bodyGroupId++;

	this->guid = guid;
    this->description = description;
    this->epoch = epoch;
    this->offset = offest;

	this->startTime = 0.0;

	inserted = false;
}

BodyGroup::BodyGroup(std::string guid, std::string description, std::string epoch, double offest, std::string referenceFrame)
{
	this->_id = BodyGroup::_bodyGroupId++;

    this->guid = guid;
    this->description = description;
    this->epoch = epoch;
    this->offset = offest;
    this->referenceFrame = referenceFrame;

	this->startTime = 0.0;

	inserted = false;
}

/// Counts how many bodies' type is equal to the specified type.
/// The number of these bodies are returned. 
int BodyGroup::CountBy(BodyType type)
{
	int	result = 0;
	for (std::list<Body>::iterator bodyIterator = items.begin(); bodyIterator != items.end(); bodyIterator++) {
		if (bodyIterator->type == type)
			result++;
	}
	return result;
}


/// Counts how many bodies' mass is equal to or greater than the value of the parameter 'mass'.
/// The number of these bodies are returned. 
int BodyGroup::CountBy(double mass)
{
    int result = 0;
	for (std::list<Body>::iterator bodyIterator = items.begin(); bodyIterator != items.end(); bodyIterator++) {
		if (bodyIterator->type == TestParticle)
			continue;
		if (bodyIterator->characteristics->mass >= mass)
			result++;
	}

    return result;
}

/// Iterates over the Bodies and compares its type with the specified type.
/// In case of equality the reference of the Body is added to the result.
void BodyGroup::FindBy(BodyType type, std::list<Body *> &result)
{
	for (std::list<Body>::iterator it = items.begin(); it != items.end(); it++) {
		if (it->type == type) {
			result.push_back(&(*it));
		}
	}
}

/// The start time is the time instant when the bodies of the BodyGroup object became part of
/// the simulation.
/// The start time is determined by the epoch and the offset fields of the BodyGroup object.
/// 4 possible cases exists:
/// 1. both epoch and offset are defined: start time = from the epoch the julian date is computed + offset
/// 2. only epoch is defined: start time = from the epoch the julian date is computed
/// 3. only offset is defined: start time = startTimeOfMainSimulation + offset
/// 4. none of them are defined: start time = startTimeOfMainSimulation
int BodyGroup::SetStartTime(double startTimeOfMainSimulation)
{
    // The Epoch property is defined
    if (!epoch.empty())
    {
		if (Ephemeris::GetJulianDate(epoch, startTime) == 1) {
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			return 1;
		}
    }
    else
    {
        startTime = startTimeOfMainSimulation;
    }
    startTime += offset;
	
	return 0;
}

/// <summary>
/// Returns true if the BodyGroup contains at least one body with mass.
/// The central body is ignored.
/// </summary>
bool BodyGroup::ContainsMassiveBody()
{
	for (std::list<Body>::iterator it = items.begin(); it != items.end(); it++) {
		if (it->type == TestParticle)
			continue;

		if (it->type != CentralBody && it->characteristics->mass > 0)
			return true;
	}
	return false;
}

/// Iterates over the Bodies and stores their reference in the result parameter.
int BodyGroup::ToBodyList(std::list<Body *> &result)
{
	for (std::list<Body>::iterator it = items.begin(); it != items.end(); it++) {
		result.push_back(&(*it));
	}

	return 0;
}
