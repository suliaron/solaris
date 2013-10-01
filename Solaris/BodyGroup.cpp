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

int BodyGroup::CountBy(BodyType type)
{
	int	result = 0;
	for (std::list<Body>::iterator bodyIterator = items.begin(); bodyIterator != items.end(); bodyIterator++) {
		if (bodyIterator->type == type)
			result++;
	}
	return result;
}

int BodyGroup::CountBy(double mass)
{
    int result = 0;
	for (std::list<Body>::iterator bodyIterator = items.begin(); bodyIterator != items.end(); bodyIterator++) {
		if (bodyIterator->type == TestParticle)
			continue;
		if (bodyIterator->characteristics->mass > mass)
			result++;
	}

    return result;
}

void BodyGroup::FindBy(BodyType type, std::list<Body *> &result)
{
	for (std::list<Body>::iterator it = items.begin(); it != items.end(); it++) {
		if (it->type == type) {
			result.push_back(&(*it));
		}
	}
}

int BodyGroup::CalculateStartTime(double tau)
{
    // The Epoch property is defined
    if (!epoch.empty())
    {
		if (Ephemeris::GetJulianDate(epoch, this->startTime) == 1) {
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			return 1;
		}
        // The offset property is defined
        if (offset != 0.0) {
            this->startTime += offset;
        }
    }
    // The Epoch property is undefined
    else
    {
        this->startTime = tau;
        // The offset property is defined
        if (offset != 0.0) {
            this->startTime += offset;
        }
    }
	
	return 0;
}

/// <summary>
/// Returns true if the BodyGroup contains at least one body with mass. The
/// central body is not considered.
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

int BodyGroup::ToBodyList(std::list<Body *> &result)
{
	for (std::list<Body>::iterator it = items.begin(); it != items.end(); it++) {
		result.push_back(&(*it));
	}

	return 0;
}
