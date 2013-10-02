#include <sstream>

#include "BodyGroupList.h"
#include "Ephemeris.h"
#include "Error.h"

BodyGroupList::BodyGroupList()
{
	nOfDistinctStartTimes = 0;
}

/// Depending on the value of the FL the epoch parameter will receive the first or last epoch defined in the Body Groups.
/// If none of the Body Groups contain epoch it returns -1.
int BodyGroupList::GetEpoch(double &epoch, FL fl)
{
	std::list<double> distinctEpochs;
    if (DistinctEpochs(distinctEpochs) == 1) {
		Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
		return 1;
    }
    if (distinctEpochs.size() > 0) {
        if (fl == Last) {
            distinctEpochs.reverse();
        }
        epoch = distinctEpochs.front();
        return 0;
    }
    else {
        return -1;
    }
}

// Compute the start time of each BodyGroup: at the StartTime they became part of the integration
int BodyGroupList::CalculateStartTime(double startTimeOfTheSimulation)
{	
	for (std::list<BodyGroup>::iterator it = this->items.begin(); it != this->items.end(); it++) {
		if (it->CalculateStartTime(startTimeOfTheSimulation) == 1) {
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			return 1;
		}
	}

	return 0;
}

double BodyGroupList::GetStartTime(FL fl)
{
	std::list<double> distinctStartTimes;
    DistinctStartTimes(distinctStartTimes, fl == First ? true : false);
	return distinctStartTimes.front();
}

/// <summary>
/// Returns the first StartTime defined in the BodyGroups.
/// </summary>
double BodyGroupList::FirstStartTime()
{
	std::list<double> distinctStartTimes;
	// The start times of the BodyGroups will be stored in the distinctStartTimes list in increasing order
	DistinctStartTimes(distinctStartTimes, true);
	return distinctStartTimes.front();
}

/// <summary>
/// Returns the last StartTime defined in the BodyGroups.
/// </summary>
double BodyGroupList::LastStartTime()
{
	std::list<double> distinctStartTimes;
	// The start times of the BodyGroups will be stored in the distinctStartTimes list in decreasing order
	DistinctStartTimes(distinctStartTimes, false);
	return distinctStartTimes.front();
}

int BodyGroupList::CountBy(BodyType type)
{
	int	result = 0;
	std::list<BodyGroup>::iterator bodyGroupIterator;
	
	for (bodyGroupIterator = this->items.begin(); bodyGroupIterator != this->items.end(); bodyGroupIterator++) {
		result += bodyGroupIterator->CountBy(type);
	}
	return result;
}

void BodyGroupList::FindBy(double startTime, std::list<BodyGroup *> &result)
{
	for (std::list<BodyGroup>::iterator it = items.begin(); it != items.end(); it++) {
		if (it->startTime == startTime) {
			result.push_back(&(*it));
		}
	}
}

void BodyGroupList::FindBy(BodyType type, std::list<Body *> &result)
{
	for (std::list<BodyGroup>::iterator it = items.begin(); it != items.end(); it++) {
		it->FindBy(type, result);
	}
}

/**
 * Set the iterator to point to the body group containing the massive bodies.
 *
 * @param it the iterator which iterates over the body groups. It will point to the first body group in the list
 *        which contains massive bodies.
 * @return true if the list contains at least one body groups with massive bodies.
 */
bool BodyGroupList::GetBodyGroupWithMassiveBodies(std::list<BodyGroup>::iterator &it)
{
	for (it = items.begin(); it != items.end(); it++) {
		if (it->ContainsMassiveBody()) {
			return true;
		}
	}
	return false;
}

int BodyGroupList::DistinctReferenceFrame(std::list<std::string> &referenceFrames)
{
	std::list<BodyGroup>::iterator bodyGroupIterator;
	
	for (bodyGroupIterator = this->items.begin(); bodyGroupIterator != this->items.end(); bodyGroupIterator++) {
		if (!bodyGroupIterator->referenceFrame.empty()) {
			referenceFrames.push_back(bodyGroupIterator->referenceFrame);
		}
	}
	referenceFrames.sort();
	referenceFrames.unique();

	return 0;
}

/// <summary>
/// Returns a sorted list (increasing) of the distinct julian dates (epochs) of the BodyGroups.
/// </summary>
/// <returns>List of epochs in increasing or decreasing order</returns>
int BodyGroupList::DistinctEpochs(std::list<double> &epochs)
{
	for (std::list<BodyGroup>::iterator it = items.begin(); it != items.end(); it++) {
		std::string epoch = it->epoch;
		if (!epoch.empty()) {
			double jd = 0.0;
			if (Ephemeris::GetJulianDate(epoch, jd) == 1) {
				Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
				return 1;
			}
			epochs.push_back(jd);
		}
	}
	epochs.sort();
	epochs.unique();

    return 0;
}

/// <summary>
/// Returns a List of the distinct julian dates (epochs) of the BodyGroups and sort 
/// the list into increasing or decreasing order depending on the value of the increasing
/// parameter.
/// </summary>
/// <param name="increasing">The order in which the epochs will be sorted</param>
/// <returns>List of epochs in increasing or decreasing order</returns>
int BodyGroupList::DistinctEpochs(std::list<double> &epochs, bool increasing)
{
    if (DistinctEpochs(epochs) == 1) {
		Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
		return 1;
    }

	if (!increasing)
		epochs.reverse();

	return 0;
}

void BodyGroupList::DistinctStartTimes(std::list<double> &startTimes, bool increasing)
{
	for (std::list<BodyGroup>::iterator it = items.begin(); it != items.end(); it++) {
		startTimes.push_back(it->startTime);
	}
    // Sort the start times in increasing order
    startTimes.sort();
	startTimes.unique();
    if (!increasing)
    {
        startTimes.reverse();
    }

	if (nOfDistinctStartTimes == 0) nOfDistinctStartTimes = startTimes.size();
}

int BodyGroupList::CountDistinctStartTimesForMassiveBodies(int &total)
{
	std::list<double> startTimes;
	DistinctStartTimes(startTimes, true);

	std::list<BodyGroup *> list;
	for (std::list<double>::iterator it = startTimes.begin(); it != startTimes.end(); it++) {
		FindBy(*it, list);
		if (list.size() == 0)
		{
			std::ostringstream stream;
			stream << "Could not find the body group with start time = " << *it << "!";
			Error::_errMsg = stream.str();
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			return 1;
		}
		int partial = 0;
		for (std::list<BodyGroup *>::iterator bgIt = list.begin(); bgIt != list.end(); bgIt++) {
            // The number of massive bodies in BodyGroup
			partial += (*bgIt)->CountBy(0.0);
            // The number of the central body must be subtracted
            partial -= (*bgIt)->CountBy(CentralBody);
		}
		list.clear();
		total += partial > 0 ? 1 : 0;
	}

	return 0;
}

bool BodyGroupList::CompareStartTime(BodyGroup bg1, BodyGroup bg2)
{
	return bg1.startTime < bg2.startTime ? true : false;
}
