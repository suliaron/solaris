#ifndef BODYGROUPLIST_H_
#define BODYGROUPLIST_H_

#include <list>
#include "BodyGroup.h"

enum FL {
    First,
    Last
};

class BodyGroupList
{
public:
	BodyGroupList();

	int		GetEpoch(double &firstEpoch, FL fl);
	int		DistinctEpochs(std::list<double> &epochs);
	int		DistinctEpochs(std::list<double> &epochs, bool increasing);

	double	GetStartTime(FL fl);
	int		SetStartTime(double startTimeOfMainIntegrationPhase);
	void	DistinctStartTimes(std::list<double> &startTimes, bool increasing);
	int		DistinctStartTimesOfMassiveBodies(int &count);

	int		CountBy(BodyType type);

	void	FindBy(BodyType type, std::list<Body *> &result);
	void	FindBy(double startTime, std::list<BodyGroup *> &result);

	bool	GetBodyGroupWithMassiveBodies(std::list<BodyGroup>::iterator &it);

	int		DistinctReferenceFrame(std::list<std::string> &referenceFrames);

	// This is needed for the sort method of the list object
	static bool CompareStartTime(BodyGroup bg1, BodyGroup bg2);

	std::list<BodyGroup> items;
	int nOfDistinctStartTimes;
};

#endif
