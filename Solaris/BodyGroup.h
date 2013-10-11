#ifndef BODYGROUP_H_
#define BODYGROUP_H_

#include <string>
#include <list>

#include "Body.h"

class BodyGroup
{
public:
	BodyGroup();
	BodyGroup(std::string guid);
	BodyGroup(std::string guid, std::string description);
	BodyGroup(std::string guid, std::string description, std::string epoch);
	BodyGroup(std::string guid, std::string description, std::string epoch, double offest);
	BodyGroup(std::string guid, std::string description, std::string epoch, double offest, std::string referenceFrame);

	int		CountBy(BodyType type);
	int		CountBy(double mass);
	
	void	FindBy(BodyType type, std::list<Body *> &result);
	void	BodyListBy(BodyType type, std::list<Body *> &result);

	int		ToBodyList(std::list<Body *> &result);

	int		SetStartTime(double startTimeOfMainSimulation);
	bool	ContainsMassiveBody();

	// This is needed for the std::list<BodyGroup> remove() method
	bool operator==(const BodyGroup &rhs)
    {
		return (rhs._id == _id);
    }

	std::list<Body> items;

	double offset;
	double startTime;

	std::string description;
	std::string epoch;
	std::string referenceFrame;
	std::string guid;

	// This flag is needed by the pre-integration process
	bool inserted;

private:
	// id of the bodyGroup instance
	int			_id;

	// id of the class BodyGroup, its value will be assigned to the next instance of BodyGroup
	// guaranteeing the uniqueness of the _id field of each BodyGroup object.
	static int _bodyGroupId;
};

#endif
