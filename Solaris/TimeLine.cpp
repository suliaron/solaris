#include "Constants.h"
#include "TimeLine.h"

TimeLine::TimeLine()
{
	startTimeDefined = false;

	start		= 0.0;
	length		= 0.0;
	output		= 0.0;

	save		= 0.0;
	time		= 0.0;
	millenium	= 0;

	hNext		= 0.0;
	hDid		= 0.0;

	elapsedTime = 0.0;
	lastSave	= 0.0;
	lastNSteps	= 0.0;
}

std::ostream& operator<<(std::ostream& output, TimeLine timeLine)
{
	output << "Elapsed  time: " << timeLine.elapsedTime * Constants::DayToYear << " [yr]" << std::endl;
	output << "  Actual time: " << timeLine.millenium * 1000.0 + timeLine.time * Constants::DayToYear << " [yr]" << std::endl;
	output << "Relative time: " << timeLine.time * Constants::DayToYear << " [yr]" << std::endl;
	output << "    millenium: " << timeLine.millenium << std::endl;
	output << "        hNext: " << timeLine.hNext << " [day]"  << std::endl;
	output << "         hDid: " << timeLine.hDid << " [day]"  << std::endl;

	return output;
}