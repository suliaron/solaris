#include "TimeLine.h"

TimeLine::TimeLine()
{
	startTimeDefined = false;

	start = 0.0;
	length = 0.0;
	output = 0.0;

	time = 0.0;
	save = 0.0;

	hNext = 0.0;
	hDid = 0.0;
}
