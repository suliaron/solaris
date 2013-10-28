#include "Counter.h"

Counter::Counter()
{
	succededStep	= 0ull;
	failedStep		= 0ull;
	ejection		= 0ull;
	hitCentrum		= 0ull;
	collision		= 0ull;
}

std::ostream& operator<<(std::ostream& output, Counter counter)
{
	output << "succeded steps: " << counter.succededStep << std::endl;
	output << "  failed steps: " << counter.failedStep << std::endl;
	output << "      ejection: " << counter.ejection << std::endl;
	output << "   hit centrum: " << counter.hitCentrum << std::endl;
	output << "     collision: " << counter.collision << std::endl;

	return output;
}
