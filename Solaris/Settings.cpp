#include <cstdio>

#include "Settings.h"
#include "Output.h"

Settings::Settings() {
	baryCentric					= false;
	enableDistinctStartTimes	= false;

	timeLine		= 0;
	collision		= 0;
	closeEncounter	= 0;
	weakCapture		= 0;

	ejection		= 0.0;
	hitCentrum		= 0.0;
}
