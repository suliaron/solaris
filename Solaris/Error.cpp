#include <iostream>
#include <string>
#include <stdio.h>
#include <list>
#include <string>

#include "Constants.h"
#include "Error.h"

std::list<std::string> Error::_stackTrace;
std::string Error::_errMsg;

void Error::PushLocation(const char *file, const char *function, long line)
{
	char		location[512];

	sprintf(location, "File: %s, function: %s, line: %ld.", file, function, line);

	_stackTrace.push_back(std::string(location));
}

void Error::PrintStackTrace()
{
	std::cerr << Constants::CodeName << ": " << _errMsg << std::endl;
	// Note: this is not functioning under windows!!
	for (std::list<std::string>::iterator it = _stackTrace.end(); it != _stackTrace.begin(); it--) {
		std::cerr << *it << std::endl;
	}
}
