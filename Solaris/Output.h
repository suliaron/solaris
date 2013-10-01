#ifndef OUTPUT_H_
#define OUTPUT_H_

#include <string>

class Output {
public:
	Output();

	std::string GetPath(const std::string fileName);

	static std::string directory;
	static char directorySeparator;

	std::string phases;
	std::string integrals;
	std::string constantProperties;
	std::string variableProperties;
	std::string compositionProperties;
	std::string twoBodyAffair;
	std::string log;
};

#endif