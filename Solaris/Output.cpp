#include "Output.h"

std::string Output::directory;
char Output::directorySeparator;

Output::Output() {
	phases					= "Phases.dat";
	integrals				= "Integrals.dat";
    constantProperties		= "ConstantProperties.dat";
    variableProperties		= "VariableProperties.dat";
    compositionProperties	= "CompositionProperties.dat";
    twoBodyAffair			= "TwoBodyAffair.dat";
    log						= "Log.txt";
}

std::string Output::GetPath(const std::string fileName)
{
	return Output::directory + Output::directorySeparator + fileName;
}

std::string Output::GetFilenameWithoutExt(const std::string& path)
{
	std::string result;

	if (path.size() > 0)
	{
		size_t pos = path.find_last_of("/\\");
		result = path.substr(pos + 1);
		pos = result.find_last_of('.');
		result = result.substr(0, pos);
	}

	return result;
}

std::string CombinePath(std::string dir, std::string filename)
{
	if (dir.size() > 0) {
		if (*(dir.end() - 1) != '/' && *(dir.end() - 1) != '\\') {
			return dir + '/' + filename;
		}
		else {
			return dir + filename;
		}
	}
	else {
		return filename;
	}
}
