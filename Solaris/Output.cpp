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
