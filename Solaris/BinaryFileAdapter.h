#ifndef BINARYFILEADAPTER_H_
#define BINARYFILEADAPTER_H_

#include <list>

class Body;
class Output;
class TwoBodyAffair;

class BinaryFileAdapter
{
public:
	enum OutputType {
		BINARY,
		TEXT
	};

	BinaryFileAdapter(Output *output);

    bool    FileExists(const std::string& name);

	void	Log(std::string msg, bool printToScreen);
	void	LogStartParameters(int argc, char* argv[]);
	void	LogTimeSpan(std::string msg, time_t startTime);

	void	SavePhases(double time, int n, double *y, int *id, OutputType type);
	void	SavePhase(std::ofstream& writer, double *y, int *id, OutputType type);

	void	SaveIntegrals(double time, int n, double *integrals);

	void	SaveTwoBodyAffairs(std::list<TwoBodyAffair>& list);
	void	SaveTwoBodyAffair(std::ofstream& writer, TwoBodyAffair& affair);

	void	SaveBodyProperties(double time, std::list<Body *>& bodyList);

	void	SaveConstantProperty(std::ofstream& writer, Body* body);
	void	SaveConstantProperty(Body* body);

	void	SaveVariableProperty(Body* body, double time);
	void	SaveVariableProperty(std::ofstream& writer, Body* body, double time);

	void	SaveCompositionProperty(std::ofstream& writer, Body* body, int propertyId);

private:
	std::string _errMsg;
	Output		*output;
	static int	_propertyId;
	static int	_compositionId;
};

#endif
