#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cstdio>
#include <list>
#include <string>
#include <string.h>
#include <time.h>
#include <sys/stat.h>

#include "BinaryFileAdapter.h"
#include "Body.h"
#include "BodyGroup.h"
#include "Output.h"
#include "Phase.h"
#include "Tools.h"
#include "TwoBodyAffair.h"

int BinaryFileAdapter::_propertyId = 0;
int BinaryFileAdapter::_compositionId = 0;

BinaryFileAdapter::BinaryFileAdapter(Output *output) 
{
	this->output = output;
}

bool BinaryFileAdapter::FileExists(const std::string& name) {
    struct stat buffer;
    std::string path = output->GetPath(name);
    return (stat (path.c_str(), &buffer) == 0); 
}

/**
 * Write the message to the log file and if printToScreen = true than also to the screen. On error
 * the function exits to the system.
 *
 * @param msg message to write into the log file
 * @param printToScreen if true the message will be also printed to the screen
 */
void BinaryFileAdapter::Log(std::string msg, bool printToScreen)
{
	char dateTime[20];
	time_t now = time(0);
	strftime(dateTime, 20, "%Y-%m-%d %H:%M:%S", localtime(&now));

	std::string path = output->GetPath(output->log);
	std::ofstream logger(path.c_str(), std::ios::out | std::ios::app);
	if (logger) {
		logger << dateTime << " " << msg << std::endl;
		if (logger.bad()) {
			perror(("error while reading file " + path).c_str());
		}
		logger.close();
	}
	else {
		std::cerr << "The file '" << path << "' could not opened!\r\nExiting to system!" << std::endl;
		exit(1);
	}

	if (printToScreen) {
		std::cerr << dateTime << " " << msg << std::endl;
	}
}

/**
 * Writes the command line arguments to the log file.
 *
 * @param argc contains the number of arguments passed to the program
 * @param argv is a one-dimensional array of strings, each string is one of the arguments that was passed to the program
 * @param binaryFileAdapter file adapter to write the data to the log file
 */
void BinaryFileAdapter::LogStartParameters(int argc, char* argv[])
{
	std::string commandLine = "CommandLine:";
	for (int i=0; i<argc; i++) {
		commandLine += ' ';
		commandLine += argv[i];
	}
	Log(commandLine, false);
}

/**
 * Writes the elapsed wall time between the function call and the parameter t0
 *
 * @param msg the message to write to the log file
 * @param t0 the reference wall time from wgich the elapsed time is calculated
 */
void BinaryFileAdapter::LogTimeSpan(std::string msg, time_t t0)
{
	int diff = (int)(time(0) - t0);

	std::string diffTime;
	Tools::TimeDifference(diff, diffTime);
	msg += diffTime;
	std::cerr << msg;
	Log(msg, false);
}

/// <summary>
/// Save the phases of the bodies defined in the list parameter into the file defined
/// by the path parameter. The format of the output file is binary.
/// </summary>
/// <param name="path">The path of the output file</param>
/// <param name="list">The list of the bodies whose phases will be saved</param>
void BinaryFileAdapter::SavePhases(double time, int n, double *y, int *id)
{
	std::string path = output->GetPath(output->phases);
	std::ofstream writer(path.c_str(), std::ios::out | std::ios::app | std::ios::binary);
	if (writer) {
		writer.write(reinterpret_cast<char*>(&time), sizeof(time));
		writer.write(reinterpret_cast<char*>(&n),    sizeof(n));
		for (register int i=0; i<n; i++) {
			SavePhase(writer, &(y[6*i]), &(id[i]));
		}
		writer.close();
	}
	else {
		Log("The file '" + path + "' could not opened!", true);
		exit(1);
	}
}

void BinaryFileAdapter::SavePhase(std::ofstream& writer, double *y, int *id)
{
	writer.write(reinterpret_cast<char*>(id), sizeof(*id));
	writer.write(reinterpret_cast<char*>(y),  6*sizeof(*y));
	if (writer.bad()) {
		_errMsg = "An error occurred during writing the phase!";
		Log(_errMsg, true);
		perror(_errMsg.c_str());
		exit(1);
	}
}

/**
 * Saves the energy, the angular momentum vector and its length, the position vector of the baricenter and its length
 * and the velocity of the barycenter and its length.
 */
void BinaryFileAdapter::SaveIntegrals(double time, int n, double *integrals)
{
	static bool firstCall = true;

	std::string path = output->GetPath(output->integrals);
	std::ofstream writer(path.c_str(), std::ios::out | std::ios::app | std::ios::binary);
	if (writer) {
		if (firstCall)
		{
			char header[] = "time [day], mass [Solar], bary center position r: x [AU], y [AU], z [AU], bary center velocity v: x [AU/day], y [AU/day], z [AU/day], length of r [AU], length of v [AU/day], angular momentum vector c: x, y, z, length of c, kinetic energy, potential energy, total energy";
			int len = (int)strlen(header);
			writer.write(reinterpret_cast<char*>(&len), sizeof(int));
			writer.write(header, len);
			firstCall = false;
		}
		int nElement = n + 1;
		writer.write(reinterpret_cast<char*>(&nElement), sizeof(nElement));
		writer.write(reinterpret_cast<char*>(&time), sizeof(time));
		writer.write(reinterpret_cast<char*>(integrals), n * sizeof(*integrals));

		if (writer.bad()) {
			_errMsg = "An error occurred during writing the integrals!";
			Log(_errMsg, true);
			perror(_errMsg.c_str());
			exit(1);
		}
	}
	else {
		Log("The file '" + path + "' could not opened!", true);
		exit(1);
	}

	writer.close();
}

/// <summary>
/// Writes the list of TwoBodyAffairs data to the disk.
/// </summary>
/// <param name="path">The path of the output file</param>
/// <param name="list">The data of the TwoBodyAffairs</param>
void BinaryFileAdapter::SaveTwoBodyAffairs(std::list<TwoBodyAffair>& list)
{
	std::string path = output->GetPath(output->twoBodyAffair);
	std::ofstream writer(path.c_str(), std::ios::out | std::ios::app | std::ios::binary);
	if (writer) {
		for (std::list<TwoBodyAffair>::iterator it = list.begin(); it != list.end(); it++) {
			SaveTwoBodyAffair(writer, *it);
		}
		writer.close();
	}
	else {
		_errMsg = "The file '" + path + "' could not opened!";
		Log(_errMsg, true);
		perror(_errMsg.c_str());
		exit(1);
	}
}

void BinaryFileAdapter::SaveTwoBodyAffair(std::ofstream& writer, TwoBodyAffair& affair)
{
	writer.write(reinterpret_cast<char*>(&(affair.id)), sizeof(int));
	writer.write(reinterpret_cast<char*>(&(affair.type)), sizeof(int));

	writer.write(reinterpret_cast<char*>(&(affair.body1Id)), sizeof(int));
	writer.write(reinterpret_cast<char*>(&(affair.body2Id)), sizeof(int));
	writer.write(reinterpret_cast<char*>(&(affair.body1Phase)), 6*sizeof(double));
	writer.write(reinterpret_cast<char*>(&(affair.body2Phase)), 6*sizeof(double));
	writer.write(reinterpret_cast<char*>(&(affair.time)), sizeof(double));

	if (writer.bad()) {
		_errMsg = "An error occurred during writing the two body affair!";
		Log(_errMsg, true);
		perror(_errMsg.c_str());
		exit(1);
	}
}

void BinaryFileAdapter::SaveBodyProperties(double time, std::list<Body* >& bodyList)
{
	std::string constPropPath = output->GetPath(output->constantProperties);
	std::ofstream constPropWriter(constPropPath.c_str(), std::ios::out | std::ios::app | std::ios::binary);
	if (!constPropWriter) {
		_errMsg = "The file '" + constPropPath + "' could not opened!";
		Log(_errMsg, true);
		perror(_errMsg.c_str());
		exit(1);
	}

	std::string varPropPath = output->GetPath(output->variableProperties);
	std::ofstream varPropWriter(varPropPath.c_str(), std::ios::out | std::ios::app | std::ios::binary);
	if (!varPropWriter) {
		_errMsg = "The file '" + varPropPath + "' could not opened!";
		Log(_errMsg, true);
		perror(_errMsg.c_str());
		constPropWriter.close();
		exit(1);
	}

	for (std::list<Body* >::iterator it = bodyList.begin(); it != bodyList.end(); it++) {
		SaveConstantProperty(constPropWriter, *it);
		SaveVariableProperty(varPropWriter, *it, time);
	}
	constPropWriter.close();
	varPropWriter.close();
}

void BinaryFileAdapter::SaveConstantProperty(Body* body)
{
	std::string constPropPath = output->GetPath(output->constantProperties);
	std::ofstream constPropWriter(constPropPath.c_str(), std::ios::out | std::ios::app | std::ios::binary);
	if (!constPropWriter) {
		_errMsg = "The file '" + constPropPath + "' could not opened!";
		Log(_errMsg, true);
		perror(_errMsg.c_str());
		exit(1);
	}
	SaveConstantProperty(constPropWriter, body);
	constPropWriter.close();
}

void BinaryFileAdapter::SaveConstantProperty(std::ofstream& writer, Body* body)
{
	int i = body->GetId();
	writer.write(reinterpret_cast<char *>(&i), sizeof(i));

	std::string s = body->guid;
	char *guid = new char[16];
	Tools::GuidToCharArray(s, guid);

	// NOTE: this must be done in such a way that when the C# initializes a new guid form an array of 16 bytes, the correct
	// guid will be created. The constructor : Guid(int, short, short, byte, byte, byte, byte, byte, byte, byte, byte)
	// This will be the int part
	for (int j=3; j>=0; j--) writer.write(&guid[j], 1);
	// This will be the 1. short part
	for (int j=5; j>=4; j--) writer.write(&guid[j], 1);
	// This will be the 2. short part
	for (int j=7; j>=6; j--) writer.write(&guid[j], 1);
	// This will be the last 8 bytes
	writer.write(&guid[8], 8);
	//writer.write(guid, 16);

	// NOTE: The length of the strings cannot be more than 127!
	// See: http://stackoverflow.com/questions/1550560/encoding-an-integer-in-7-bit-format-of-c-sharp-binaryreader-readstring
	// because I will read this file with C#, and BinaryReader.ReadString method use a 7-bit encoding for the length of the string.
	const char *p = body->name.c_str();
	char len = (char)strlen(p);
	writer.write(&len, sizeof(char));
	writer.write(p, len);

	p = body->designation.c_str();
	len = (char)strlen(p);
	writer.write(&len, sizeof(char));
	writer.write(p, len);

	p = body->provisionalDesignation.c_str();
	len = (char)strlen(p);
	writer.write(&len, sizeof(char));
	writer.write(p, len);

	p = body->reference.c_str();
	len = (char)strlen(p);
	writer.write(&len, sizeof(char));
	writer.write(p, len);

	p = body->opposition.c_str();
	len = (char)strlen(p);
	writer.write(&len, sizeof(char));
	writer.write(p, len);

	i = body->ln;
	writer.write(reinterpret_cast<char *>(&i), sizeof(i));

	i = body->type;
	writer.write(reinterpret_cast<char *>(&i), sizeof(i));

	i = body->mPCOrbitType;
	writer.write(reinterpret_cast<char *>(&i), sizeof(i));

	i = body->migrationType;
	writer.write(reinterpret_cast<char *>(&i), sizeof(i));

	if (body->characteristics != 0) {
		double d = body->characteristics->absVisMag;
		writer.write(reinterpret_cast<char *>(&d), sizeof(d));

		d = body->characteristics->stokes;
		writer.write(reinterpret_cast<char *>(&d), sizeof(d));
	}
	else {
		double d = 0.0;
		writer.write(reinterpret_cast<char *>(&d), sizeof(d));
		writer.write(reinterpret_cast<char *>(&d), sizeof(d));
	}

	if (writer.bad()) {
		_errMsg = "An error occurred during writing the constant property!";
		Log(_errMsg, true);
		perror(_errMsg.c_str());
		exit(1);
	}
}

void BinaryFileAdapter::SaveVariableProperty(Body* body, double time)
{
	std::string varPropPath = output->GetPath(output->variableProperties);
	std::ofstream varPropWriter(varPropPath.c_str(), std::ios::out | std::ios::app | std::ios::binary);
	if (!varPropWriter) {
		_errMsg = "The file '" + varPropPath + "' could not opened!";
		Log(_errMsg, true);
		perror(_errMsg.c_str());
		exit(1);
	}
	SaveVariableProperty(varPropWriter, body, time);
	varPropWriter.close();
}

void BinaryFileAdapter::SaveVariableProperty(std::ofstream& writer, Body* body, double time)
{
	int id = _propertyId++;
	writer.write(reinterpret_cast<char *>(&id), sizeof(id));

	int bodyId = body->GetId();
	writer.write(reinterpret_cast<char *>(&bodyId), sizeof(bodyId));

	double d = time;
	writer.write(reinterpret_cast<char *>(&d), sizeof(d));

	if (body->characteristics != 0) {
		d = body->characteristics->mass;
		writer.write(reinterpret_cast<char *>(&d), sizeof(d));

		d = body->characteristics->radius;
		writer.write(reinterpret_cast<char *>(&d), sizeof(d));

		d = body->characteristics->density;
		writer.write(reinterpret_cast<char *>(&d), sizeof(d));
		if (body->characteristics->componentList.size() > 0) {
			std::string compPropPath = output->GetPath(output->compositionProperties);
			std::ofstream compPropWriter(compPropPath.c_str(), std::ios::out | std::ios::app | std::ios::binary);
			if (!compPropWriter) {
				_errMsg = "The file '" + compPropPath + "' could not opened!";
				Log(_errMsg, true);
				perror(_errMsg.c_str());
				exit(1);
			}
			SaveCompositionProperty(compPropWriter, body, id);
			compPropWriter.close();
		}
	}
	else {
		double d = 0.0;
		writer.write(reinterpret_cast<char *>(&d), sizeof(d));
		writer.write(reinterpret_cast<char *>(&d), sizeof(d));
		writer.write(reinterpret_cast<char *>(&d), sizeof(d));
		if (writer.bad()) {
			_errMsg = "An error occurred during writing the variable property!";
			Log(_errMsg, true);
			perror(_errMsg.c_str());
			exit(1);
		}
	}
}

void BinaryFileAdapter::SaveCompositionProperty(std::ofstream& writer, Body* body, int propertyId)
{
	for (std::list<Component>::iterator it = body->characteristics->componentList.begin(); 
		it != body->characteristics->componentList.end(); it++) {

		int id = _compositionId++;
		writer.write(reinterpret_cast<char *>(&id), sizeof(id));

		id = propertyId;
		writer.write(reinterpret_cast<char *>(&id), sizeof(id));

		const char *p = it->name.c_str();
		char len = (char)strlen(p);
		writer.write(&len, sizeof(char));
		writer.write(p, len);

		double d = it->ratio;
		writer.write(reinterpret_cast<char *>(&d), sizeof(d));
		if (writer.bad()) {
			_errMsg = "An error occurred during writing the composition property!";
			Log(_errMsg, true);
			perror(_errMsg.c_str());
			exit(1);
		}
	}
}
