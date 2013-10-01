#define _CRT_SECURE_NO_WARNINGS
#include <algorithm>
#include <cctype>
#include <cmath>
#include <iostream>
#include <list>
#include <string>
#include <string.h>
#include <stdio.h>

#include "Tools.h"
#include "Component.h"
#include "Constants.h"
#include "Phase.h"
#include "Error.h"

enum OS {
	Windows,
	Linux
} Os;

void  Tools::ToPhase(double *y, Phase *phase)
{
	phase->position.x = y[0];
	phase->position.y = y[1];
	phase->position.z = y[2];

	phase->velocity.x = y[3];
	phase->velocity.y = y[4];
	phase->velocity.z = y[5];
}

void Tools::CheckAgainstSmallestNumber(const int n, double *y)
{
	for (int i=0; i<n; i++)
	{
		if (fabs(y[i]) < Constants::SmallestNumber)
			y[i] = 0.0;
	}
}

bool Tools::IsNumber(const std::string& str)
{
   for (size_t i = 0; i < str.length(); i++) {
	   if (std::isdigit(str[i]) || str[i] == 'e' || str[i] == 'E' || str[i] == '.' || str[i] == '-' || str[i] == '+')
           continue;
	   else
		   return false;
   }
   return true;
}

void Tools::TrimRight(std::string& str)
{
	// trim trailing spaces
	size_t endpos = str.find_last_not_of(" \t");
	if (std::string::npos != endpos ) {
		str = str.substr( 0, endpos+1 );
	}
}

void Tools::TrimLeft(std::string& str)
{
	// trim leading spaces
	size_t startpos = str.find_first_not_of(" \t");
	if (std::string::npos != startpos ) {
		str = str.substr( startpos );
	}
}

void Tools::Trim(std::string& str)
{
	TrimRight(str);
	TrimLeft(str);
}

void Tools::GuidToCharArray(std::string& guid, char *result)
{
	const char *s = guid.c_str();
	char *g = result;
	while (*s) {
		if (*s == '-') {
			s++;
			continue;
		}
		char hexByte[3] = {*s, *(s+1), 0};
		sscanf(hexByte, "%x", g);
		s += 2;
		g++;
	}
}

void Tools::SplitPath(const std::string path, const char directorySeparator, std::string &directory, std::string &fileName)
{
	size_t found;
	if ( path.find_first_of(directorySeparator) >= path.length() ) {
		fileName = path;
		return;
	}
	found = path.find_last_of(directorySeparator);
	directory = path.substr(0, found);
	fileName = path.substr(found+1);
}

int Tools::GetDirectorySeparator(char *c) {
	if (GetOs() == 1) {
		char *os = getenv("OS");
		if (os == 0) {
			std::cerr << "Missing OS environment parameter.\n";
		}
		else {
			std::cerr << "Unrecognized operating system. OS=" << os << "\n";
		}
		std::cerr << "Please enter the directory separator character: ";
		std::cin >> *c;
	}
	else {
		switch (Os) {
		case Windows:
			*c = '\\';
			break;
		case Linux:
			*c = '/';
			break;
		default:
			std::cerr << "Failed to determine the operating system type.\n";
			std::cerr << "Please enter the directory separator character: ";
			std::cin >> *c;
			break;
		}
	}

	return 0;
}

std::string Tools::GetWorkingDirectory()
{
	std::string wd;

	if (GetOs() == 1) {
		char *os = getenv("OS");
		if (os == 0) {
			std::cerr << "Missing OS environment parameter.\n";
		}
		else {
			std::cerr << "Unrecognized operating system. OS=" << os << "\n";
		}
		std::cerr << "Please enter the working directory: ";
		std::cin >> wd;
		return wd;
	}
	else {
		char *c = 0;
		switch (Os) {
		case Windows:
			// Note: This is not functioning!
			c = getenv("PATH");
			if (c == 0) {
				std::cerr << "Missing CD environment parameter.\n";
				std::cerr << "Please enter the working directory: ";
				std::cin >> wd;
				return wd;
			}
			else {
				wd = std::string(c);
				return wd;
			}
			break;
		case Linux:
			c = getenv("PWD");
			if (c == 0) {
				std::cerr << "Missing PWD environment parameter.\n";
				std::cerr << "Please enter the working directory: ";
				std::cin >> wd;
				return wd;
			}
			else {
				wd = std::string(c);
				return wd;
			}
			break;
		default:
			std::cerr << "Failed to determine the operating system type.\n";
			std::cerr << "Please enter the working directory: ";
			std::cin >> wd;
			return wd;
			break;
		}
	}
	std::cerr << "Please enter the working directory: ";
	std::cin >> wd;
	return wd;
}

int	Tools::GetOs()
{
	// Windows
	char *os = getenv("OS");
	if (os == 0 ) {
		// Linux
		os = getenv("OSTYPE");
	}
	if (os == 0 ) {
		Error::_errMsg = "Could not find the environmental variable for the operating system!";
		Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
		return 1;
	}

	char *p = os;
	do {
		*p = (char)tolower(*p);
	}while ( *p++ != 0 );


	if (	 strcmp(os, "windows_nt") == 0) Os = Windows;
	else if (strcmp(os, "linux")	  == 0) Os = Linux;
	else {
		Error::_errMsg = "Could not determine the operating system!";
		Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
		return 1;
	}
	return 0;
}

int	Tools::CreatePath(const std::string &directory, const std::string &fileName, const char separator, char **path)
{
	// Add 2 to length: 1 for the separator char and 1 for the string terminating char
	int length = directory.length() + fileName.length() + 2;
	*path = (char *)malloc(length);
	if (*path == 0) {
		Error::_errMsg = "host memory allocation";
		Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
		return 1;
	}
	strcpy(*path, directory.c_str());
	char sep[2] = {separator, '\0'};
	strcat(*path, sep);
	strcat(*path, fileName.c_str());

	return 0;
}

void Tools::TimeDifference(int diff, std::string& result)
{
	int day = diff / (24*3600);
	diff -= day*(24*3600);
	int hour = diff / 3600;
	diff -= hour * 3600;
	int minute = diff / 60;
	int second = diff - minute * 60;
	char buffer[50];
	memset(buffer, 0, 50*sizeof(char));
	sprintf(buffer, "%03dd:%02dh:%02dm:%02ds", day, hour, minute, second);
	result = buffer;
}

int Tools::MergeComponentList(std::list<Component> *list1, std::list<Component> *list2, std::list<Component> *result)
{
	if (list1->size() == 0 && list2->size() == 0) {
        return 0;
    }
	if (list1->size() == 0) {
		// TODO: the assignment copies the content or just the address???
        *result = *list2;
        return 0;
    }
	if (list2->size() == 0) {
		// TODO: the assignment copies the content or just the address???
        *result = *list1;
        return 0;
    }

	std::list<Component>::iterator it1 = result->begin();
	result->insert(it1, list1->begin(), list1->end());

	for (std::list<Component>::iterator it2 = list2->begin(); it2 != list2->end(); it2++) {
		if (Tools::Contains(result, *it2) == false) {
			result->push_back(*it2);
		}
	}

	return 0;
}

bool Tools::Contains(std::list<Component> *list, Component &component)
{
	std::string name = component.name;
	std::transform(name.begin(), name.end(), name.begin(), ::tolower);
	for (std::list<Component>::iterator it = list->begin(); it != list->end(); it++) {
		std::string n = it->name;
		std::transform(n.begin(), n.end(), n.begin(), ::tolower);
		if (n == name)
			return true;
	}
	return false;
}

Component* Tools::FindByName(std::list<Component> *list, std::string name)
{
	std::transform(name.begin(), name.end(), name.begin(), ::tolower);
	for (std::list<Component>::iterator it = list->begin(); it != list->end(); it++) {
		std::string n = it->name;
		std::transform(n.begin(), n.end(), n.begin(), ::tolower);
		if (n == name)
			return &(*it);
	}
	return 0;
}
