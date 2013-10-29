#pragma once

#include <exception>
#include <string>

using namespace std;

class nbody_exception : public exception
{
private:
	string message;
public:
	nbody_exception(string message);
	~nbody_exception() throw();

	const char* what() const throw();
};