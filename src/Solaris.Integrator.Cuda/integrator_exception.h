#pragma once

#include <exception>
#include <string>

using namespace std;

class integrator_exception : public exception
{
private:
	string message;
public:
	integrator_exception(string message);
	~integrator_exception() throw();

	const char* what() const throw();
};