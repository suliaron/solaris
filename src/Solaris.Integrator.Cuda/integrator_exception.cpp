#include "integrator_exception.h"

integrator_exception::integrator_exception(string message) :
	message(message)
{
}

integrator_exception::~integrator_exception() throw()
{
}

const char* integrator_exception::what() const throw()
{
	return message.c_str();
}