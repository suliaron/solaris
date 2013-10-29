#include "nbody_exception.h"

nbody_exception::nbody_exception(string message) :
	message(message)
{
}

nbody_exception::~nbody_exception() throw()
{
}

const char* nbody_exception::what() const throw()
{
	return message.c_str();
}