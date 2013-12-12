#include "nbody_exception.h"
#include "cuda_runtime.h"

nbody_exception::nbody_exception(string message) :
	message(message),
	cuda_error(cudaSuccess)
{
}

nbody_exception::nbody_exception(string message, cudaError_t cuda_error) :
	message(message),
	cuda_error(cuda_error)
{
}

nbody_exception::~nbody_exception() throw()
{
}

const char* nbody_exception::what() const throw()
{
	return message.c_str();
}