#pragma once

#include "cuda_runtime.h"

#include <exception>
#include <string>

using namespace std;

class nbody_exception : public exception
{
private:
	string message;
	cudaError_t cuda_error;
public:
	nbody_exception(string message);
	nbody_exception(string message, cudaError_t cuda_error);
	~nbody_exception() throw();

	const char* what() const throw();
};