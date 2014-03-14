// includes, system 
#include <ctime>
#include <iomanip>
#include <iostream>
#include <fstream>

// includes CUDA
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

// includes Thrust
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/reduce.h>
#include <thrust/transform.h>

#include "config.h"
#include "Constants.h" 
#include "gas_disk.h"
#include "nbody.h"
#include "nbody_exception.h"
#include "ode.h"
#include "options.h"

using namespace std;

static cudaError_t HandleError(cudaError_t cudaStatus, const char *file, int line)
{
    if (cudaSuccess != cudaStatus) {
        printf( "%s in %s at line %d\n", cudaGetErrorString( cudaStatus ), file, line );
        return cudaStatus;
    }
	return cudaStatus;
}
#define HANDLE_ERROR(cudaStatus) (HandleError(cudaStatus, __FILE__, __LINE__))


__global__
void print_gas_disc(gas_disk *gasDisk)
{
	printf("eta: %10lf, %10lf\n", gasDisk->eta.x, gasDisk->eta.y);
	printf("rho: %10lf, %10lf\n", gasDisk->rho.x, gasDisk->rho.y);
	printf("sch: %10lf, %10lf\n", gasDisk->sch.x, gasDisk->sch.y);
	printf("tau: %10lf, %10lf\n", gasDisk->tau.x, gasDisk->tau.y);
}

cudaError_t unit_test_cpy_gas_disc_to_dev()
{
	cudaError_t cudaStatus = cudaSuccess;

	bool	succeeded = true;
	char	func_name[256];
	char	err_msg[1024];

	{
		bool	failed = false;
		strcpy(func_name, "unit_test_cpy_gas_disc_to_dev");

		var2_t eta = {2.0e-3, 1.0/2.0	};
		var2_t rho = {1.0e-9, -11.0/4.0	};		// g / cm^3
		var2_t sch = {5.0e-2, 5.0/4.0	};
		var2_t tau = {2.0/3.0, 2.0		};
		rho.x	*= Constants::GramPerCm3ToSolarPerAu3; // M_sun / AU^3

		gas_disk*	gasDisk;
		gas_disk*	d_gasDisk;
		gasDisk = new gas_disk(rho, sch, eta, tau);

		cout << "gasDisk: " << endl;
		cout << *gasDisk;

		cudaStatus = HANDLE_ERROR(cudaMalloc((void**)&d_gasDisk, sizeof(gas_disk)));
		if (cudaStatus != cudaSuccess) {
			sprintf(err_msg, "\t%30s() function failed at line %d.", func_name, __LINE__);
			cerr << err_msg << endl;
			failed = true;
		}

		cudaStatus = HANDLE_ERROR(cudaMemcpy(d_gasDisk, gasDisk, sizeof(gas_disk), cudaMemcpyHostToDevice ));
		if (cudaStatus != cudaSuccess) {
			sprintf(err_msg, "\t%30s() function failed at line %d.", func_name, __LINE__);
			cerr << err_msg << endl;
			failed = true;
		}

		print_gas_disc<<<1,1>>>(d_gasDisk);
		cudaStatus = HANDLE_ERROR(cudaGetLastError());
		if (cudaStatus != cudaSuccess) {
			sprintf(err_msg, "\t%30s() function failed at line %d.", func_name, __LINE__);
			cerr << err_msg << endl;
			failed = true;
		}

		cudaFree(d_gasDisk);
		delete gasDisk;
	}

	return cudaStatus;
}

// a = a + b
__global__
void add_two_vector(int_t n, var_t *a, const var_t *b)
{
	int tid = blockIdx.x * blockDim.x + threadIdx.x;
	if (n > tid) {
		a[tid] += b[tid];
	}
}

cudaError_t unit_test_transform_plus()
{
	cudaError_t cudaStatus = cudaSuccess;

	bool	succeeded = true;
	char	func_name[256];
	char	err_msg[1024];

	{
		bool	failed = false;
		strcpy(func_name, "add_two_vector");

		h_var_t h_acce;
		h_var_t	h_acceGasDrag;
		h_acce.resize(10 * 4);
		h_acceGasDrag.resize(3 * 4);

		for (int i = 0; i < 10*4; i++ ) {
			h_acce[i] = 0.0;
		}

		for (int i = 0; i < 3*4; i++ ) {
			h_acceGasDrag[i] = 1.0;
		}

		d_var_t acce = h_acce;
		d_var_t	acceGasDrag = h_acceGasDrag;

		int_t n = acceGasDrag.size();
		// 1 star + 1 gp + 3 rp
		int offset = 5 * 4;

		add_two_vector<<<1, n>>>(n, (var_t*)(acce.data().get() + offset), (var_t*)acceGasDrag.data().get());
		cudaStatus = HANDLE_ERROR(cudaGetLastError());
		if (cudaStatus != cudaSuccess) {
			sprintf(err_msg, "\t%30s() function failed at line %d.", func_name, __LINE__);
			cerr << err_msg << endl;
			failed = true;
		}

		h_acce = acce;
		for (int i = 0; i < 10; i++ ) {
			int idx = 4*i;
			printf("h_acce[%d] = %10lf, %10lf, %10lf, %10lf\n", idx, h_acce[idx], h_acce[idx+1], h_acce[idx+2], h_acce[idx+3]);
		}

	}

	return cudaStatus;
}

int main(int argc, const char** argv)
{
	cudaError_t cudaStatus = cudaSuccess;
	int		result = 0;
	char	func_name[256];
	char	err_msg[1024];

	{
		strcpy(func_name, "unit_test_cpy_gas_disc_to_dev");

		cudaStatus = unit_test_cpy_gas_disc_to_dev();
		if (cudaSuccess == cudaStatus) {
			sprintf(err_msg, "The unit test(s) of the %s() function passed.", func_name);
			cout << endl << err_msg << endl;
		}
		else {
			sprintf(err_msg, "The unit test(s) of the %s() function failed.", func_name);
			cout << endl << err_msg << endl;
		}
	}

	{
		strcpy(func_name, "unit_test_transform_plus");

		cudaStatus = unit_test_transform_plus();
		if (cudaSuccess == cudaStatus) {
			sprintf(err_msg, "The unit test(s) of the %s() function passed.", func_name);
			cout << endl << err_msg << endl;
		}
		else {
			sprintf(err_msg, "The unit test(s) of the %s() function failed.", func_name);
			cout << endl << err_msg << endl;
		}
	}

	return result;
}
