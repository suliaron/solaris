// includes CUDA
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include "integrator_exception.h"
#include "opt_midpoint_method.h"

#define THREADS_PER_BLOCK	256

static cudaError_t HandleError(cudaError_t cudaStatus, const char *file, int line)
{
    if (cudaSuccess != cudaStatus) {
        printf( "%s in %s at line %d\n", cudaGetErrorString( cudaStatus ), file, line );
        return cudaStatus;
    }
	return cudaStatus;
}
#define HANDLE_ERROR(cudaStatus) (HandleError(cudaStatus, __FILE__, __LINE__))


// ytemp = y0 + a21 * k1
static __global__
void calc_2nd_arg_of_k2_kernel(int_t n, var_t *ytemp, const var_t *y0, const var_t *k1, var_t a21)
{
	int tid = blockIdx.x * blockDim.x + threadIdx.x;

	if (n > tid) {
		ytemp[tid] = y0[tid] + a21 * k1[tid];
	}
}

// y_n+1 = y_n + k2
static __global__
void calc_final_dependent_variables_kernel(int_t n, var_t *y, const var_t *y0, const var_t *k2)
{
	int tid = blockIdx.x * blockDim.x + threadIdx.x;

	if (n > tid) {
		y[tid] = y0[tid] + k2[tid];
	}
}

var_t opt_midpoint_method::a[] = {1.0/2.0};
var_t opt_midpoint_method::b[] = {0.0, 1.0};
ttt_t opt_midpoint_method::c[] = {0.0, 1.0/2.0};


void opt_midpoint_method::calculate_grid(int nData, int threads_per_block)
{
	int	nThread = std::min(threads_per_block, nData);
	int	nBlock = (nData + nThread - 1)/nThread;
	grid.x  = nBlock;
	block.x = nThread;
}

void opt_midpoint_method::call_calc_2nd_arg_of_k2_kernel()
{
	for (int i = 0; f.get_order(); i++) {
		var_t *y0 = f.d_y[i].data().get();
		var_t *k1 = d_k[i][0].data().get();

		calculate_grid(f.d_y[i].size(), THREADS_PER_BLOCK);
		calc_2nd_arg_of_k2_kernel<<<grid, block>>>(f.d_y[i].size(), d_ytemp[i].data().get(), y0, k1, a[0] * dt);
		cudaError cudaStatus = HANDLE_ERROR(cudaGetLastError());
		if (cudaSuccess != cudaStatus) {
			throw integrator_exception("calc_2nd_arg_of_k2_kernel failed");
		}
	}
}

void opt_midpoint_method::call_calc_final_dependent_variables_kernel()
{
	for (int i = 0; f.get_order(); i++) {
		var_t *y0 = f.d_yout[i].data().get();
		var_t *k2 = d_k[i][1].data().get();

		calculate_grid(f.d_y[i].size(), THREADS_PER_BLOCK);
		calc_final_dependent_variables_kernel<<<grid, block>>>(f.d_y[i].size(), f.d_yout[i].data().get(), y0, k2);
		cudaError cudaStatus = HANDLE_ERROR(cudaGetLastError());
		if (cudaSuccess != cudaStatus) {
			throw integrator_exception("calc_2nd_arg_of_k2_kernel failed");
		}
	}
}

opt_midpoint_method::opt_midpoint_method(ode& f, ttt_t dt, bool adaptive, var_t tolerance) :
		integrator(f, dt),
		adaptive(adaptive),
		tolerance(tolerance),
		d_k(f.get_order()),
		d_ytemp(f.get_order(), d_var_t())
{
	RKOrder = 2;
	int forder = f.get_order();

	for (int i = 0; i < forder; i++) {
		d_ytemp[i].resize(f.d_y[i].size());
		d_k[i].resize(RKOrder);
		for (int r = 0; r < RKOrder; r++) {
			d_k[i][r].resize(f.d_y[i].size());
		} 
	}
}

ttt_t	opt_midpoint_method::step()
{
	cudaError cudaStatus = cudaSuccess;

	int	forder = f.get_order();
	int r = 0;

	ttt_t ttemp = f.t + c[r] * dt;
	for (int i = 0; i < forder; i++) {
		f.calculate_dy(i, r, ttemp, f.d_p, f.d_y, d_k[i][r]);
	}

	r++;
	call_calc_2nd_arg_of_k2_kernel();
	ttemp = f.t + c[r] * dt;
	for (int i = 0; i < forder; i++) {
		f.calculate_dy(i, r, ttemp, f.d_p, d_ytemp, d_k[i][r]); 
	}

	call_calc_final_dependent_variables_kernel();
	f.tout = f.t + dt;
	f.swap_in_out();

	return dt;
}
