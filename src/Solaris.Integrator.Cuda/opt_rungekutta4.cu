// include CUDA
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

// include project
#include "integrator_exception.h"
#include "opt_rungekutta4.h"
#include "util.h"

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

#define	LAMBDA	1.0/10.0

ttt_t opt_rungekutta4::c[] =  {0.0, 1.0/2.0, 1.0/2.0, 1.0, 1.0};
var_t opt_rungekutta4::a[] =  {0.0, 1.0/2.0, 1.0/2.0, 1.0, 1.0/6.0, 1.0/3.0, 1.0/3.0, 1.0/6.0};
var_t opt_rungekutta4::bh[] = {1.0/6.0, 1.0/3.0, 1.0/3.0, 1.0/6.0, 0.0};
var_t opt_rungekutta4::b[] =  {1.0/6.0, 1.0/3.0, 1.0/3.0, 1.0/6.0 - LAMBDA, LAMBDA};


// ytemp = y_n + a*kr, r = 2, 3, 4
static __global__
void calc_ytemp_for_kr_kernel(int_t n, var_t *ytemp, const var_t *y_n, const var_t *kr, var_t a)
{
	int tid = blockIdx.x * blockDim.x + threadIdx.x;
	int stride = gridDim.x * blockDim.x;

	while (n > tid) {
		ytemp[tid] = y_n[tid] + a * kr[tid];
		tid += stride;
	}
}

static __global__
void calc_yHat_kernel(int_t n, var_t *y_hat, const var_t *y_n, const var_t *k1, const var_t *k2, const var_t *k3, const var_t *k4, var_t b0, var_t b1)
{
	int tid = blockIdx.x * blockDim.x + threadIdx.x;
	int stride = gridDim.x * blockDim.x;

	while (n > tid) {
		y_hat[tid] = y_n[tid] + b0 * (k1[tid] + k4[tid]) + b1 * (k2[tid] + k3[tid]);
		tid += stride;
	}
}

static __global__
void calc_k4_sub_k5_kernel(int_t n, var_t *result, const var_t *k4, const var_t* k5)
{
	int tid = blockIdx.x * blockDim.x + threadIdx.x;
	int stride = gridDim.x * blockDim.x;

	while (n > tid) {
		result[tid] = k4[tid] - k5[tid];
		tid += stride;
	}
}

void opt_rungekutta4::calculate_grid(int nData, int threads_per_block)
{
	int	nThread = std::min(threads_per_block, nData);
	int	nBlock = (nData + nThread - 1)/nThread;
	grid.x  = nBlock;
	block.x = nThread;
}

void opt_rungekutta4::call_calc_ytemp_for_kr_kernel(int r)
{
	for (int i = 0; i < f.get_order(); i++) {
		int n		= f.d_y[i].size();
		var_t *y_n	= f.d_y[i].data().get();
		var_t *kr	= d_k[i][r-1].data().get();

		calculate_grid(f.d_y[i].size(), THREADS_PER_BLOCK);
		calc_ytemp_for_kr_kernel<<<grid, block>>>(n, d_ytemp[i].data().get(), y_n, kr, a[r] * dt_try);
		cudaError cudaStatus = HANDLE_ERROR(cudaGetLastError());
		if (cudaSuccess != cudaStatus) {
			throw integrator_exception("calc_ytemp_for_kr_kernel failed");
		}
	}
}

void opt_rungekutta4::call_calc_yHat_kernel()
{
	for (int i = 0; i < f.get_order(); i++) {
		int n = f.d_y[i].size();
		var_t *y_n   = f.d_y[i].data().get();
		var_t *y_Hat = f.d_yout[i].data().get();
		var_t *k1	 = d_k[i][0].data().get();
		var_t *k2	 = d_k[i][1].data().get();
		var_t *k3	 = d_k[i][2].data().get();
		var_t *k4	 = d_k[i][3].data().get();

		calculate_grid(n, THREADS_PER_BLOCK);
		calc_yHat_kernel<<<grid, block>>>(n, y_Hat, y_n, k1, k2, k3, k4, b[0] * dt_try, b[1] * dt_try);
		cudaError cudaStatus = HANDLE_ERROR(cudaGetLastError());
		if (cudaSuccess != cudaStatus) {
			throw integrator_exception("calc_yHat_kernel failed");
		}
	}
}

void opt_rungekutta4::call_calc_k4_sub_k5_kernel()
{
	for (int i = 0; i < f.get_order(); i++) {
		int n = f.d_y[i].size();
		var_t *err = d_err[i].data().get();
		var_t *k4  = d_k[i][3].data().get();
		var_t *k5  = d_k[i][4].data().get();

		calculate_grid(n, THREADS_PER_BLOCK);
		calc_k4_sub_k5_kernel<<<grid, block>>>(n, err, k4, k5);
		cudaError cudaStatus = HANDLE_ERROR(cudaGetLastError());
		if (cudaSuccess != cudaStatus) {
			throw integrator_exception("calc_k5_sub_k4_kernel failed");
		}
	}
}

opt_rungekutta4::opt_rungekutta4(ode& f, ttt_t dt, bool adaptive, var_t tolerance) :
		integrator(f, dt),
		adaptive(adaptive),
		tolerance(tolerance),
		d_k(f.get_order()),
		d_ytemp(f.get_order(), d_var_t()),
		d_err(f.get_order(), d_var_t())
{
	RKOrder = 4;
	r_max = adaptive ? RKOrder + 1 : RKOrder;
	int	forder = f.get_order();

	for (int i = 0; i < forder; i++) {
		d_ytemp[i].resize(f.d_y[i].size());
		if (adaptive) {
			d_err[i].resize(f.d_y[i].size());
		}
		d_k[i].resize(r_max);
		for (int r = 0; r < r_max; r++) {
			d_k[i][r].resize(f.d_y[i].size());
		}
	}
}

ttt_t opt_rungekutta4::step()
{
	cudaError cudaStatus = cudaSuccess;

	int	forder = f.get_order();

	int r = 0;
	// Calculate k1 = f(tn, yn) = d_k[][0]
	ttt_t ttemp = f.t + c[r] * dt;
	for (int i = 0; i < forder; i++) {
		f.calculate_dy(i, r, ttemp, f.d_p, f.d_y, d_k[i][r]);
	}

	dt_try = dt;
	var_t max_err = 0.0;
	int_t iter = 0;
	do {
		// Calculate k2 = f(tn + c2 * dt, yn + a21 * dt * k1) = d_k[][1]
		// Calculate k3 = f(tn + c3 * dt, yn + a31 * dt * k2) = d_k[][2]
		// Calculate k4 = f(tn + c4 * dt, yn + a41 * dt * k3) = d_k[][3]
		for (r = 1; r < RKOrder; r++) {
			ttemp = f.t + c[r] * dt;
			call_calc_ytemp_for_kr_kernel(r);
			for (int i = 0; i < forder; i++) {
				f.calculate_dy(i, r, ttemp, f.d_p, d_ytemp, d_k[i][r]);
			}
		}

		// yHat_(n+1) = y_n + dt*(1/6*k1 + 1/3*k2 + 1/3*k3 + 1/6*k4) + O(dt^5)
		// f.d_yout = yHat_(n+1)
		call_calc_yHat_kernel();

		dt_did = dt_try;
		if (adaptive) {
			r = 4;
			ttemp = f.t + c[r] * dt;
			// Calculate k5 = f(tn + c5 * dt, yn + a41 * dt * k3) = d_k[][3]
			for (int i = 0; i < forder; i++) {
				f.calculate_dy(i, r, ttemp, f.d_p, f.d_yout, d_k[i][r]);
			}
			// calculate: d_err = k4 - k5
			call_calc_k4_sub_k5_kernel();
			max_err = fabs(dt_try*LAMBDA*std::max(max_vec(d_err[0]), max_vec(d_err[1])));
			dt_try *= 0.9 * pow(tolerance / max_err, 1.0/4.0);
		}
		iter++;
	} while (adaptive && max_err > tolerance);
	n_failed_step += (iter - 1);
	n_step++;
	// Set the next step size
	dt = dt_try;

	f.tout = f.t + dt_did;
	f.swap_in_out();

	return dt_did;
}

#undef LAMBDA
