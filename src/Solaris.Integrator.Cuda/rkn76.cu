// includes system 
#include <sstream>      // std::ostringstream

// include CUDA
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

// include project
#include "integrator_exception.h"
#include "rkn76.h"
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

#define	LAMBDA	1.0/20.0
#define sQ sqrt(21.0)
ttt_t rkn76::c[] = { 0.0, 1.0/10.0, 1.0/5.0, 3.0/8.0, 1.0/2.0, (7.0-sQ)/14.0, (7.0+sQ)/14.0, 1.0, 1.0 };
var_t rkn76::a[] = { 1.0/200.0,
						 1.0/150.0,                  1.0/75.0,
						 171.0/8192.0,              45.0/4096.0,                  315.0/8192.0,
						 5.0/288.0,                 25.0/528.0,                    25.0/672.0,                       16.0/693.0,
						 (1003.0-205.0*sQ)/12348.0,-25.0*(751.0-173.0*sQ)/90552.0, 25.0*(624.0-137.0*sQ)/43218.0,  -128.0*(361.0-79.0*sQ)/237699.0,      (3411.0-745.0*sQ)/24696.0,
						 (793.0+187.0*sQ)/12348.0, -25.0*(331.0+113.0*sQ)/90552.0, 25.0*(1044.0+247.0*sQ)/43218.0, -128.0*(14885.0+3779.0*sQ)/9745659.0, (3327.0+797.0*sQ)/24696.0,   -(581.0+127.0*sQ)/1722.0,
						-(157.0-3.0*sQ)/378.0,      25.0*(143.0-10.0*sQ)/2772.0,  -25.0*(876.0+55.0*sQ)/3969.0,    1280.0*(913.0+18.0*sQ)/596673.0,     -(1353.0+26.0*sQ)/2268.0,  7.0*(1777.0+377.0*sQ)/4428.0, 7.0*(5.0-sQ)/36.0,
						 1.0/20.0,                   0.0,                           0.0,                              0.0,                               8.0/45.0,                 7.0*(7.0+sQ)/360.0,           7.0*(7.0-sQ)/360.0, 0.0 };
var_t rkn76::bh[]= { 1.0/20.0, 0.0, 0.0, 0.0, 8.0/45.0, 7.0*(7.0+sQ)/360.0, 7.0*(7.0-sQ)/360.0,     0.0,    0.0 };
var_t rkn76::b[] = { 1.0/20.0, 0.0, 0.0, 0.0, 8.0/45.0, 7.0*(7.0+sQ)/360.0, 7.0*(7.0-sQ)/360.0, -LAMBDA, LAMBDA };
#undef sQ

// ytemp = y_n + dt*(a21*k1)
static __global__
void calc_ytemp_for_k2_kernel(int_t n, var_t *ytemp, const var_t *y_n, const var_t *k1, var_t k1f)
{
	int tid = blockIdx.x * blockDim.x + threadIdx.x;
	int stride = gridDim.x * blockDim.x;

	while (n > tid) {
		ytemp[tid] = y_n[tid] + k1f * k1[tid];
		tid += stride;
	}
}

// ytemp = y_n + dt*(a31*k1 + a32*k2)
static __global__
void calc_ytemp_for_k3_kernel(int_t n, var_t *ytemp, const var_t *y_n, const var_t *k1, const var_t *k2, var_t k1f, var_t k2f)
{
	int tid = blockIdx.x * blockDim.x + threadIdx.x;
	int stride = gridDim.x * blockDim.x;

	while (n > tid) {
		ytemp[tid] = y_n[tid] + k1f * k1[tid] + k2f * k2[tid];
		tid += stride;
	}
}

// ytemp = y_n + dt*(a41*k1 + a42*k2 + a43*k3)
static __global__
void calc_ytemp_for_k4_kernel(int_t n, var_t *ytemp, const var_t *y_n, const var_t *k1, const var_t *k2, const var_t *k3, var_t k1f, var_t k2f, var_t k3f)
{
	int tid = blockIdx.x * blockDim.x + threadIdx.x;
	int stride = gridDim.x * blockDim.x;

	while (n > tid) {
		ytemp[tid] = y_n[tid] + k1f * k1[tid] + k2f * k2[tid] + k3f * k3[tid];
		tid += stride;
	}
}

// ytemp = y_n + dt*(a51*k1 + a52*k2 + a53*k3 + a54*k4)
static __global__
void calc_ytemp_for_k5_kernel(int_t n, var_t *ytemp, const var_t *y_n, const var_t *k1, const var_t *k2, const var_t *k3, const var_t *k4, var_t k1f, var_t k2f, var_t k3f, var_t k4f)
{
	int tid = blockIdx.x * blockDim.x + threadIdx.x;
	int stride = gridDim.x * blockDim.x;

	while (n > tid) {
		ytemp[tid] = y_n[tid] + k1f * k1[tid] + k2f * k2[tid] + k3f * k3[tid] + k4f * k4[tid];
		tid += stride;
	}
}

// ytemp = y_n + dt*(a61*k1 + a62*k2 + a63*k3 + a64*k4 + a65*k5)
static __global__
void calc_ytemp_for_k6_kernel(int_t n, var_t *ytemp, const var_t *y_n, const var_t *k1, const var_t *k2, const var_t *k3, const var_t *k4, const var_t *k5, var_t k1f, var_t k2f, var_t k3f, var_t k4f, var_t k5f)
{
	int tid = blockIdx.x * blockDim.x + threadIdx.x;
	int stride = gridDim.x * blockDim.x;

	while (n > tid) {
		ytemp[tid] = y_n[tid] + k1f * k1[tid] + k2f * k2[tid] + k3f * k3[tid] + k4f * k4[tid] + k5f * k5[tid];
		tid += stride;
	}
}

// ytemp = y_n + dt*(a71*k1 + a72*k2 + a73*k3 + a74*k4 + a75*k5 + a76*k6)
static __global__
void calc_ytemp_for_k7_kernel(int_t n, var_t *ytemp, const var_t *y_n, const var_t *k1, const var_t *k2, const var_t *k3, const var_t *k4, const var_t *k5, const var_t *k6, var_t k1f, var_t k2f, var_t k3f, var_t k4f, var_t k5f, var_t k6f)
{
	int tid = blockIdx.x * blockDim.x + threadIdx.x;
	int stride = gridDim.x * blockDim.x;

	while (n > tid) {
		ytemp[tid] = y_n[tid] + k1f * k1[tid] + k2f * k2[tid] + k3f * k3[tid] + k4f * k4[tid] + k5f * k5[tid] + k6f * k6[tid];
		tid += stride;
	}
}

// ytemp = y_n + dt*(a81*k1 + a82*k2 + a83*k3 + a84*k4 + a85*k5 + a86*k6 + a87*k7)
static __global__
void calc_ytemp_for_k8_kernel(int_t n, var_t *ytemp, const var_t *y_n, const var_t *k1, const var_t *k2, const var_t *k3, const var_t *k4, const var_t *k5, const var_t *k6, const var_t *k7, var_t k1f, var_t k2f, var_t k3f, var_t k4f, var_t k5f, var_t k6f, var_t k7f)
{
	int tid = blockIdx.x * blockDim.x + threadIdx.x;
	int stride = gridDim.x * blockDim.x;

	while (n > tid) {
		ytemp[tid] = y_n[tid] + k1f * k1[tid] + k2f * k2[tid] + k3f * k3[tid] + k4f * k4[tid] + k5f * k5[tid] + k6f * k6[tid] + k7f * k7[tid];
		tid += stride;
	}
}

// ytemp = y_n + dt*(a91*k1 + a95*k5 + a96*k6 + a97*k7)
static __global__
void calc_ytemp_for_k9_kernel(int_t n, var_t *ytemp, const var_t *y_n, const var_t *k1, const var_t *k5, const var_t *k6, const var_t *k7, var_t k1f, var_t k5f, var_t k6f, var_t k7f)
{
	int tid = blockIdx.x * blockDim.x + threadIdx.x;
	int stride = gridDim.x * blockDim.x;

	while (n > tid) {
		ytemp[tid] = y_n[tid] + k1f * k1[tid] + k5f * k5[tid] + k6f * k6[tid] + k7f * k7[tid];
		tid += stride;
	}
}

// y = y_n + dt*(bh1*k1 + bh5*k5 + bh6*k6 + bh7*k7)
static __global__
void calc_y_kernel(int_t n, var_t *y, const var_t *y_n, const var_t *k1, const var_t *k5, const var_t *k6, const var_t *k7, var_t k1f, var_t k5f, var_t k6f, var_t k7f)
{
	int tid = blockIdx.x * blockDim.x + threadIdx.x;
	int stride = gridDim.x * blockDim.x;

	while (n > tid) {
		y[tid] = y_n[tid] + k1f * k1[tid] + k5f * k5[tid] + k6f * k6[tid] + k7f * k7[tid];
		tid += stride;
	}
}

static __global__
void calc_f8_sub_f9_kernel(int_t n, var_t* result, const var_t* f8, const var_t* f9)
{
	int tid = blockIdx.x * blockDim.x + threadIdx.x;
	int stride = gridDim.x * blockDim.x;

	while (n > tid) {
		result[tid] = f8[tid] - f9[tid];
		tid += stride;
	}
}

void rkn76::call_calc_k8_sub_k9_kernel()
{
	for (int i = 0; i < f.get_order(); i++) {
		int n		= f.d_y[i].size();
		var_t *err = d_err[i].data().get();
		var_t* k8	= d_k[i][7].data().get();
		var_t* k9	= d_k[i][8].data().get();

		calculate_grid(n, THREADS_PER_BLOCK);
		calc_f8_sub_f9_kernel<<<grid, block>>>(n, err, k8, k9);
		cudaError cudaStatus = HANDLE_ERROR(cudaGetLastError());
		if (cudaSuccess != cudaStatus) {
			throw integrator_exception("calc_f8_sub_f9_kernel failed");
		}
	}
}

void rkn76::call_calc_ytemp_for_kr_kernel(int r)
{
	int idx = 0;

	for (int i = 0; i < f.get_order(); i++) {
		int n		= f.d_y[i].size();
		calculate_grid(n, THREADS_PER_BLOCK);

		var_t* y_n= f.d_y[i].data().get();
		var_t* k1 = d_k[i][0].data().get();
		var_t* k2 = d_k[i][1].data().get();
		var_t* k3 = d_k[i][2].data().get();
		var_t* k4 = d_k[i][3].data().get();
		var_t* k5 = d_k[i][4].data().get();
		var_t* k6 = d_k[i][5].data().get();
		var_t* k7 = d_k[i][6].data().get();
		var_t* k8;
		if (adaptive) {
			k8 = d_k[i][7].data().get();
		}
		switch (r) {
		case 1:
			idx = 0;		
			calc_ytemp_for_k2_kernel<<<grid, block>>>(n, d_ytemp[i].data().get(), y_n, k1, a[idx]*dt_try);
			break;
		case 2:
			idx = 1;
			calc_ytemp_for_k3_kernel<<<grid, block>>>(n, d_ytemp[i].data().get(), y_n, k1, k2, a[idx]*dt_try, a[idx+1]*dt_try);
			break;
		case 3:
			idx = 3;
			calc_ytemp_for_k4_kernel<<<grid, block>>>(n, d_ytemp[i].data().get(), y_n, k1, k2, k3, a[idx]*dt_try, a[idx+1]*dt_try, a[idx+2]*dt_try);
			break;
		case 4:
			idx = 6;
			calc_ytemp_for_k5_kernel<<<grid, block>>>(n, d_ytemp[i].data().get(), y_n, k1, k2, k3, k4, a[idx]*dt_try, a[idx+1]*dt_try, a[idx+2]*dt_try, a[idx+3]*dt_try);
			break;
		case 5:
			idx = 10;
			calc_ytemp_for_k6_kernel<<<grid, block>>>(n, d_ytemp[i].data().get(), y_n, k1, k2, k3, k4, k5, a[idx]*dt_try, a[idx+1]*dt_try, a[idx+2]*dt_try, a[idx+3]*dt_try, a[idx+4]*dt_try);
			break;
		case 6:
			idx = 15;
			calc_ytemp_for_k7_kernel<<<grid, block>>>(n, d_ytemp[i].data().get(), y_n, k1, k2, k3, k4, k5, k6, a[idx]*dt_try, a[idx+1]*dt_try, a[idx+2]*dt_try, a[idx+3]*dt_try, a[idx+4]*dt_try, a[idx+5]*dt_try);
			break;
		case 7:
			idx = 21;
			calc_ytemp_for_k8_kernel<<<grid, block>>>(n, d_ytemp[i].data().get(), y_n, k1, k2, k3, k4, k5, k6, k7, a[idx]*dt_try, a[idx+1]*dt_try, a[idx+2]*dt_try, a[idx+3]*dt_try, a[idx+4]*dt_try, a[idx+5]*dt_try, a[idx+6]*dt_try);
			break;
		case 8:
			idx = 28;
			calc_ytemp_for_k9_kernel<<<grid, block>>>(n, d_ytemp[i].data().get(), y_n, k1, k5, k6, k7, a[idx]*dt_try, a[idx+4]*dt_try, a[idx+5]*dt_try, a[idx+6]*dt_try);
			break;
		default:
			ostringstream msg("call_calc_ytemp_for_kr_kernel() function was called with invalid parameter: ", ostringstream::ate);
			msg << r+1 << "!";
			throw integrator_exception(msg.str());
		}
		cudaError cudaStatus = HANDLE_ERROR(cudaGetLastError());
		if (cudaSuccess != cudaStatus) {
			ostringstream msg("calc_ytemp_for_k", ostringstream::ate);
			msg << r+1 << "_kernel failed";
			throw integrator_exception(msg.str());
		}
	}
}

void rkn76::call_calc_y_kernel()
{
	for (int i = 0; i < f.get_order(); i++) {
		int n		= f.d_y[i].size();
		calculate_grid(n, THREADS_PER_BLOCK);

		var_t* y_n= f.d_y[i].data().get();
		var_t *y  = f.d_yout[i].data().get();
		var_t* k1 = d_k[i][0].data().get();
		var_t* k5 = d_k[i][4].data().get();
		var_t* k6 = d_k[i][5].data().get();
		var_t* k7 = d_k[i][6].data().get();
		calc_y_kernel<<<grid, block>>>(n, y, y_n, k1, k5, k6, k7, b[0]*dt_try, b[4]*dt_try, b[5]*dt_try, b[6]*dt_try);
		cudaError cudaStatus = HANDLE_ERROR(cudaGetLastError());
		if (cudaSuccess != cudaStatus) {
			throw integrator_exception("calc_y_kernel failed");
		}
	}
}

rkn76::rkn76(ode& f, ttt_t dt, bool adaptive, var_t tolerance) :
		integrator(f, dt),
		adaptive(adaptive),
		tolerance(tolerance),
		d_k(f.get_order()),
		d_ytemp(f.get_order(), d_var_t()),
		d_err(f.get_order(), d_var_t())
{
	RKOrder = 7;
	r_max = adaptive ? RKOrder + 2 : RKOrder;
	int forder = f.get_order();

	for (int i = 0; i < forder; i++) {
		int size = f.d_y[i].size();
		d_ytemp[i].resize(size);
		if (adaptive) {
			d_err[i].resize(size);
		}
		d_k[i].resize(r_max);
		for (int r = 0; r < r_max; r++) {
			d_k[i][r].resize(size);
		}
	}
}

void rkn76::calculate_grid(int nData, int threads_per_block)
{
	int	nThread = std::min(threads_per_block, nData);
	int	nBlock = (nData + nThread - 1)/nThread;
	grid.x  = nBlock;
	block.x = nThread;
}

ttt_t rkn76::step()
{
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
		// Calculate k3 = f(tn + c3 * dt, yn + a31 * dt * k1 + ...) = d_k[][2]
		// Calculate k4 = f(tn + c4 * dt, yn + a41 * dt * k1 + ...) = d_k[][3]
		// ...
		// Calculate k7 = f(tn + c7 * dt, yn + a71 * dt * k1 + ...) = d_k[][6]
		for (r = 1; r < RKOrder; r++) {
			ttemp = f.t + c[r] * dt;
			call_calc_ytemp_for_kr_kernel(r);
			for (int i = 0; i < forder; i++) {
				f.calculate_dy(i, r, ttemp, f.d_p, d_ytemp, d_k[i][r]);
			}
		}

		dt_did = dt_try;
		if (adaptive) {
			// Calculate k8 = f(tn + c8 * dt, yn + a81 * dt * k1 + ...) = d_k[][7]
			// Calculate k9 = f(tn + c9 * dt, yn + a91 * dt * k1 + ...) = d_k[][8]
			for (r = RKOrder; r < r_max; r++) {
				ttemp = f.t + c[r] * dt;
				call_calc_ytemp_for_kr_kernel(r);
				for (int i = 0; i < forder; i++) {
					f.calculate_dy(i, r, ttemp, f.d_p, r == r_max - 1 ? f.d_yout : d_ytemp, d_k[i][r]);
				}
			}
			// calculate d_err = f8 - f9
			call_calc_k8_sub_k9_kernel();
			max_err = fabs(dt_try*LAMBDA*std::max(max_vec(d_err[0]), max_vec(d_err[1])));
			dt_try *= 0.9 * pow(tolerance / max_err, 1.0/8.0);
		}
		else {
			call_calc_y_kernel();
		}
		iter++;
	} while(adaptive && max_err > tolerance);
	n_failed_step += (iter - 1);
	n_step++;
	// Set the next step size
	dt = dt_try;

	f.tout = f.t + dt_did;
	f.swap_in_out();

	return dt_did;
}

#undef LAMBDA
