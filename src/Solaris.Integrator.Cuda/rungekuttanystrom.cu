#include "rungekuttanystrom.h"
#include "integrator_exception.h"
#include <math.h>
#include "util.h"

#define sQ sqrt(21.0)

template<> var_t rungekuttanystrom<9>::a[] = { 1.0/200.0,
								 1.0/150.0, 1.0/75.0,
								 171.0/8192.0, 45.0/4096.0, 315.0/8192.0,
								 5.0/288.0, 25.0/528.0, 25.0/672.0, 16.0/693.0,
								 (1003.0-205.0*sQ)/12348.0, -25.0*(751.0-173.0*sQ)/90552.0, 25.0*(624.0-137.0*sQ)/43218.0, -128.0*(361.0-79.0*sQ)/237699.0, (3411.0-745.0*sQ)/24696.0,
								 (793.0+187.0*sQ)/12348.0, -25.0*(331.0+113.0*sQ)/90552.0, 25.0*(1044.0+247.0*sQ)/43218.0, -128.0*(14885.0+3779.0*sQ)/9745659.0, (3327.0+797.0*sQ)/24696.0, -(581.0+127.0*sQ)/1722.0,
								 -(157.0-3.0*sQ)/378.0, 25.0*(143.0-10.0*sQ)/2772.0, -25.0*(876.0+55.0*sQ)/3969.0, 1280.0*(913.0+18.0*sQ)/596673.0, -(1353.0+26.0*sQ)/2268.0, 7.0*(1777.0+377.0*sQ)/4428.0, 7.0*(5.0-sQ)/36.0,
								 1.0/20.0, 0, 0, 0, 8.0/45.0, 7.0*(7.0+sQ)/360.0, 7.0*(7.0-sQ)/360.0, 0 };
template<> var_t rungekuttanystrom<9>::bh[] = { 1.0/20.0, 0.0, 0.0, 0.0, 8.0/45.0, 7.0*(7.0+sQ)/360.0, 7.0*(7.0-sQ)/360.0, -1.0/20.0, 1.0/20.0 };
template<> var_t rungekuttanystrom<9>::b[] = { 1.0/20.0, 0.0, 0.0, 0.0, 8.0/45.0, 7.0*(7.0+sQ)/360.0, 7.0*(7.0-sQ)/360.0, 0.0, 0.0 };
template<> var_t rungekuttanystrom<9>::bd[] = { 1.0/20.0, 0.0, 0.0, 0.0, 16.0/45.0, 49.0/180.0, 49.0/180.0, 1.0/20.0, 0.0 };
template<> ttt_t rungekuttanystrom<9>::c[] = { 0.0, 1.0/10.0, 1.0/5.0, 3.0/8.0, 1.0/2.0, (7.0-sQ)/14.0, (7.0+sQ)/14.0, 1.0, 1.0 };
template<> var_t rungekuttanystrom<9>::err = 1 / 20.0;

template <int RKOrder>
rungekuttanystrom<RKOrder>::rungekuttanystrom(ode& f, ttt_t dt, bool adaptive, var_t tolerance) :
	integrator(f, dt),
	adaptive(adaptive),
	tolerance(tolerance),
	d_k(f.get_order()),
	d_ytemp(f.get_order())
{
	int forder = f.get_order();

	if (forder != 2)
	{
		throw integrator_exception("Only second order equations are supported.");
	}

	size_t size = f.h_y[0].size();

	d_k.resize(RKOrder);
	for (int r = 0; r < RKOrder; r++)
	{
		d_k[r].resize(size);
	}
	d_ytemp[0].resize(size);
	d_yhtemp.resize(size);
	d_dytemp.resize(size);
	d_err.resize(size);
}

template <int RKOrder>
rungekuttanystrom<RKOrder>::~rungekuttanystrom()
{
}

template <int RKOrder>
ttt_t rungekuttanystrom<RKOrder>::step()
{
	// Symbols used: (in the paper)
	// d_y[0]: coordinates (X)
	// d_y[1]: velocities  (X')
	// d_ytemp: intermediate coordinates (X)
	// d_yhtemp: intermediate coordinates (X^)
	// d_dytemp: intermediate coordinates (X'^)
	// d_k: intermediate differentials (f)

	var_t max_err;
	ttt_t dttemp;

	do
	{
		ttt_t ttemp;
		int rr = 0;

		// Calculate the acelleration matrix first
		for (int r = 0; r < RKOrder; r++)
		{
			ttemp = f.t + rungekuttanystrom::c[r] * dt;
		
			sum_vec(d_ytemp[0], f.d_y[0], f.d_y[1], (var_t)(c[r] * dt));

			// Calculate temporary values of x
			for (int s = 0; s < r; s++)
			{
				sum_vec(d_ytemp[0], d_ytemp[0], d_k[s], a[rr] * (var_t)(dt * dt));
				rr++;
			}

			// Only calculate acceleration here, not velocities
			f.calculate_dy(1, r, ttemp, f.d_p, d_ytemp, d_k[r]);
		}

		// Now the acceleration matrix is done, update variables and
		// determine error
		sum_vec(d_ytemp[0], f.d_y[0], f.d_y[1], (var_t)dt);
		copy_vec(d_dytemp, f.d_y[1]);
		if (adaptive)
		{
			sum_vec(d_yhtemp, f.d_y[0], f.d_y[1], (var_t)dt);
		}

		// Sum up orders
		for (int r = 0; r < RKOrder; r++)
		{
			sum_vec(d_ytemp[0], d_ytemp[0], d_k[r], b[r] * (var_t)(dt * dt));
			sum_vec(d_dytemp, d_dytemp, d_k[r], bd[r] * (var_t)dt);
			if (adaptive)
			{
				sum_vec(d_yhtemp, d_yhtemp, d_k[r], bh[r] * (var_t)(dt * dt));
			}
		}

		// Save current step length
		dttemp = dt;

		if (adaptive)
		{
			absdiff_vec(d_err, d_k[RKOrder - 2], d_k[RKOrder - 1], (var_t)(dt * dt) * err);

			// Find max of error and calculate new dt estimate
			max_err = max_vec(d_err);

			// Update step size
			dt *= 1.0e-15 > max_err ? 2.0 : 0.9 * pow(tolerance / max_err, (var_t)1.0 / RKOrder);
		}
	}
	while (adaptive && max_err > tolerance);

	// Time step is now accurate enough
	// Propagate variables
	copy_vec(f.d_yout[0], d_ytemp[0]);
	copy_vec(f.d_yout[1], d_dytemp);

	// Step time
	f.tout = f.t + dttemp;

	f.swap_in_out();

	return dt;
}

template class rungekuttanystrom<9>;
