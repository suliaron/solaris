#include "rungekutta.h"
#include "util.h"

template<> var_t rungekutta<4>::a[] = { 1.0/2.0, 0.0, 1.0/2.0, 0.0, 0.0, 1.0 };
template<> var_t rungekutta<4>::b[] = { 1.0/6.0, 1.0/3.0, 1.0/3.0, 1.0/6.0 };
template<> ttt_t rungekutta<4>::c[] = { 0.0, 1.0/2.0, 1.0/2.0, 1.0 };

template<> var_t rungekutta<2>::a[] = { 1.0/2.0 };
template<> var_t rungekutta<2>::b[] = { 0.0, 1.0 };
template<> ttt_t rungekutta<2>::c[] = { 0.0, 1.0/2.0 };

template <int RKOrder>
rungekutta<RKOrder>::rungekutta(ode& f, ttt_t dt, bool adaptive, var_t tolerance) :
	integrator(f, dt),
	adaptive(adaptive),
	tolerance(tolerance),
	d_ytemp(f.get_order(), d_var_t()),
	d_k(f.get_order())
{
	int forder = f.get_order();

	// Allocate storage for differentials on the device	
	for (int i = 0; i < forder; i++)
	{
		d_ytemp[i].resize(f.h_y[i].size());
		d_k[i].resize(RKOrder);
		for (int r = 0; r < RKOrder; r++)
		{
			d_k[i][r].resize(f.h_y[i].size());
		}
	}
}

template <int RKOrder>
rungekutta<RKOrder>::~rungekutta()
{
}

template <int RKOrder>
ttt_t rungekutta<RKOrder>::step()
{
	int forder = f.get_order();

	int rr = 0;
	ttt_t ttemp;

	for (int r = 0; r < RKOrder; r++) {
		ttemp = f.t + rungekutta::c[r] * dt;

		for (int i = 0; i < forder; i++) {
			copy_vec(d_ytemp[i], f.d_y[i]);
		}
		
		// Calculate temporary values of the dependent variables
		for (int s = 0; s < r; s++) {
			for (int i = 0; i < forder && a[rr] != 0.0; i ++) {
				sum_vec(d_ytemp[i], d_ytemp[i], d_k[i][s], (var_t)(a[rr] * dt));
			}
			rr++;
		}

		for (int i = 0; i < forder; i++) {
			f.calculate_dy(i, r, ttemp, f.d_p, d_ytemp, d_k[i][r]);
		}
	}

	// Advance dependent variables
	for (int i = 0; i < forder; i++) {
		copy_vec(f.d_yout[i], f.d_y[i]);
		for (int r = 0; r < RKOrder; r++) {
			if (0.0 == b[r]) {
				continue;
			}
			sum_vec(f.d_yout[i], f.d_yout[i], d_k[i][r], (var_t)(b[r] * dt));
		}
	}
	// Advance time
	f.tout = f.t + dt;
	f.swap_in_out();

	return dt;
}

template class rungekutta<2>;
template class rungekutta<4>;
