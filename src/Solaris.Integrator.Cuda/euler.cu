#include "euler.h"
#include "util.h"

euler::euler(ode& f, ttt_t dt) :
	integrator(f, dt),
	d_dy(f.get_order(), d_var_t())
{
	int forder = f.get_order();

	// Allocate storage for differentials on the device
	for (int i = 0; i < forder; i++) {
		d_dy[i].resize(f.h_y[i].size());
	}
}

euler::~euler()
{
}

ttt_t euler::step()
{
	int forder = f.get_order();

	// Calculate differentials
	for (int i = 0; i < forder; i++) {
		f.calculate_dy(i, 0, f.t, f.d_p, f.d_y, d_dy[i]);
	}

	// Propagate variables
	for (int i = 0; i < forder; i ++) {
		sum_vec(f.d_yout[i], f.d_y[i], d_dy[i], (var_t)dt);
	}
	// Propagate time
	f.tout = f.t + dt;
	f.swap_in_out();

	return dt;
}
