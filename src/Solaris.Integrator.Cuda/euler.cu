#include "euler.h"
#include "util.h"

euler::euler(ode& f, ttt_t dt) :
	integrator(f, dt),
	d_dy(f.get_order(), device_var_t())
{
	// Allocate storage for differentials on the device
	for (int i = 0; i < f.h_y.size(); i ++)
	{
		d_dy[i].resize(f.h_y[i].size());
	}
}

euler::~euler()
{
}

ttt_t euler::step()
{
	int order = f.get_order();

	// Calculate differentials
	for (int i = 0; i < order; i++)
	{
		f.calculate_dy(i, 0, f.t, f.d_p, f.d_y, d_dy[i]);
	}

	// Propagate variables
	for (int i = 0; i < order; i ++)
	{
		sum_vec(f.d_yout[i], f.d_y[i], d_dy[i], (var_t)dt);
	}

	// Propagate time
	f.tout = f.t + dt;

	f.swap_in_out();

	return dt;
}
