#include "ode.h"

ode::ode(int order, ttt_t t) :
	t(t),
	h_y(order, h_var_t()),
	d_y(order, d_var_t()),
	d_yout(order, d_var_t())
{
}

ode::~ode()
{
}

int ode::get_order()
{
	return h_y.size();
}

void ode::copy_to_host()
{
	// Copy parameters to the host
	thrust::copy(d_p.begin(), d_p.end(), h_p.begin());

	// Copy variables to the host
	for (int i = 0; i < h_y.size(); i++)
	{
		thrust::copy(d_y[i].begin(), d_y[i].end(), h_y[i].begin());
	}
}

void ode::copy_to_device()
{
	// Copy parameters to the device
	d_p.resize(h_p.size());
	thrust::copy(h_p.begin(), h_p.end(), d_p.begin());

	// Copy variables to the device
	for (int i = 0; i < h_y.size(); i++)
	{
		d_y[i].resize(h_y[i].size());
		thrust::copy(h_y[i].begin(), h_y[i].end(), d_y[i].begin());

		d_yout[i].resize(h_y[i].size());
	}
}

void ode::swap_in_out()
{
	// Swap values
	for (int i = 0; i < d_y.size(); i ++)
	{
		d_yout[i].swap(d_y[i]);
	}

	// Swap time
	std::swap(t, tout);
}