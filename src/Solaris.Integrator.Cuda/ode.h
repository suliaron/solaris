#pragma once

#include <thrust/host_vector.h>
#include <thrust/device_vector.h>

#include "config.h"

class ode
{
public:
	// Time
	ttt_t t;
	ttt_t tout;

	// Parameters
	host_var_t h_p;
	device_var_t d_p;

	// Variables
	std::vector<host_var_t> h_y;
	std::vector<device_var_t> d_y;
	std::vector<device_var_t> d_yout;

public:

	ode(int order);
	~ode();

	int get_order();

	void copy_to_host();
	void copy_to_device();

	virtual void calculate_dy(int i, ttt_t t, const device_var_t& p, const std::vector<device_var_t>& y, device_var_t& dy) = 0;
	void swap_in_out();
};