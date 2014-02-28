#pragma once

#include <thrust/host_vector.h>
#include <thrust/device_vector.h>

#include "config.h"

class ode
{
public:
	//! Current time
	ttt_t t;
	//! Time at the end of step
	ttt_t tout;
	//! Host vector of ODE parameters
	h_var_t h_p;
	//! Device vector of ODE parameters
	d_var_t d_p;
	//! Host vector of ODE variables
	std::vector<h_var_t> h_y;
	//! Device vectors of ODE variables at the beginning of the step (at time t)
	std::vector<d_var_t> d_y;
	//! Device vectors of ODE variables at the end of the step (at time tout)
	std::vector<d_var_t> d_yout;

	//! Constructs an ordinary differential equation.
	/*!
		\param order The order of the ODE.
	*/
	ode(int order, ttt_t t);

	//! Frees up resources used by the ODE
	~ode();

	//! Returns the order of the ODE
	int get_order();

	//! Copies ODE parameters and variables from the cuda device to the host
	void copy_to_host();

	//! Copies ODE parameters and variables from the host to the cuda device
	void copy_to_device();

	//! Swaps in and out variables
	/*!
		To preserve memory, two sets of variables are used when calculating values
		of variables for the next iteration.
	*/
	void swap_in_out();

	//! Calculates the differentials of variables
	/*!
		This function is called by the integrator when calculation of the differentials is necessary
		\param i Order of the variables to calculate (e.g. 0: velocities, 1: acceleration ...)
		\param r Number of substep, used with higher order methods, like Runge-Kutta
		\param t Time
		\param p Device vector with ODE parameters
		\param y Device vector with ODE variables
		\param dy Device vector that will hold the differentials
	*/
	virtual void calculate_dy(int i, int r, ttt_t t, const d_var_t& p, const std::vector<d_var_t>& y, d_var_t& dy) = 0;
	
};