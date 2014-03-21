#include <stdio.h>
#include <stdlib.h>
#include <cmath>

#include "RungeKutta4.h"
#include "Acceleration.h"
#include "BodyData.h"
#include "Error.h"
#include "TimeLine.h"

RungeKutta4::RungeKutta4()
{
	name		= "The classical Runge–Kutta method";
	reference	= "";//"http\://en.wikipedia.org/wiki/Runge-kutta#Explicit_Runge.E2.80.93Kutta_methods";

	accuracy = -10.0;
	epsilon	 = pow(10, accuracy);
}

int RungeKutta4::Driver(BodyData *bodyData, Acceleration *acceleration, TimeLine *timeLine)
{
	int	result = 0;
	int	nVar = bodyData->nBodies.NOfVar();

	bodyData->time	= timeLine->time;
	bodyData->h		= timeLine->hNext;

	acceleration->evaluateGasDrag			= true;
	acceleration->evaluateTypeIMigration	= true;
	acceleration->evaluateTypeIIMigration	= true;

	// Calculate the acceleration in the initial point
	acceleration->Compute(timeLine->time, bodyData->y0, bodyData->accel);

	//acceleration->evaluateGasDrag			= false;
	acceleration->evaluateTypeIMigration	= false;
	acceleration->evaluateTypeIIMigration	= false;

	if ((result = Step(bodyData, acceleration))== 1) {
		Error::_errMsg = "An error occurred during Runge-Kutta4 step!";
		Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
	}
	if (result == 0)
	{
		timeLine->hDid = bodyData->h;
		// Update time
		timeLine->time += timeLine->hDid;
		bodyData->time = timeLine->time;
		// Calculate the next stepsize
		timeLine->hNext = bodyData->h;
		// Update the phases of the system
		std::swap(bodyData->y0, bodyData->y);
	}

	return result;
}

int RungeKutta4::Step(BodyData *bodyData, Acceleration *acceleration)
{
	static double	a21 = 1.0/2.0;
	static double	a32 = 1.0/2.0;
	static double	a43 = 1.0;
	static double	b1 = 1.0/6.0;
	static double	b2 = 1.0/3.0;
	static double	b3 = 1.0/3.0;
	static double	b4 = 1.0/6.0;
	static double	c1 = 0.0;
	static double	c2 = 1.0/2.0;
	static double	c3 = 1.0/2.0;
	static double	c4 = 1.0;

	// These arrays will contain the accelerations computed along the trajectory of the current step
	double	*fk[4] = {0, 0, 0, 0};
	// Contains the approximation of the solution
	double	*yTemp = 0;
	int		nVar = bodyData->nBodies.NOfVar();

	// i=0 changed to i=1
	for (int i = 1; i < 4; i++) {
		fk[i] = new double[nVar];
		if (fk[i] == 0) {
			Error::_errMsg = "host memory allocation";
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			return 1;
		}
	}
	yTemp = new double[nVar];
	if (yTemp == 0) {
		Error::_errMsg = "host memory allocation";
		Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
		return 1;
	}

	double	h = bodyData->h;
	double	t = bodyData->time;

	// k1:
	fk[0] = bodyData->accel;

	// k2:
	for (int i = 0; i < nVar; i++) {
		yTemp[i] = bodyData->y0[i] + h*(a21*fk[0][i]);
	}
	acceleration->Compute(t + c2*h, yTemp, fk[1]);

	// k3:
	for (int i = 0; i < nVar; i++) {
		yTemp[i] = bodyData->y0[i] + h*(a32*fk[1][i]);
	}
	acceleration->Compute(t + c3*h, yTemp, fk[2]);

	// k4:
	for (int i = 0; i < nVar; i++) {
		yTemp[i] = bodyData->y0[i] + h*(a43*fk[2][i]);
	}
	acceleration->Compute(t + c4*h, yTemp, fk[3]);

	// The result of the step
	for (int i = 0; i < nVar; i++) {
		bodyData->y[i] = bodyData->y0[i] + h*(b1*fk[0][i] + b2*fk[1][i] + b3*fk[2][i] + b4*fk[3][i]);
	}

	// i=0 changed to i=1
	for (int i = 1; i < 4; i++) {
		delete[] fk[i];
	}
	delete[] yTemp;

	return 0;
}