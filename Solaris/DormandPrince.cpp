#include <stdio.h>
#include <stdlib.h>
#include <cmath>

#include "DormandPrince.h"
#include "Acceleration.h"
#include "BodyData.h"
#include "Error.h"
#include "TimeLine.h"

#define HANDLE_NULL(a) \
	if ((a) == NULL) { \
		Error::_errMsg = "host memory allocation"; \
		Error::PushLocation(__FILE__, __FUNCTION__, __LINE__); \
        return 1; \
	}

#define HANDLE_RESULT(result) \
	if (result == 1) { \
		Error::PushLocation(__FILE__, __FUNCTION__, __LINE__); \
		return 1; \
	}

#define CUBE(a)		((a)*(a)*(a))

DormandPrince::DormandPrince()
{
	name		= "Dormand-Prince 7(6)";
	reference	= "New Runge-Kutta Algorithms for Numerical Simulation in Dynamical Astronomy, Celestial Mechanics, Vol. 18(1978), 223-232.";

	accuracy = -10.0;
	epsilon	 = pow(10, accuracy);

	maxIter			= 10;
	sizeHeightRKD	= 9;

	double sQ = sqrt(21.0); 
       
	D1[0] = 1.0/20.0; 
	D1[1] = 0.0; 
	D1[2] = 0.0; 
	D1[3] = 0.0; 
	D1[4] = 8.0/45.0; 
	D1[5] = 7.0*(7.0 + sQ)/360.0; 
	D1[6] = 7.0*(7.0 - sQ)/360.0; 
	D1[7] =-1.0/20.0; 
	D1[8] = 1.0/20.0; 
     
	D2[0] = 1.0/20.0; 
	D2[1] = 0.0; 
	D2[2] = 0.0; 
	D2[3] = 0.0; 
	D2[4] = 8.0/45.0; 
	D2[5] = 7.0*(7.0 + sQ)/360.0; 
	D2[6] = 7.0*(7.0 - sQ)/360.0; 
	D2[7] = 0.0, 
	D2[8] = 0.0;        
       
	D3[0] = 1.0/20.0; 
	D3[1] = 0.0; 
	D3[2] = 0.0; 
	D3[3] = 0.0; 
	D3[4] = 16.0/45.0; 
	D3[5] = 49.0/180.0; 
	D3[6] = 49.0/180.0; 
	D3[7] = 1.0/20.0; 
	D3[8] = 0.0; 
     
	D4[0] = 0.0; 
	D4[1] = 1.0/10.0; 
	D4[2] = 1.0/5.0; 
	D4[3] = 3.0/8.0; 
	D4[4] = 1.0/2.0; 
	D4[5] = (7.0 - sQ)/14.0; 
	D4[6] = (7.0 + sQ)/14.0; 
	D4[7] = 1.0; 
	D4[8] = 1.0; 

	for (int i=0; i<9; i++)
		for (int j=0; j<8; j++)
			D[i][j] = 0.0;

   	D[1][0] = 1.0/200.0; 
       
	D[2][0] = 1.0/150.0; 
	D[2][1] = 1.0/75.0; 
       
	D[3][0] = 171.0/8192.0; 
	D[3][1] = 45.0/4096.0; 
	D[3][2] = 315.0/8192.0; 
     
	D[4][0] = 5.0/288.0; 
	D[4][1] = 25.0/528.0; 
	D[4][2] = 25.0/672.0; 
	D[4][3] = 16.0/693.0; 
       
	D[5][0] = (1003.0 - 205.0*sQ)/12348.0; 
	D[5][1] =-25.0*(751.0 - 173.0*sQ)/90552.0; 
	D[5][2] = 25.0*(624.0 - 137.0*sQ)/43218.0; 
	D[5][3] =-128.0*(361.0 - 79.0*sQ)/237699.0; 
	D[5][4] = (3411.0 - 745.0*sQ)/24696.0; 
       
	D[6][0] = (793.0 + 187.0*sQ)/12348.0; 
	D[6][1] =-25.0*(331.0 + 113.0*sQ)/90552.0; 
	D[6][2] = 25.0*(1044.0 + 247.0*sQ)/43218.0; 
	D[6][3] =-128.0*(14885.0 + 3779.0*sQ)/9745659.0; 
	D[6][4] = (3327.0 + 797.0*sQ)/24696.0; 
	D[6][5] =-(581.0 + 127.0*sQ)/1722.0; 
     
	D[7][0] =-(157.0 - 3.0*sQ)/378.0; 
	D[7][1] = 25.0*(143.0 - 10.0*sQ)/2772.0; 
	D[7][2] =-25.0*(876.0 + 55.0*sQ)/3969.0; 
	D[7][3] = 1280.0*(913.0 + 18.0*sQ)/596673.0; 
	D[7][4] =-(1353.0 + 26.0*sQ)/2268.0; 
	D[7][5] = 7.0*(1777.0 + 377.0*sQ)/4428.0; 
	D[7][6] = 7.0*(5.0 - sQ)/36.0; 
       
	D[8][0] = 1.0/20.0; 
	D[8][4] = 8.0/45.0; 
	D[8][5] = 7.0*(7.0 + sQ)/360.0; 
	D[8][6] = 7.0*(7.0 - sQ)/360.0; 
}

int DormandPrince::Driver(BodyData *bodyData, Acceleration *acceleration, TimeLine *timeLine)
{
	bodyData->time	= timeLine->time;

	acceleration->evaluateGasDrag			= true;
	acceleration->evaluateTypeIMigration	= true;
	acceleration->evaluateTypeIIMigration	= true;

	// Calculate the acceleration in the initial point
	acceleration->Compute(timeLine->time, bodyData->y0, bodyData->accel);

	// NOTE: Kikapcsolom a GasDrag erők kiszámítását, gyorsítva ezzel az integrálást.
	// Készíteni összehasonlításokat, és értékelni az eredményeket, abbol a szempontbol, hogy így mennyire pontos az integralas.
	acceleration->evaluateGasDrag			= false;
	acceleration->evaluateTypeIMigration	= false;
	acceleration->evaluateTypeIIMigration	= false;

	int		iter = 0;
	double	errorMax = 0.0;
	do {
		iter++;
		bodyData->h		= timeLine->hNext;
		if (Step(bodyData, acceleration) == 1) {
			Error::_errMsg = "An error occurred during Prince-Dormand step!";
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			return 1;
		}
		errorMax = GetErrorMax(bodyData->nBodies.NOfVar(), bodyData->error);
		timeLine->hDid = bodyData->h;
		timeLine->hNext = errorMax < 1.0e-20 ? 2.0*bodyData->h : 0.9*bodyData->h*pow(epsilon / errorMax, 1.0/7.0);
	} while (errorMax > epsilon && iter <= maxIter);
	if (iter > maxIter) {
		Error::_errMsg = "An error occurred during Prince-Dormand driver: iteration number exceeded maxIter!";
		Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
		return 1;
	}
	timeLine->time += timeLine->hDid;
	bodyData->time = timeLine->time;

	// Update the phases of the system
	std::swap(bodyData->y0, bodyData->y);

	return 0;
}

int DormandPrince::Step(BodyData *bodyData, Acceleration *acceleration)
{
	// These arrays will contain the accelerations computed along the trajectory of the current step
	double	*fMatrix[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
	// Contains the approximation of the solution
	double	*yTemp = 0;
	int		nVar = bodyData->nBodies.NOfVar();

	fMatrix[0] = bodyData->accel;
	for (int i=1; i<9; i++) {
		fMatrix[i] = new double[nVar];
		HANDLE_NULL(fMatrix[i]);
		memset(fMatrix[i], 0, nVar*sizeof(double));
	}

	yTemp = new double[nVar];
	HANDLE_NULL(yTemp);
	memset(yTemp, 0, nVar*sizeof(double));

	double	h = bodyData->h;
	double	h2 = h*h;

	for (int k = 1; k < sizeHeightRKD; k++) {
		double	t = bodyData->time + D4[k] * h;
		for (int i = 0; i < bodyData->nBodies.total; i++) {
			int i0 = 6*i;
			for (int j = 0; j < 3; j++) {
				double	xi = bodyData->y0[i0 + j];
				double	vi = bodyData->y0[i0 + j + 3];
				// Compute the coordinates
				yTemp[i0 + j] = xi + D4[k]*h*vi;
				// Copy the velocities
				yTemp[i0 + j + 3] = vi;
				double sum_a = 0.0;
				for (int l = 0; l < k; l++) {
					sum_a += D[k][l]*fMatrix[l][i0 + j + 3];
				}
				// Compute the coordinates
				yTemp[i0 + j] += h2*sum_a;
				// Compute the velocities
				yTemp[i0 + j + 3] += h*sum_a;
			}
		}
		acceleration->Compute(t, yTemp, fMatrix[k]);
	}

	for (int i = 0; i < bodyData->nBodies.total; i++) {
		int i0 = 6*i;
		for (int j = 0; j<3; j++) {
			int n = i0 + j;
			// Compute the coordinates after the set
			bodyData->y[n]			 = bodyData->y0[n] + h*bodyData->y0[n + 3];
			bodyData->yBetterEst[n]  = bodyData->y[n];
			// Compute the velocities after the set
			bodyData->y[n + 3]		 = bodyData->y0[n + 3];
			for (int k = 0; k < sizeHeightRKD; k++)
			{
				double	a = fMatrix[k][n + 3];
				bodyData->y[n]			+= h2*D1[k]*a;
				bodyData->yBetterEst[n]	+= h2*D2[k]*a;
				bodyData->y[n + 3]		+=  h*D3[k]*a;
			}
			bodyData->error[n] = h2 * fabs(fMatrix[7][n + 3] - fMatrix[8][n + 3]) / 20.0;
			bodyData->error[n + 3] = 0.0;
		}
	}

	// i=0 changed to i=1
	for (int i=1; i<9; i++)
		delete[] fMatrix[i];
	delete[] yTemp;

	return 0;
}

double DormandPrince::GetErrorMax(const int n, const double *yerr)
{
	double errorMax = 0.0;

	for (int i = 0; i<n; i++) {
		double error = fabs(yerr[i]);
		if (error > errorMax)
			errorMax = error;
	}
	return errorMax;
}
