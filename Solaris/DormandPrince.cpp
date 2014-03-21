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

	b[0] = 1.0/20.0; 
	b[1] = 0.0; 
	b[2] = 0.0; 
	b[3] = 0.0; 
	b[4] = 8.0/45.0; 
	b[5] = 7.0*(7.0 + sQ)/360.0; 
	b[6] = 7.0*(7.0 - sQ)/360.0; 
	b[7] =-1.0/20.0; 
	b[8] = 1.0/20.0; 

	bh[0] = 1.0/20.0; 
	bh[1] = 0.0; 
	bh[2] = 0.0; 
	bh[3] = 0.0; 
	bh[4] = 8.0/45.0; 
	bh[5] = 7.0*(7.0 + sQ)/360.0; 
	bh[6] = 7.0*(7.0 - sQ)/360.0; 
	bh[7] = 0.0, 
	bh[8] = 0.0;        

	bdh[0] = 1.0/20.0; 
	bdh[1] = 0.0; 
	bdh[2] = 0.0; 
	bdh[3] = 0.0; 
	bdh[4] = 16.0/45.0; 
	bdh[5] = 49.0/180.0; 
	bdh[6] = 49.0/180.0; 
	bdh[7] = 1.0/20.0; 
	bdh[8] = 0.0; 

	c[0] = 0.0; 
	c[1] = 1.0/10.0; 
	c[2] = 1.0/5.0; 
	c[3] = 3.0/8.0; 
	c[4] = 1.0/2.0; 
	c[5] = (7.0 - sQ)/14.0; 
	c[6] = (7.0 + sQ)/14.0; 
	c[7] = 1.0; 
	c[8] = 1.0; 

	for (int i = 0; i < 9; i++) {
		for (int j = 0; j < 8; j++) {
			a[i][j] = 0.0;
		}
	}

	a[1][0] = 1.0/200.0; 
       
	a[2][0] = 1.0/150.0; 
	a[2][1] = 1.0/75.0; 
       
	a[3][0] = 171.0/8192.0; 
	a[3][1] = 45.0/4096.0; 
	a[3][2] = 315.0/8192.0; 
     
	a[4][0] = 5.0/288.0; 
	a[4][1] = 25.0/528.0; 
	a[4][2] = 25.0/672.0; 
	a[4][3] = 16.0/693.0; 
       
	a[5][0] = (1003.0 - 205.0*sQ)/12348.0; 
	a[5][1] =-25.0*(751.0 - 173.0*sQ)/90552.0; 
	a[5][2] = 25.0*(624.0 - 137.0*sQ)/43218.0; 
	a[5][3] =-128.0*(361.0 - 79.0*sQ)/237699.0; 
	a[5][4] = (3411.0 - 745.0*sQ)/24696.0; 
       
	a[6][0] = (793.0 + 187.0*sQ)/12348.0; 
	a[6][1] =-25.0*(331.0 + 113.0*sQ)/90552.0; 
	a[6][2] = 25.0*(1044.0 + 247.0*sQ)/43218.0; 
	a[6][3] =-128.0*(14885.0 + 3779.0*sQ)/9745659.0; 
	a[6][4] = (3327.0 + 797.0*sQ)/24696.0; 
	a[6][5] =-(581.0 + 127.0*sQ)/1722.0; 
     
	a[7][0] =-(157.0 - 3.0*sQ)/378.0; 
	a[7][1] = 25.0*(143.0 - 10.0*sQ)/2772.0; 
	a[7][2] =-25.0*(876.0 + 55.0*sQ)/3969.0; 
	a[7][3] = 1280.0*(913.0 + 18.0*sQ)/596673.0; 
	a[7][4] =-(1353.0 + 26.0*sQ)/2268.0; 
	a[7][5] = 7.0*(1777.0 + 377.0*sQ)/4428.0; 
	a[7][6] = 7.0*(5.0 - sQ)/36.0; 
       
	a[8][0] = 1.0/20.0; 
	a[8][4] = 8.0/45.0; 
	a[8][5] = 7.0*(7.0 + sQ)/360.0; 
	a[8][6] = 7.0*(7.0 - sQ)/360.0; 
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
	//acceleration->evaluateGasDrag			= false;
	acceleration->evaluateTypeIMigration	= false;
	acceleration->evaluateTypeIIMigration	= false;

	int		iter = 0;
	double	errorMax = 0.0;
	do {
		iter++;
		bodyData->h		= timeLine->hNext;
		//if (Step(bodyData, acceleration) == 1) {
		if (Step2(bodyData, acceleration) == 1) {
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
	double	*f[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
	// Contains the approximation of the solution
	double	*yTemp = 0;
	int		nVar = bodyData->nBodies.NOfVar();

	f[0] = bodyData->accel;
	for (int i = 1; i < 9; i++) {
		f[i] = new double[nVar];
		HANDLE_NULL(f[i]);
		memset(f[i], 0, nVar*sizeof(double));
	}

	yTemp = new double[nVar];
	HANDLE_NULL(yTemp);
	memset(yTemp, 0, nVar*sizeof(double));

	double	h = bodyData->h;
	double	h2 = h*h;

	for (int k = 1; k < sizeHeightRKD; k++) {
		double	ttemp = bodyData->time + c[k] * h;
		for (int i = 0; i < bodyData->nBodies.total; i++) {
			int i0 = 6*i;
			for (int j = 0; j < 3; j++) {
				int n = i0 + j;
				// Compute the coordinates
				yTemp[n]	= bodyData->y0[n] + c[k]*h*bodyData->y0[n+3];
				// Copy the velocities
				yTemp[n+3]	= bodyData->y0[n+3];
				double sum_a = 0.0;
				for (int l = 0; l < k; l++) {
					sum_a += a[k][l]*f[l][n+3];
				}
				// Compute the coordinates
				yTemp[n]	+= h2*sum_a;
				// Compute the velocities
				yTemp[n+3]	+= h*sum_a;
			}
		}
		acceleration->Compute(ttemp, yTemp, f[k]);
	}

	for (int i = 0; i < bodyData->nBodies.total; i++) {
		int i0 = 6*i;
		for (int j = 0; j<3; j++) {
			int n = i0 + j;
			// Compute the coordinates after the set
			bodyData->y[n]			 = bodyData->y0[n] + h*bodyData->y0[n+3];
			bodyData->yBetterEst[n]  = bodyData->y[n];
			// Compute the velocities after the set
			bodyData->y[n+3]		 = bodyData->y0[n+3];
			for (int k = 0; k < sizeHeightRKD; k++)
			{
				double	a = f[k][n+3];
				bodyData->y[n]			+= h2*b[k]*a;
				bodyData->yBetterEst[n]	+= h2*bh[k]*a;
				bodyData->y[n+3]		+=  h*bdh[k]*a;
			}
			bodyData->error[n] = h2 * fabs(f[7][n+3] - f[8][n+3]) / 20.0;
			bodyData->error[n+3] = 0.0;
		}
	}

	// i=0 changed to i=1
	for (int i = 1; i < 9; i++)
		delete[] f[i];
	delete[] yTemp;

	return 0;
}

int DormandPrince::Step2(BodyData *bodyData, Acceleration *acceleration)
{
	// This array will contain the accelerations computed along the trajectory of the current step
	double	*f[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
	// Contains the approximation of the solution
	double	*yTemp = 0;
	int		nVar = bodyData->nBodies.NOfVar();

	for (int i = 1; i < 9; i++) {
		f[i] = new double[nVar];
		HANDLE_NULL(f[i]);
		memset(f[i], 0, nVar*sizeof(double));
	}

	yTemp = new double[nVar];
	HANDLE_NULL(yTemp);
	memset(yTemp, 0, nVar*sizeof(double));

	int		n_total = bodyData->nBodies.total;
	double	h = bodyData->h;
	double	h2 = h*h;

	int k = 0;
	// f1 = f[0]
	double	ttemp = bodyData->time + c[k] * h;
	f[k] = bodyData->accel;

	// compute yTemp in order to compute f2
	k = 1;
	ttemp = bodyData->time + c[k] * h;
	for (int i = 0; i < n_total; i++) {
		int i0 = 6*i;
		for (int j = 0; j < 3; j++) {
			int n = i0 + j;
			double var	= a[k][0]*f[0][n+3];
			// Compute the new position
			yTemp[n]	= bodyData->y0[n] + c[k]*h*bodyData->y0[n+3] + h2*(var);
			// Compute the new velocity
			yTemp[n+3]	= bodyData->y0[n+3] + h*(var);
		}
	} // yTemp is computed
	// f2 = f[1]
	acceleration->Compute(ttemp, yTemp, f[k]);

	// compute yTemp in order to compute f3
	k = 2;
	ttemp = bodyData->time + c[k] * h;
	for (int i = 0; i < n_total; i++) {
		int i0 = 6*i;
		for (int j = 0; j < 3; j++) {
			int n = i0 + j;
			double var	= a[k][0]*f[0][n+3] + a[k][1]*f[1][n+3];
			// Compute the new position
			yTemp[n]	= bodyData->y0[n] + c[k]*h*bodyData->y0[n+3] + h2*(var);
			// Compute the new velocity
			yTemp[n+3]	= bodyData->y0[n+3] + h*(var);
		}
	} // yTemp is computed
	// f3 = f[2]
	acceleration->Compute(ttemp, yTemp, f[k]);

	// compute yTemp in order to compute f4
	k = 3;
	ttemp = bodyData->time + c[k] * h;
	for (int i = 0; i < n_total; i++) {
		int i0 = 6*i;
		for (int j = 0; j < 3; j++) {
			int n = i0 + j;
			double var	= a[k][0]*f[0][n+3] + a[k][1]*f[1][n+3] + a[k][2]*f[2][n+3];
			// Compute the new position
			yTemp[n]	= bodyData->y0[n] + c[k]*h*bodyData->y0[n+3] + h2*(var);
			// Compute the new velocity
			yTemp[n+3]	= bodyData->y0[n+3] + h*(var);
		}
	} // yTemp is computed
	// f4 = f[3]
	acceleration->Compute(ttemp, yTemp, f[k]);

	// compute yTemp in order to compute f5
	k = 4;
	ttemp = bodyData->time + c[k] * h;
	for (int i = 0; i < n_total; i++) {
		int i0 = 6*i;
		for (int j = 0; j < 3; j++) {
			int n = i0 + j;
			double var	= a[k][0]*f[0][n+3] + a[k][1]*f[1][n+3] + a[k][2]*f[2][n+3] + a[k][3]*f[3][n+3];
			// Compute the new position
			yTemp[n]	= bodyData->y0[n] + c[k]*h*bodyData->y0[n+3] + h2*(var);
			// Compute the new velocity
			yTemp[n+3]	= bodyData->y0[n+3] + h*(var);
		}
	} // yTemp is computed
	// f5 = f[4]
	acceleration->Compute(ttemp, yTemp, f[k]);

	// compute yTemp in order to compute f6
	k = 5;
	ttemp = bodyData->time + c[k] * h;
	for (int i = 0; i < n_total; i++) {
		int i0 = 6*i;
		for (int j = 0; j < 3; j++) {
			int n = i0 + j;
			double var	= a[k][0]*f[0][n+3] + a[k][1]*f[1][n+3] + a[k][2]*f[2][n+3] + a[k][3]*f[3][n+3] + a[k][4]*f[4][n+3];
			// Compute the new position
			yTemp[n]	= bodyData->y0[n] + c[k]*h*bodyData->y0[n+3] + h2*(var);
			// Compute the new velocity
			yTemp[n+3]	= bodyData->y0[n+3] + h*(var);
		}
	} // yTemp is computed
	// f6 = f[5]
	acceleration->Compute(ttemp, yTemp, f[k]);

	// compute yTemp in order to compute f7
	k = 6;
	ttemp = bodyData->time + c[k] * h;
	for (int i = 0; i < n_total; i++) {
		int i0 = 6*i;
		for (int j = 0; j < 3; j++) {
			int n = i0 + j;
			double var	= a[k][0]*f[0][n+3] + a[k][1]*f[1][n+3] + a[k][2]*f[2][n+3] + a[k][3]*f[3][n+3] + a[k][4]*f[4][n+3] +
						  a[k][5]*f[5][n+3];
			// Compute the new position
			yTemp[n]	= bodyData->y0[n] + c[k]*h*bodyData->y0[n+3] + h2*(var);
			// Compute the new velocity
			yTemp[n+3]	= bodyData->y0[n+3] + h*(var);
		}
	} // yTemp is computed
	// f7 = f[6]
	acceleration->Compute(ttemp, yTemp, f[k]);

	// compute yTemp in order to compute f8
	k = 7;
	ttemp = bodyData->time + c[k] * h;
	for (int i = 0; i < n_total; i++) {
		int i0 = 6*i;
		for (int j = 0; j < 3; j++) {
			int n = i0 + j;
			double var	= a[k][0]*f[0][n+3] + a[k][1]*f[1][n+3] + a[k][2]*f[2][n+3] + a[k][3]*f[3][n+3] + a[k][4]*f[4][n+3] +
						  a[k][5]*f[5][n+3] + a[k][6]*f[6][n+3];
			// Compute the new position
			yTemp[n]	= bodyData->y0[n] + c[k]*h*bodyData->y0[n+3] + h2*(var);
			// Compute the new velocity
			yTemp[n+3]	= bodyData->y0[n+3] + h*(var);
		}
	} // yTemp is computed
	// f8 = f[7]
	acceleration->Compute(ttemp, yTemp, f[k]);

	// compute yTemp in order to compute f9
	k = 8;
	ttemp = bodyData->time + c[k] * h;
	for (int i = 0; i < n_total; i++) {
		int i0 = 6*i;
		for (int j = 0; j < 3; j++) {
			int n = i0 + j;
			double var	= a[k][0]*f[0][n+3] + a[k][4]*f[4][n+3] + a[k][5]*f[5][n+3] + a[k][6]*f[6][n+3];
			// Compute the new position
			yTemp[n]	= bodyData->y0[n] + c[k]*h*bodyData->y0[n+3] + h2*(var);
			// Compute the new velocity
			yTemp[n+3]	= bodyData->y0[n+3] + h*(var);
		}
	} // yTemp is computed
	// f9 = f[8]
	acceleration->Compute(ttemp, yTemp, f[k]);

	//// yTemp to compute f3
	//for (int i = 0; i < nVar; i++) {
	//	yTemp[i] = bodyData->y0[i] + h*(a[2][0]*f[0][i] + a[2][1]*f[1][i]);
	//}

	//// f3 = f[2]
	//ttemp = bodyData->time + c[2] * h;
	//acceleration->Compute(ttemp, yTemp, f[2]);

	//for (int i = 0; i < nVar; i++) {
	//	yTemp[i] = bodyData->y0[i] + h*(a[3][0]*f[0][i] + a[3][1]*f[1][i] + a[3][2]*f[2][i]);
	//}

	//// f4 = f[3]
	//ttemp = bodyData->time + c[3] * h;
	//acceleration->Compute(ttemp, yTemp, f[3]);

	//for (int i = 0; i < nVar; i++) {
	//	yTemp[i] = bodyData->y0[i] + h*(a[4][0]*f[0][i] + a[4][1]*f[1][i] + a[4][2]*f[2][i] + a[4][3]*f[3][i]);
	//}

	//// f5 = f[4]
	//ttemp = bodyData->time + c[4] * h;
	//acceleration->Compute(ttemp, yTemp, f[4]);

	//for (int i = 0; i < nVar; i++) {
	//	yTemp[i] = bodyData->y0[i] + h*(a[5][0]*f[0][i] + a[5][1]*f[1][i] + a[5][2]*f[2][i] + a[5][3]*f[3][i] + a[5][4]*f[4][i]);
	//}

	//// f6 = f[5]
	//ttemp = bodyData->time + c[5] * h;
	//acceleration->Compute(ttemp, yTemp, f[5]);

	//for (int i = 0; i < nVar; i++) {
	//	yTemp[i] = bodyData->y0[i] + h*(a[6][0]*f[0][i] + a[6][1]*f[1][i] + a[6][2]*f[2][i] + a[6][3]*f[3][i] + a[6][4]*f[4][i] +
	//									a[6][5]*f[5][i]);
	//}

	//// f7 = f[6]
	//ttemp = bodyData->time + c[6] * h;
	//acceleration->Compute(ttemp, yTemp, f[6]);

	//for (int i = 0; i < nVar; i++) {
	//	yTemp[i] = bodyData->y0[i] + h*(a[7][0]*f[0][i] + a[7][1]*f[1][i] + a[7][2]*f[2][i] + a[7][3]*f[3][i] + a[7][4]*f[4][i] + 
	//									a[7][5]*f[5][i] + a[7][6]*f[6][i]);
	//}

	//// f8 = f[7]
	//ttemp = bodyData->time + c[7] * h;
	//acceleration->Compute(ttemp, yTemp, f[7]);

	//for (int i = 0; i < nVar; i++) {
	//	yTemp[i] = bodyData->y0[i] + h*(a[8][0]*f[0][i] + a[8][1]*f[1][i] + a[8][2]*f[2][i] + a[8][3]*f[3][i] + a[8][4]*f[4][i] + 
	//									a[8][5]*f[5][i] + a[8][6]*f[6][i]);
	//}

	//// f9 = f[8]
	//ttemp = bodyData->time + c[8] * h;
	//acceleration->Compute(ttemp, yTemp, f[8]);

	for (int i = 0; i < bodyData->nBodies.total; i++) {
		int i0 = 6 * i;
		for (int j = 0; j < 3; j++) {
			int n = i0 + j;
			bodyData->y[n] = bodyData->y0[n] + h*bodyData->y0[n+3] + h2*(b[0]*f[0][n+3] + b[4]*f[4][n+3] + b[5]*f[5][n+3] + 
																	     b[6]*f[6][n+3] + b[7]*f[7][n+3] + b[8]*f[8][n+3]);
			bodyData->error[n] = h2 * fabs(f[7][n+3] - f[8][n+3]) / 20.0;

			bodyData->y[n+3] = bodyData->y0[n+3] + h*(bdh[0]*f[0][n+3] + bdh[4]*f[4][n+3] + bdh[5]*f[5][n+3] + 
													  bdh[6]*f[6][n+3] + bdh[7]*f[7][n+3]);
			bodyData->error[n+3] = 0.0;
		}
	}

	// i=0 changed to i=1
	for (int i = 1; i < 9; i++)
		delete[] f[i];
	delete[] yTemp;

	return 0;
}

double DormandPrince::GetErrorMax(int n, const double *yerr)
{
	double errorMax = 0.0;

	for (int i = 0; i < n; i++) {
		double error = fabs(yerr[i]);
		if (error > errorMax)
			errorMax = error;
	}
	return errorMax;
}
