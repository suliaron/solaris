/*
 * RungeKuttaFehlberg.cpp
 *
 *  Created on: 2013.04.03.
 *      Author: Dr. Áron Süli (suli)
 */

#include <stdio.h>
#include <stdlib.h>
#include <cmath>

#include "RungeKuttaFehlberg78.h"
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

RungeKuttaFehlberg78::RungeKuttaFehlberg78()
{
	name		= "Classical eight-order Runge-Kutta-Nyström method with stepsize control.";
	reference	= "NASA Technical Reports R-381, by Erwin Fehlberg, 1972.";

	accuracy = -10.0;
	epsilon	 = pow(10, accuracy);

	D1_0 = 41.0/840.0, D1_1 = 0.0, D1_2 = 0.0, D1_3 = 0.0, D1_4 = 0.0, D1_5 = 34.0/105.0;
	D1_6 = 9.0/35.0, D1_7 = 9.0/35.0, D1_8 = 9.0/280.0, D1_9 = 9.0/280.0, D1_10 = 41.0/840.0;
	D1_11 = 0.0, D1_12 = 0.0;
	D_1_0 = 2.0/27.0, D_2_0 = 1.0/36.0, D_3_0 = 1.0/24.0, D_4_0 = 5.0/12.0;
	D_5_0 = 1.0/20.0, D_6_0 = -25.0/108.0, D_7_0 = 31.0/300.0, D_8_0 = 2.0;
	D_9_0 = -91.0/108.0, D_10_0 = 2383.0/4100.0, D_11_0 = 3.0/205.0, D_12_0 = -1777.0/4100.0;
	D_2_1 = 1.0/12.0;
	D_3_2 = 1.0/8.0, D_4_2 = -25.0/16.0;
	D_4_3 = 25.0/16.0, D_5_3 = 1.0/4.0, D_6_3 = 125.0/108.0, D_8_3 = -53.0/6.0, D_9_3 = 23.0/108.0, D_10_3 = -341.0/164.0, D_12_3 = -341.0/164.0;
	D_5_4 = 1.0/5.0, D_6_4 = -65.0/27.0, D_7_4 = 61.0/225.0, D_8_4 = 704.0/45.0, D_9_4 = -976.0/135.0, D_10_4 = 4496.0/1025.0, D_12_4 = 4496.0/1025.0;
	D_6_5 = 125.0/54.0, D_7_5 = -2.0/9.0, D_8_5 = -107.0/9.0, D_9_5 = 311.0/54.0, D_10_5 = -301.0/82.0, D_11_5 = -6.0/41.0, D_12_5 = -289.0/82.0;
	D_7_6 = 13.0/900.0, D_8_6 = 67.0/90.0, D_9_6 = -19.0/60.0, D_10_6 = 2133.0/4100.0, D_11_6 = -3.0/205.0, D_12_6 = 2193.0/4100.0;
	D_8_7 = 3.0, D_9_7 = 17.0/6.0, D_10_7 = 45.0/82.0, D_11_7 = -3.0/41.0, D_12_7 = 51.0/82.0;
	D_9_8 = -1.0/12.0, D_10_8 = 45.0/164.0, D_11_8 = 3.0/41.0, D_12_8 = 33.0/164.0;
	D_10_9 = 18.0/41.0, D_11_9 = 6.0/41.0, D_12_9 = 12.0/41.0;
	D_12_11 = 1.0;
}

// constants for the Runge-Kutta-Fehlberg7(8) integrator
#define SAFETY	 0.9
#define PGROW	-0.2
#define PSHRNK	-0.25
#define ERRCON	 1.89e-4
// For the scaling used to monitor accuracy
#define TINY	1.0e-30
int RungeKuttaFehlberg78::Driver(BodyData *bodyData, Acceleration *acceleration, TimeLine *timeLine)
{
#ifdef _DEBUG
//	fprintf(stderr, "File: %40s, Function: %40s, Line: %10d\n", __FILE__, __FUNCTION__, __LINE__);
#endif
	int	result = 0;
	int	nVar = bodyData->nBodies.NOfVar();

	bodyData->time	= timeLine->time;
	bodyData->h		= timeLine->hNext;

	acceleration->evaluateGasDrag			= true;
	acceleration->evaluateTypeIMigration	= true;
	acceleration->evaluateTypeIIMigration	= true;

	// Calculate the acceleration in the initial point
	acceleration->Compute(timeLine->time, bodyData->y0, bodyData->accel);

	// TODO: ezt a ciklust a végtelen ciklus belsejebe kell mozgatni, hiszen ha a lépésköz változik akkor annak
	// megfelelően ennek is változnia kell. Ezt átgondolni.
	// Most megnézem. Nem javított a helyzeten, ezért visszateszem a cikoson kívulre.
	for (int i=0; i<nVar; i++) {
		bodyData->yscale[i] = fabs(bodyData->y0[i]) + fabs(bodyData->h * bodyData->accel[i]) + TINY;
	}

	// NOTE: Kikapcsolom a GasDrag erők kiszámítását, gyorsítva ezzel az integrálást.
	// Készíteni összehasonlításokat, és értékelni az eredményeket, abbol a szempontbol, hogy így mennyire pontos az integralas.
	//acceleration->evaluateGasDrag			= false;
	acceleration->evaluateTypeIMigration	= false;
	acceleration->evaluateTypeIIMigration	= false;

	double	errorMax = 0.0;
	for ( ; ; ) {

		//if (Step(bodyData, acceleration, bodyData->accel, bodyData->time, bodyData->h, bodyData->y, bodyData->error) == 1) {
		if (Step(bodyData, acceleration) == 1) {
			Error::_errMsg = "An error occurred during Runge-Kutta-Fehlberg7(8) step!";
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			result = 1;
			break;
		}
		errorMax = GetErrorMax(nVar, bodyData->error, bodyData->yscale);
		if (errorMax < 1.0) {
			timeLine->hDid = bodyData->h;
			result = 0;
			break; 		/* step succeeded, exit from the infinite loop */
		}
		double hTemp = SAFETY * bodyData->h * pow(errorMax, PSHRNK);
		bodyData->h = fabs(hTemp) > fabs(0.1*bodyData->h) ? hTemp : 0.1*bodyData->h;

		double tNew = timeLine->time + bodyData->h;
		if (tNew == timeLine->time) {
			Error::_errMsg = "Stepsize-underflow occurred during Runge-Kutta-Fehlberg7(8) step!";
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			result = 1;
			break;
		}
	}
	if (result == 0) {
		// Update time
		timeLine->time += timeLine->hDid;
		bodyData->time = timeLine->time;
		// Calculate the next stepsize
		timeLine->hNext = errorMax > ERRCON ? (SAFETY * bodyData->h * pow(errorMax, PGROW)) : (5.0 * bodyData->h);
		bodyData->h = timeLine->hNext;
		// Update the phases of the system
		std::swap(bodyData->y0, bodyData->y);
	}

	//delete[] accel;
	//delete[] yscale;
	//delete[] yerr;

	return result;
}
#undef SAFETY
#undef PGROW
#undef PSHRNK
#undef ERRCON
#undef TINY

int RungeKuttaFehlberg78::Step(BodyData *bodyData, Acceleration *acceleration)
{
	// These arrays will contain the accelerations computed along the trajectory of the current step
	double	*fk[13] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
	// Contains the approximation of the solution
	double	*yTemp = 0;
	int		nVar = bodyData->nBodies.NOfVar();

	// i=0 changed to i=1
	for (int i=1; i<13; i++) {
		fk[i] = new double[nVar];
		HANDLE_NULL(fk[i]);
	}
	yTemp = new double[nVar];
	HANDLE_NULL(yTemp);

	// Copy the initial acceleration into fk[0]
	// NOTE: this copy can be avoided if a is used instead of fk[0], than we do not need to allocate/free fk[0]
	//memcpy(fk[0], accel, nVar*sizeof(double));
	fk[0] = bodyData->accel;
	double	h = bodyData->h;
	double	t = bodyData->time;
//1. substep
	for (int i=0; i<nVar; i++)
		yTemp[i] = bodyData->y0[i] + h*(D_1_0*fk[0][i]);

	acceleration->Compute(t, yTemp, fk[1]);
//2. substep
	for (int i=0; i<nVar; i++)
		yTemp[i] = bodyData->y0[i] + h*(D_2_0*fk[0][i] + D_2_1*fk[1][i]);

	acceleration->Compute(t, yTemp, fk[2]);
//3. substep
	for (int i=0; i<nVar; i++)
		yTemp[i] = bodyData->y0[i] + h*(D_3_0*fk[0][i] + D_3_2*fk[2][i]);

	acceleration->Compute(t, yTemp, fk[3]);
//4. substep
	for (int i=0; i<nVar; i++)
		yTemp[i] = bodyData->y0[i] + h*(D_4_0*fk[0][i] + D_4_2*fk[2][i] + D_4_3*fk[3][i]);

	acceleration->Compute(t, yTemp, fk[4]);
//5. substep
	for (int i=0; i<nVar; i++)
		yTemp[i] = bodyData->y0[i] + h*(D_5_0*fk[0][i] + D_5_3*fk[3][i] + D_5_4*fk[4][i]);

	acceleration->Compute(t, yTemp, fk[5]);
//6. substep
	for (int i=0; i<nVar; i++)
		yTemp[i] = bodyData->y0[i] + h*(D_6_0*fk[0][i] + D_6_3*fk[3][i] + D_6_4*fk[4][i] + D_6_5*fk[5][i]);

	acceleration->Compute(t, yTemp, fk[6]);
//7. substep
	for (int i=0; i<nVar; i++)
		yTemp[i] = bodyData->y0[i] + h*(D_7_0*fk[0][i] + D_7_4*fk[4][i] + D_7_5*fk[5][i] + D_7_6*fk[6][i]);

	acceleration->Compute(t, yTemp, fk[7]);
//8. substep
	for (int i=0; i<nVar; i++)
		yTemp[i] = bodyData->y0[i] + h*(D_8_0*fk[0][i] + D_8_3*fk[3][i] + D_8_4*fk[4][i] +
						  D_8_5*fk[5][i] + D_8_6*fk[6][i] + D_8_7*fk[7][i]);

	acceleration->Compute(t, yTemp, fk[8]);
//9. substep
	for (int i=0; i<nVar; i++)
		yTemp[i] = bodyData->y0[i] + h*(D_9_0*fk[0][i] + D_9_3*fk[3][i] + D_9_4*fk[4][i] +
						  D_9_5*fk[5][i] + D_9_6*fk[6][i] + D_9_7*fk[7][i] + D_9_8*fk[8][i]);

	acceleration->Compute(t, yTemp, fk[9]);
//10. substep
	for (int i=0; i<nVar; i++)
		yTemp[i] = bodyData->y0[i] + h*(D_10_0*fk[0][i] + D_10_3*fk[3][i] + D_10_4*fk[4][i] + D_10_5*fk[5][i] +
						  D_10_6*fk[6][i] + D_10_7*fk[7][i] + D_10_8*fk[8][i] + D_10_9*fk[9][i]);

	acceleration->Compute(t, yTemp, fk[10]);
//11. substep
	for (int i=0; i<nVar; i++)
		yTemp[i] = bodyData->y0[i] + h*(D_11_0*fk[0][i] + D_11_5*fk[5][i] + D_11_6*fk[6][i] +
						  D_11_7*fk[7][i] + D_11_8*fk[8][i] + D_11_9*fk[9][i]);

	acceleration->Compute(t, yTemp, fk[11]);
//12. substep
	for (int i=0; i<nVar; i++)
		yTemp[i] = bodyData->y0[i] + h*(D_12_0*fk[0][i] + D_12_3*fk[3][i] + D_12_4*fk[4][i] + D_12_5*fk[5][i] +
						  D_12_6*fk[6][i] + D_12_7*fk[7][i] + D_12_8*fk[8][i] + D_12_9*fk[9][i] +
						  D_12_11*fk[11][i]);

	acceleration->Compute(t, yTemp, fk[12]);
	// The result of the step
	for (int i=0; i<nVar; i++)
		bodyData->y[i] = bodyData->y0[i] + h*(D1_0*fk[0][i] + D1_5*fk[5][i] + D1_6*(fk[6][i] + fk[7][i]) +
							D1_8*(fk[8][i] + fk[9][i]) + D1_10*fk[10][i]);

	// Error estimation
	for (int i=0; i<nVar; i++)
		bodyData->error[i] = h*fabs(fk[0][i] + fk[10][i]-fk[11][i]-fk[12][i])*41.0/840.0;

	// i=0 changed to i=1
	for (int i=1; i<13; i++)
		delete[] fk[i];
	delete[] yTemp;

	return 0;
}

double RungeKuttaFehlberg78::GetErrorMax(const int n, const double *yerr, const double *yscale)
{
	double errorMax = 0.0;

	for (int i=0; i<n; i++) {
		double err = fabs(yerr[i] / yscale[i]);
		if (err > errorMax)
			errorMax = err;
	}
	return errorMax / epsilon;
}
