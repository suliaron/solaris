#include <cstdio>
#include <iostream>
#include <string.h>

#include "Acceleration.h"
#include "BodyData.h"
#include "Calculate.h"
#include "Ephemeris.h"
#include "Error.h"
#include "Event.h"
#include "EventCondition.h"
#include "IntegratorType.h"
#include "Nebula.h"
#include "Tools.h"
#include "TwoBodyAffair.h"

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

#define SQR(a)		((a)*(a))
#define CUBE(a)		((a)*(a)*(a))
#define FORTH(a)	((a)*(a)*(a)*(a))
#define FIFTH(a)	((a)*(a)*(a)*(a)*(a))

Acceleration::Acceleration(IntegratorType iType, bool baryCentric, BodyData *bD, Nebula *n)
{
	_integratorType			= iType;
	_baryCentric			= baryCentric;
	evaluateGasDrag			= true;
	evaluateTypeIMigration	= true;
	evaluateTypeIIMigration = true;

	bodyData	= bD;
	nebula		= n;

	rm3					= 0;
	accelGasDrag		= 0;
	accelMigrationTypeI	= 0;
	accelMigrationTypeII= 0;
}

Acceleration::~Acceleration()
{
	delete[] rm3;
	delete[] accelGasDrag;
	delete[] accelMigrationTypeI;
	delete[] accelMigrationTypeII;
}

int	Acceleration::Compute(double t, double *y, double *totalAccel)
{
	int	result = 0;
	// The inverse of the cube of the distance from the central body
	// nBodies contains the number of the different bodies
	if (rm3 == 0 /*or reallocation must performed*/ ) {
		rm3	= new double[bodyData->nBodies.total];
		HANDLE_NULL(rm3);
		memset(rm3, 0, bodyData->nBodies.total*sizeof(double));
	}

	if (_baryCentric) {
		result = ComputeBaryCentric(t, y, totalAccel);
		HANDLE_RESULT(result);
	}
	else {
		result = ComputeAstroCentric(t, y, totalAccel);
		HANDLE_RESULT(result);
	}

	return 0;
}

int Acceleration::ComputeBaryCentric(double t, double *y, double *totalAccel)
{
	int result = GravityBC(t, y, totalAccel);
	HANDLE_RESULT(result);

	if (nebula != 0) {
		if (accelGasDrag == 0 && bodyData->nBodies.NOfPlAndSpl() > 0) {
			accelGasDrag = new double[3*bodyData->nBodies.NOfPlAndSpl()];
			HANDLE_NULL(accelGasDrag);
			memset(accelGasDrag, 0, 3*bodyData->nBodies.NOfPlAndSpl()*sizeof(double));
		}
		if (evaluateGasDrag && accelGasDrag != 0) {
			result = GasDragBC(t, y, accelGasDrag);
			HANDLE_RESULT(result);
		}
		int lower = bodyData->nBodies.NOfMassive();
		int upper = lower + bodyData->nBodies.NOfPlAndSpl();
		for (int i=lower; i<upper; i++) {
			int i0 = 6*i;
			int j0 = 3*(i-lower);
			totalAccel[i0 + 3] += accelGasDrag[j0 + 0];
			totalAccel[i0 + 4] += accelGasDrag[j0 + 1];
			totalAccel[i0 + 5] += accelGasDrag[j0 + 2];
		}

		if (accelMigrationTypeI == 0 && bodyData->nBodies.protoPlanet) {
			accelMigrationTypeI = new double[3*bodyData->nBodies.protoPlanet];
			HANDLE_NULL(accelMigrationTypeI);
			memset(accelMigrationTypeI, 0, 3*bodyData->nBodies.protoPlanet*sizeof(double));
		}
		if (evaluateTypeIMigration && accelMigrationTypeI != 0) {
			result = MigrationTypeIBC(t, y, accelMigrationTypeI);
			HANDLE_RESULT(result);
		}
		lower = bodyData->nBodies.centralBody + bodyData->nBodies.giantPlanet;
		upper = bodyData->nBodies.NOfMassive();
		for (int i=lower; i<upper; i++) {
			if (bodyData->migType[i] == No) {
				continue;
			}
			int	i0 = 6*i;
			int j0 = 3*(i-lower);
			totalAccel[i0 + 3] += accelMigrationTypeI[j0 + 0];
			totalAccel[i0 + 4] += accelMigrationTypeI[j0 + 1];
			totalAccel[i0 + 5] += accelMigrationTypeI[j0 + 2];
		}

		if (accelMigrationTypeII == 0 && bodyData->nBodies.giantPlanet > 0) {
			accelMigrationTypeII = new double[3*bodyData->nBodies.giantPlanet];
			HANDLE_NULL(accelMigrationTypeII);
			memset(accelMigrationTypeII, 0, 3*bodyData->nBodies.giantPlanet*sizeof(double));
		}
		if (evaluateTypeIIMigration && accelMigrationTypeII != 0) {
			result = MigrationTypeIIBC(t, y, accelMigrationTypeII);
			HANDLE_RESULT(result);
		}
		lower = bodyData->nBodies.centralBody;
		upper = bodyData->nBodies.centralBody + bodyData->nBodies.giantPlanet;
		for (int i=lower; i<upper; i++) {
			if (bodyData->migType[i] == No) {
				continue;
			}
			int	i0 = 6*i;
			int j0 = 3*(i-lower);
			totalAccel[i0 + 3] += accelMigrationTypeII[j0 + 0];
			totalAccel[i0 + 4] += accelMigrationTypeII[j0 + 1];
			totalAccel[i0 + 5] += accelMigrationTypeII[j0 + 2];
		}
		// TODO: check this threshold and delete statements
		if (nebula->gasComponent.ReductionFactor(t) < 1.0e-5) {
			evaluateGasDrag			= false;
			evaluateTypeIMigration	= false;
			evaluateTypeIIMigration = false;
			delete[] accelGasDrag;
			delete[] accelMigrationTypeI;
			delete[] accelMigrationTypeII;
			delete nebula;
		}
	}

	return 0;
}

int Acceleration::ComputeAstroCentric(double t, double *y, double *totalAccel)
{
	int result = GravityAC(t, y, totalAccel);
	HANDLE_RESULT(result);

#ifdef TEST_DUSTPARTICLE
	double a_g = sqrt(SQR(totalAccel[9]) + SQR(totalAccel[10]) + SQR(totalAccel[11]));
	printf("a_g: %30.15lf AU/d^2\n", a_g);
#endif

	if (nebula != 0) {
		if (accelGasDrag == 0 && bodyData->nBodies.NOfPlAndSpl() > 0) {
			accelGasDrag = new double[3*bodyData->nBodies.NOfPlAndSpl()];
			HANDLE_NULL(accelGasDrag);
			memset(accelGasDrag, 0, 3*bodyData->nBodies.NOfPlAndSpl()*sizeof(double));
		}
		if (evaluateGasDrag && accelGasDrag != 0) {
			result = GasDragAC(t, y, accelGasDrag);
			HANDLE_RESULT(result);
#ifdef TEST_DUSTPARTICLE
	double a_d = sqrt(SQR(accelGasDrag[0]) + SQR(accelGasDrag[1]) + SQR(accelGasDrag[2]));
	printf("a_d: %30.15lf AU/d^2\n", a_d);
#endif
        }
		int lower = bodyData->nBodies.NOfMassive();
		int upper = lower + bodyData->nBodies.NOfPlAndSpl();
		for (int i=lower; i<upper; i++) {
			int i0 = 6*i;
			int j0 = 3*(i-lower);
			totalAccel[i0 + 3] += accelGasDrag[j0 + 0];
			totalAccel[i0 + 4] += accelGasDrag[j0 + 1];
			totalAccel[i0 + 5] += accelGasDrag[j0 + 2];
		}

		if (accelMigrationTypeI == 0 && bodyData->nBodies.protoPlanet) {
			accelMigrationTypeI = new double[3*bodyData->nBodies.protoPlanet];
			HANDLE_NULL(accelMigrationTypeI);
			memset(accelMigrationTypeI, 0, 3*bodyData->nBodies.protoPlanet*sizeof(double));
		}
		if (evaluateTypeIMigration && accelMigrationTypeI != 0) {
			result = MigrationTypeIAC(t, y, accelMigrationTypeI);
			HANDLE_RESULT(result);
		}
		lower = bodyData->nBodies.centralBody + bodyData->nBodies.giantPlanet;
		upper = bodyData->nBodies.NOfMassive();
		for (int i=lower; i<upper; i++) {
			if (bodyData->migType[i] == No) {
				continue;
			}
			int	i0 = 6*i;
			int j0 = 3*(i-lower);
			totalAccel[i0 + 3] += accelMigrationTypeI[j0 + 0];
			totalAccel[i0 + 4] += accelMigrationTypeI[j0 + 1];
			totalAccel[i0 + 5] += accelMigrationTypeI[j0 + 2];
		}

		if (accelMigrationTypeII == 0 && bodyData->nBodies.giantPlanet > 0) {
			accelMigrationTypeII = new double[3*bodyData->nBodies.giantPlanet];
			HANDLE_NULL(accelMigrationTypeII);
			memset(accelMigrationTypeII, 0, 3*bodyData->nBodies.giantPlanet*sizeof(double));
		}
		if (evaluateTypeIIMigration && accelMigrationTypeII != 0) {
			result = MigrationTypeIIAC(t, y, accelMigrationTypeII);
			HANDLE_RESULT(result);
		}
		lower = bodyData->nBodies.centralBody;
		upper = bodyData->nBodies.centralBody + bodyData->nBodies.giantPlanet;
		for (int i=lower; i<upper; i++) {
			if (bodyData->migType[i] == No) {
				continue;
			}
			int	i0 = 6*i;
			int j0 = 3*(i-lower);
			totalAccel[i0 + 3] += accelMigrationTypeII[j0 + 0];
			totalAccel[i0 + 4] += accelMigrationTypeII[j0 + 1];
			totalAccel[i0 + 5] += accelMigrationTypeII[j0 + 2];
		}
	}

	return 0;
}

int	Acceleration::GravityAC(double t, double *y, double *accel)
{
	double	rij = 0, rij2 = 0, rijm3 = 0; 
	double	xij = 0, yij = 0, zij = 0; 

	bodyData->indexOfNN[0] = -1;
	bodyData->distanceOfNN[0] = 0.0;
	// Note: i=1, since y(0,1,2,3,4,5) = 0, central body
	for (int i=1; i<bodyData->nBodies.total; i++) {
		// y: (x, y, z), (vx, vy, vz), etc.
		int i0 = 6*i; 
		double r2 = SQR(y[i0 + 0]) + SQR(y[i0 + 1]) + SQR(y[i0 + 2]); 
		double r = sqrt(r2); 
		rm3[i]= 1.0 / (r2 * r);
		bodyData->indexOfNN[i] = -1;
		bodyData->distanceOfNN[i] = 0.0;
	}

	accel[0] = accel[1] = accel[2] = accel[3] = accel[4] = accel[5] = 0.0;
	int NOfMassive = bodyData->nBodies.NOfMassive();
	for (register int i=1; i<bodyData->nBodies.total; i++) {
		double rMin = 1.0e10;
		double ax = 0.0, ay = 0.0, az = 0.0; 

		double mu = Constants::Gauss2*(bodyData->mass[0] + bodyData->mass[i]); 
		register int i0 = 6*i; 
		accel[i0 + 0] = y[i0 + 3]; 
		accel[i0 + 1] = y[i0 + 4]; 
		accel[i0 + 2] = y[i0 + 5]; 
       /* First the Kepler term is calculated and stored. This is done, because the Kepler term is much larger than 
		* the terms from the body-body interaction, and in this way the rounding error is
        * reduced. This hypothesis should be checked!
        */
		accel[i0 + 3] = -mu*rm3[i]*y[i0 + 0]; 
		accel[i0 + 4] = -mu*rm3[i]*y[i0 + 1]; 
		accel[i0 + 5] = -mu*rm3[i]*y[i0 + 2]; 

		if (bodyData->type[i] <= ProtoPlanet) {
			NOfMassive = bodyData->nBodies.NOfMassive() + bodyData->nBodies.superPlanetsimal;
		} else {
			NOfMassive = bodyData->nBodies.NOfMassive();
		}
#ifdef _DEBUG
		//if (print) fprintf(stderr, "NOfMassive: %d\t bodies: ", NOfMassive);
		//if (print) fprintf(stderr, "\n[%d t: %s] ", i, BodyTypeStr[type[i]]);
#endif
		for (register int j=1; j<NOfMassive; j++) {
			if (j == i)  // the bodies do not interact gravitationally with themselves
				continue;
#ifdef _DEBUG
//		//if (print) fprintf(stderr, "(%d t: %s) ", j, BodyTypeStr[type[j]]);
#endif
			register int j0 = 6*j; 
			xij = y[j0 + 0] - y[i0 + 0]; 
			yij = y[j0 + 1] - y[i0 + 1]; 
			zij = y[j0 + 2] - y[i0 + 2]; 
			rij2 = SQR(xij) + SQR(yij) + SQR(zij); 
			rij = sqrt(rij2); 
			rijm3 = 1.0/(rij2*rij); 

			if (rij < rMin) {
				rMin = rij;
				bodyData->indexOfNN[i] = j;
				bodyData->distanceOfNN[i] = rij;
			}
			double Gmj = Constants::Gauss2*bodyData->mass[j];
			ax += Gmj*( xij*rijm3 - y[j0 + 0]*rm3[j] ); 
			ay += Gmj*( yij*rijm3 - y[j0 + 1]*rm3[j] ); 
			az += Gmj*( zij*rijm3 - y[j0 + 2]*rm3[j] ); 
		}
	   /*
        * TODO: The acceleration from the body-body interaction (ax, ay, az) is added to the
        * Kepler-term. I hope, that the rounding error has been reduced. This
        * should be checked!
        */
		accel[i0 + 3] += ax;
		accel[i0 + 4] += ay; 
		accel[i0 + 5] += az; 
	}

	return 0;
}

int	Acceleration::GasDragAC(double t, double *y, double *accel)
{
	static bool _epstein = false;
    static bool _stokes = false;
    static bool _transition = false;

	double factor = nebula->gasComponent.ReductionFactor(t);

	// Gas drag is experienced only by planetesimals and super-planetesimlas.
	int lower = bodyData->nBodies.NOfMassive();
	int upper = lower + bodyData->nBodies.planetsimal + bodyData->nBodies.superPlanetsimal;
	for (int i=lower; i<upper; i++) {
		int i0 = 6*i;
		int j0 = 3*(i-lower);
		double r = sqrt(SQR(y[i0 + 0]) + SQR(y[i0 + 1]) + SQR(y[i0 + 2]));

        double C = 0.0;
        {
			Vector rVec(y[i0 + 0], y[i0 + 1], y[i0 + 2]);
			Vector vVec(y[i0 + 3], y[i0 + 4], y[i0 + 5]);
	        // Calculate the gas velocity at rVec
			Vector vGas = nebula->gasComponent.gas_velocity(Constants::Gauss2*bodyData->mass[0], &rVec);
			// Vector vGas = nebula->gasComponent.GasVelocity(Constants::Gauss2*bodyData->mass[0], r, atan2(y[i0 + 1], y[i0 + 0]));
		    // The relative velocity
			Vector u(vVec.x - vGas.x, vVec.y - vGas.y, vVec.z - vGas.z);
			double rhoGas = factor * nebula->gasComponent.gas_density_at(&rVec);
			//double rhoGas = factor * nebula->gasComponent.GasDensityAt(r, y[i0 + 2]);

            double lambda = this->nebula->gasComponent.meanFreePath.Evaluate(r);
            // Epstein regime
            if (     bodyData->radius[i] <= 0.1 * lambda) {
                if (!_epstein) {
					std::cout << "t: " << t*Constants::DayToYear << "[yr] " << "Epstein regime\n";
                    _epstein = true;
                }
	            // Calculate the mean thermal velocity of the gas molecules
                double vth = nebula->gasComponent.MeanThermalSpeed_CMU(bodyData->mass[0], r);
                C = bodyData->gammaEpstein[i] * vth  * rhoGas;
            }
            // Stokes regime
            else if (bodyData->radius[i] >= 10.0 * lambda) {
                if (!_stokes) {
                    std::cout << "t: " << t*Constants::DayToYear << "[yr] " << "Stokes regime\n";
                    _stokes = true;
                }
                // magnitude of the relative velocity
		        double uLength = sqrt(u.x*u.x + u.y*u.y + u.z*u.z); 
		        C = bodyData->gammaStokes[i] * uLength * rhoGas;
            }
            // Transition regime
            else {
                if (!_transition) {
                    std::cout << "t: " << t*Constants::DayToYear << "[yr] " << "Transition regime\n";
                    _transition = true;
                }
                double lambda1 = 0.1 * lambda;
                double lambda2 = 10.0 * lambda;
                double gammaE = 1.0/(bodyData->density[i] * lambda1);
                double gammaS = 3.0/8.0 * bodyData->cD[i] / (bodyData->density[i] * lambda2);

                double vth = nebula->gasComponent.MeanThermalSpeed_CMU(bodyData->mass[0], r);
                double K = gammaS * u.Length() / (gammaE * vth);
                double eta = lambda2 / lambda1;
                double kappa = log10(K) / log10(eta);
                double gamma = gammaE * vth * pow(lambda1, -kappa);
                C = gamma * pow(bodyData->radius[i], kappa) * rhoGas;
            }
		    accel[j0 + 0] = -C * u.x;
		    accel[j0 + 1] = -C * u.y;
		    accel[j0 + 2] = -C * u.z;
        }
#if 0
		Vector rVec(y[i0 + 0], y[i0 + 1], y[i0 + 2]);
		Vector vGas = nebula->gasComponent.gas_velocity(Constants::Gauss2*bodyData->mass[0], &rVec);
		// The relative velocity
		Vector u(y[i0 + 3] - vGas.x, y[i0 + 4] - vGas.y, y[i0 + 5] - vGas.z);
		double uLength = sqrt(u.x*u.x + u.y*u.y + u.z*u.z); 
	
		double rhoGas = factor * nebula->gasComponent.gas_density_at(&rVec);
		C = bodyData->gammaStokes[i] * rhoGas * uLength;

		// TODO: implement the Epstein acceleration form!
		accel[j0 + 0] = -C * u.x;
		accel[j0 + 1] = -C * u.y;
		accel[j0 + 2] = -C * u.z;

		//printf("gasAccel[%d]: (%le, %le, %le)\n", i, accel[j0 + 0], accel[j0 + 1], accel[j0 + 2]);
#endif
	}

	return 0;
}

int Acceleration::MigrationTypeIAC(double t, double *y, double *accel)
{
	double factor = this->nebula->gasComponent.ReductionFactor(t);

	int lower = bodyData->nBodies.centralBody + bodyData->nBodies.giantPlanet;
	int upper = bodyData->nBodies.NOfMassive();
	for (int i=lower; i<upper;  i++) {
		if (bodyData->migType[i] != TypeI) {
			continue;
		}

		int	i0 = 6*i;
		int j0 = 3*(i-lower);
		double r2 = SQR(y[i0 + 0]) + SQR(y[i0 + 1]) + SQR(y[i0 + 2]);
		double r = sqrt(r2);
		if (r <= bodyData->migStopAt[i]) {
			accel[j0 + 0] = accel[j0 + 1] = accel[j0 + 2] = 0.0;
			bodyData->migType[i] = No;
			continue;
		}

		double m = this->bodyData->mass[i];
		double mc= this->bodyData->mass[0];

		double a = 0.0, e = 0.0;
		Phase phase;
		Tools::ToPhase(&y[i0], &phase);
		double mu = Constants::Gauss2*(mc + m);
		Ephemeris::CalculateOrbitalElement(mu, &phase, &a, &e);
		// Orbital frequency: (note, that this differs from the formula of Fogg & Nelson 2005)
		double O = Constants::Gauss * sqrt((mc + m)/CUBE(a));
		double C = SQR(mc)/(m*this->nebula->gasComponent.MidplaneDensity(r)*a*a);
		// Aspect ratio:
		double h = this->nebula->gasComponent.scaleHeight.Evaluate(r);
		double ar = h/r;
		double er = e*r;

		/*
		 *  When e > 1.1 h/r, inward migration halts as $t_{\rm m}$ becomes negative and only
		 *  resumes when eccentricity is damped to lower values. We note that under certain
		 *  circumstances, such as there being a surface density jump, or an optically thick disk,
		 *  or MHD turbulence, type I migration may be substantially modified or reversed
		 *  (Papaloizou & Nelson 2005; Paardekooper & Mellema 2006; Nelson 2005; Masset et al. 2006).
		 */
		double tm = 0.0;
		if (e < 1.1*h/r) {
			tm = TypeIMigrationTime(C, O, ar, er, h);
			tm = 1.0/tm;
		}
		double te = TypeIEccentricityDampingTime(C, O, ar, er, h);
		double ti = te;
		double vr = y[i0 + 0]*y[i0 + 3] + y[i0 + 1]*y[i0 + 4] + y[i0 + 2]*y[i0 + 5];
		te = 2.0*vr/(r2*te);
		ti = 2.0/ti;

		accel[j0 + 0] = -factor*(tm*y[i0 + 3] + te*y[i0 + 0]);
		accel[j0 + 1] = -factor*(tm*y[i0 + 4] + te*y[i0 + 1]);
		accel[j0 + 2] = -factor*(tm*y[i0 + 5] + te*y[i0 + 2] + ti*y[i0 + 5]);
	}

	return 0;
}

int Acceleration::MigrationTypeIIAC(double t, double *y, double *accel)
{
	double factor = this->nebula->gasComponent.ReductionFactor(t);

	int lower = bodyData->nBodies.centralBody;
	int upper = bodyData->nBodies.centralBody + bodyData->nBodies.giantPlanet;
	for (int i=lower; i<upper; i++) {
		if (bodyData->migType[i] != TypeII) {
			continue;
		}

		int	i0 = 6*i;
		int j0 = 3*(i-lower);
		double r2 = SQR(y[i0 + 0]) + SQR(y[i0 + 1]) + SQR(y[i0 + 2]);
		double r = sqrt(r2);
		if (r <= bodyData->migStopAt[i]) {
			bodyData->migType[i] = No;
			accel[j0 + 0] = accel[j0 + 1] = accel[j0 + 2] = 0.0;
			continue;
		}

		double m = this->bodyData->mass[i];
		double mc= this->bodyData->mass[0];
		double a = 0.0, e = 0.0;
		Phase phase;
		Tools::ToPhase(&y[i0], &phase);
		double mu = Constants::Gauss2*(mc + m);
		Ephemeris::CalculateOrbitalElement(mu, &phase, &a, &e);
		// Orbital frequency: (note, that this differs from the formula of Fogg & Nelson 2005)
		double O = Constants::Gauss * sqrt((mc + m)/CUBE(a));

		double c0 = 1.0/TauNu(r, O);
		double vr = y[i0 + 3]*y[i0 + 0] + y[i0 + 4]*y[i0 + 1] + y[i0 + 5]*y[i0 + 2];
		double c1 = vr/r2;

		accel[j0 + 0] = -factor*(c0 * ( 0.5*y[i0 + 3] + 50 * (c1*y[i0 + 0]) ));
		accel[j0 + 1] = -factor*(c0 * ( 0.5*y[i0 + 4] + 50 * (c1*y[i0 + 1]) ));
		accel[j0 + 2] = -factor*(c0 * ( 0.5*y[i0 + 5] + 50 * (c1*y[i0 + 2]) + y[i0 + 5] ));
	}

	return 0;
}

int Acceleration::GravityBC(double t, double *y, double *accel)
{
	int result = GravityBC_SelfInteracting(t, y, accel);
	HANDLE_RESULT(result);

	result = GravityBC_NonSelfInteracting(t, y, accel);
	HANDLE_RESULT(result);

	return result;
}

int Acceleration::GravityBC_SelfInteracting(double t, double *y, double *accel)
{
	int	nMassive = bodyData->nBodies.NOfMassive();

	for (register int i = 0; i < nMassive; i++) {
		double rMin = 1.0e10;
		register int i0 = 6*i;

		bodyData->indexOfNN[i] = -1;
		bodyData->distanceOfNN[i] = 0.0;
		accel[i0 + 0] = y[i0 + 3]; 
		accel[i0 + 1] = y[i0 + 4]; 
		accel[i0 + 2] = y[i0 + 5];
		memset(&accel[i0 + 3], 0, 3*sizeof(double));

		// The bodies are sorted with increasing mass, therefore to
		// increase the accuracy the lightest ones are added first
		for (register int j = nMassive - 1; j >= 0; j--) {
		//for (register int j = 0; j < nMassive; j++) {
			if (j == i)
				continue;
			register int j0 = 6*j; 
			double dxij = y[j0 + 0] - y[i0 + 0];
			double dyij = y[j0 + 1] - y[i0 + 1];
			double dzij = y[j0 + 2] - y[i0 + 2];
			double rij2 = SQR(dxij) + SQR(dyij) + SQR(dzij);
			double rij  = sqrt(rij2);

			// Select the nearest neighbour
			if (rij < rMin) {
				rMin = rij;
				bodyData->indexOfNN[i] = j;
				bodyData->distanceOfNN[i] = rij;
			}
			// c = m_j/rij^3
			double c = bodyData->mass[j] * 1.0/(rij2*rij);
			accel[i0 + 3] += c*dxij;
			accel[i0 + 4] += c*dyij;
			accel[i0 + 5] += c*dzij;
		}
		accel[i0 + 3] *= Constants::Gauss2;
		accel[i0 + 4] *= Constants::Gauss2;
		accel[i0 + 5] *= Constants::Gauss2;
	}

	return 0;
}

int Acceleration::GravityBC_NonSelfInteracting(double t, double *y, double *accel)
{
	int	nMassive = bodyData->nBodies.NOfMassive();

	for (register int i = nMassive; i < bodyData->nBodies.total; i++) {
		double rMin = 1.0e10;
		register int i0 = 6*i;

		bodyData->indexOfNN[i] = -1;
		bodyData->distanceOfNN[i] = 0.0;
		accel[i0 + 0] = y[i0 + 3]; 
		accel[i0 + 1] = y[i0 + 4]; 
		accel[i0 + 2] = y[i0 + 5]; 
		// The bodies are sorted with increasing mass, therefore to
		// increase the accuracy the lightest ones are added first
		for (register int j = nMassive - 1; j >= 0; j--) {
		//for (register int j = 0; j < nMassive; j++) {
			if (j == i)
				continue;
			register int j0 = 6*j; 
			double dxij = y[j0 + 0] - y[i0 + 0];
			double dyij = y[j0 + 1] - y[i0 + 1];
			double dzij = y[j0 + 2] - y[i0 + 2];
			double rij2 = SQR(dxij) + SQR(dyij) + SQR(dzij);
			double rij  = sqrt(rij2);

			// Select the nearest neighbour
			if (rij < rMin) {
				rMin = rij;
				bodyData->indexOfNN[i] = j;
				bodyData->distanceOfNN[i] = rij;
			}
			// c = m_j/rij^3
			double c = bodyData->mass[j] * 1.0/(rij2*rij);

			accel[i0 + 3] += c*dxij;
			accel[i0 + 4] += c*dyij;
			accel[i0 + 5] += c*dzij;
		}
		accel[i0 + 3] *= Constants::Gauss2;
		accel[i0 + 4] *= Constants::Gauss2;
		accel[i0 + 5] *= Constants::Gauss2;
	}

	return 0;
}

/**
 * This is identical to the astrocentric implementation.
 */
int	Acceleration::GasDragBC(double t, double *y, double *a)
{
	return GasDragAC(t, y, a);
}

// TODO: check
int	Acceleration::MigrationTypeIBC(double t, double *y, double *accel)
{
	double factor = this->nebula->gasComponent.ReductionFactor(t);

	int lower = bodyData->nBodies.centralBody + bodyData->nBodies.giantPlanet;
	int upper = bodyData->nBodies.NOfMassive();
	for (int i=lower; i<upper;  i++) {
		if (bodyData->migType[i] != TypeI) {
			continue;
		}

		int	i0 = 6*i;
		int j0 = 3*(i-lower);
		double r2 = SQR(y[i0 + 0]) + SQR(y[i0 + 1]) + SQR(y[i0 + 2]);
		double r = sqrt(r2);
		if (r <= bodyData->migStopAt[i]) {
			accel[j0 + 0] = accel[j0 + 1] = accel[j0 + 2] = 0.0;
			bodyData->migType[i] = No;
			continue;
		}

		// rv contains the coordinate and velocity components of the body.
		double rv[6];
		memcpy(rv, &y[i0], 6*sizeof(double));
		// Transform rv to the astro centric reference frame.
		// The first 6 element of the y array is the coordinate and velocitiy components of the star.
		Calculate::TransformToAC(rv, y);
		Phase phase;
		Tools::ToPhase(rv, &phase);

		double m = this->bodyData->mass[i];
		double mc= this->bodyData->mass[0];
		double a = 0.0, e = 0.0;
		double mu = Constants::Gauss2*(mc + m);
		Ephemeris::CalculateOrbitalElement(mu, &phase, &a, &e);
		// Orbital frequency: (note, that this differs from the formula of Fogg & Nelson 2005)
		double O = Constants::Gauss * sqrt((mc + m)/CUBE(a));
		double C = SQR(mc)/(m*this->nebula->gasComponent.MidplaneDensity(r)*a*a);
		// Aspect ratio:
		double h = this->nebula->gasComponent.scaleHeight.Evaluate(r);
		double ar = h/r;
		double er = e*r;

		/*
		 *  When e > 1.1 h/r, inward migration halts as $t_{\rm m}$ becomes negative and only
		 *  resumes when eccentricity is damped to lower values. We note that under certain
		 *  circumstances, such as there being a surface density jump, or an optically thick disk,
		 *  or MHD turbulence, type I migration may be substantially modified or reversed
		 *  (Papaloizou & Nelson 2005; Paardekooper & Mellema 2006; Nelson 2005; Masset et al. 2006).
		 */
		double tm = 0.0;
		// To avoid negative values when 1.0 - FORTH(er/(1.1*h)) is negative
		if (e < 1.1*h/r) {
			tm = TypeIMigrationTime(C, O, ar, er, h);
			tm = 1.0/tm;
		}
		double te = TypeIEccentricityDampingTime(C, O, ar, er, h);
		double ti = te;
		double vr = y[i0 + 0]*y[i0 + 3] + y[i0 + 1]*y[i0 + 4] + y[i0 + 2]*y[i0 + 5];
		te = 2.0*vr/(r2*te);
		ti = 2.0/ti;

		accel[j0 + 0] = -factor*(tm*y[i0 + 3] + te*y[i0 + 0]);
		accel[j0 + 1] = -factor*(tm*y[i0 + 4] + te*y[i0 + 1]);
		accel[j0 + 2] = -factor*(tm*y[i0 + 5] + te*y[i0 + 2] + ti*y[i0 + 5]);
	}

	return 0;
}

// TODO: check
int	Acceleration::MigrationTypeIIBC(double t, double *y, double *accel)
{
	double factor = nebula->gasComponent.ReductionFactor(t);

	int lower = bodyData->nBodies.centralBody;
	int upper = bodyData->nBodies.centralBody + bodyData->nBodies.giantPlanet;
	for (int i=lower; i<upper; i++) {
		if (bodyData->migType[i] != TypeII) {
			continue;
		}

		int	i0 = 6*i;
		int j0 = 3*(i-lower);
		double r2 = SQR(y[i0 + 0]) + SQR(y[i0 + 1]) + SQR(y[i0 + 2]);
		double r = sqrt(r2);
		if (r <= bodyData->migStopAt[i]) {
			bodyData->migType[i] = No;
			accel[j0 + 0] = accel[j0 + 1] = accel[j0 + 2] = 0.0;
			continue;
		}

		// rv contains the coordinate and velocity components of the body.
		double rv[6];
		memcpy(rv, &y[i0], 6*sizeof(double));
		// Transform rv to the astro centric reference frame.
		// The first 6 element of the y array is the coordinate and velocitiy components of the star.
		Calculate::TransformToAC(rv, y);
		Phase phase;
		Tools::ToPhase(rv, &phase);

		double m = bodyData->mass[i];
		double mc= bodyData->mass[0];
		double a = 0.0, e = 0.0;
		double mu = Constants::Gauss2*(mc + m);
		Ephemeris::CalculateOrbitalElement(mu, &phase, &a, &e);
		// Orbital frequency: (note, that this differs from the formula of Fogg & Nelson 2005)
		double O = Constants::Gauss * sqrt((mc + m)/CUBE(a));

		double c0 = TauNu(r, O);
		double vr = y[i0 + 3]*y[i0 + 0] + y[i0 + 4]*y[i0 + 1] + y[i0 + 5]*y[i0 + 2];
		double c1 = vr/r2;

		accel[j0 + 0] = -factor*(c0 * ( 0.5*y[i0 + 3] + 50 * (c1*y[i0 + 0]) ));
		accel[j0 + 1] = -factor*(c0 * ( 0.5*y[i0 + 4] + 50 * (c1*y[i0 + 1]) ));
		accel[j0 + 2] = -factor*(c0 * ( 0.5*y[i0 + 5] + 50 * (c1*y[i0 + 2]) + y[i0 + 5] ));
	}

	return 0;
}

double Acceleration::TypeIMigrationTime(const double C, const double O, const double ar, const double er, const double h)
{
	double result = 0.0;

	double Cm = 2.0/(2.7 + 1.1*abs(this->nebula->gasComponent.density.index))/O;
	double er1 = er/(1.3*h);
	double er2 = er/(1.1*h);
	double frac = (1.0 + FIFTH(er1) )/(1.0 - FORTH(er2) );
	result = Cm*C*SQR(ar)*frac;

	return result;
}

double Acceleration::TypeIEccentricityDampingTime(const double C, const double O, const double ar, const double er, const double h)
{
	static const double Q = 0.78;
	double result = 0.0;

	double Ce = 0.1/(Q*O);
	double frac = 1.0 + 0.25*CUBE(er/h);
	result = Ce*C*FORTH(ar)*frac;

	return result;
}

//Vector Acceleration::CircularVelocity(double mu, double r, double alpha)
//{
//	Vector result;
//
//	double v = sqrt( mu/r );
//	result.x = -v*sin(alpha);
//	result.y = v*cos(alpha);
//	result.z = 0.0;
//
//	return result;
//}
//
//Vector Acceleration::GasVelocity(double mu, double r, double alpha)
//{
//	double f = sqrt(1.0 - 2.0*this->nebula->gasComponent.eta.Evaluate(r));
//	Vector v = CircularVelocity(mu, r, alpha);
//	v.x *= f;
//	v.y *= f;
//	v.z *= f;
//
//	return v;
//}
//
//double Acceleration::GasDensityAt(double r, double z)
//{
//	static double innerEdge = nebula->gasComponent.innerEdge;
//	static double a = nebula->gasComponent.density.Evaluate(innerEdge) / SQR(SQR(innerEdge));
//
//	double h = nebula->gasComponent.scaleHeight.Evaluate(r);
//	double arg = SQR(z/h);
//	double density = 0.0;
//	if ( r > innerEdge) {
//		density = nebula->gasComponent.density.Evaluate(r) * exp(-arg);
//	}
//	else {
//		density = a * SQR(SQR(r)) * exp(-arg);
//	}
//	return density;
//}

double Acceleration::TauNu(double r, double O)
{
	double index = nebula->gasComponent.tau.index;
	double c = nebula->gasComponent.tau.c;

	double result = 0;

	double h = nebula->gasComponent.scaleHeight.Evaluate(r);
	if (index == 2) 
	{
		result = c * SQR(r/h) / (nebula->gasComponent.alpha*O);
	}
	else 
	{
		result = c * pow(r/h, index) / (nebula->gasComponent.alpha*O);
	}

	return result;
}
