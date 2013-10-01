#ifndef ACCELERATION_H_
#define ACCELERATION_H_

#include "IntegratorType.h"

class BodyData;
class Nebula;
class Vector;

/**
 * This class contains the accelerations emanating from various forces.
 */

class Acceleration {
public:
	Acceleration(IntegratorType iType, bool baryCentric, BodyData *bD, Nebula *n);
	~Acceleration();

	int	Compute(            double t, double *y, double *totalAccel);
	int ComputeAstroCentric(double t, double *y, double *totalAccel);
	int ComputeBaryCentric( double t, double *y, double *totalAccel);

	int	GravityAC(          double t, double *y, double *a);
	int	GasDragAC(          double t, double *y, double *a);
	int MigrationTypeIAC(   double t, double *y, double *a);
	int MigrationTypeIIAC(  double t, double *y, double *a);

	int	GravityBC(                   double t, double *y, double *a);
	int	GravityBC_SelfInteracting(   double t, double *y, double *a);
	int GravityBC_NonSelfInteracting(double t, double *y, double *a);
	int	GasDragBC(                   double t, double *y, double *a);
	int MigrationTypeIBC(            double t, double *y, double *a);
	int MigrationTypeIIBC(           double t, double *y, double *a);

	int	GravityBC_RKDP(                   double t, double *y, double *a);
	int	GravityBC_RKDP_SelfInteracting(   double t, double *y, double *a);
	int GravityBC_RKDP_NonSelfInteracting(double t, double *y, double *a);
/*	int	GasDragBC_RKDP(                   double t, double *y, double *a);
	int MigrationTypeIBC_RKDP(            double t, double *y, double *a);
	int MigrationTypeIIBC_RKDP(           double t, double *y, double *a);*/

	double TypeIMigrationTime(          const double C, const double O, const double ar, const double er, const double h);
	double TypeIEccentricityDampingTime(const double C, const double O, const double ar, const double er, const double h);

	Vector CircularVelocity(double mu, double r, double alpha);
	Vector GasVelocity(double mu, double r, double alpha);
	double GasDensityAt(double r, double z);
	double TauNu(double r, double O);

	BodyData	*bodyData;
	Nebula		*nebula;

	double		*accelGasDrag;
	double		*accelMigrationTypeI;
	double		*accelMigrationTypeII;

	// The inverse of the cube of the distance from the central body
	// nBodies contains the number of the different bodies
	double		*rm3;

	bool evaluateGasDrag;
	bool evaluateTypeIMigration;
	bool evaluateTypeIIMigration;

private:
	IntegratorType	_integratorType;
	bool			_baryCentric;
};

#endif
