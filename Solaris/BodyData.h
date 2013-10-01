#ifndef BODYDATA_H_
#define BODYDATA_H_

#include "NBodies.h"
#include "Phase.h"
#include "Vector.h"

class BodyData {
public:

	BodyData();
	~BodyData();

	int Allocate();
	void Free();

	NBodies nBodies;

	double	time;
	double	h;

	int		*id;
	int		*type;
	int		*migType;
	double	*migStopAt;
	// The index of the nearest neighbor, i.e. the id of the jth body's nearest neighbor is: id[indexOfNN[j]]
	int		*indexOfNN;
	double	*distanceOfNN;

	double	*mass;
	double	*radius;
	double	*density;
    double  *cD;
	double	*gammaStokes;
	double	*gammaEpstein;

	double	*y0;
	double	*y;
	double	*yBetterEst;
	double	*yscale;
	double	*accel;
	double	*error;

	double	bc[6];
	Phase	phaseOfBC;
	Vector	angularMomentum[2];
	double	kineticEnergy[2];
	double	potentialEnergy[2];
	double	totalEnergy[2];
	double	integrals[16];
};

#endif
