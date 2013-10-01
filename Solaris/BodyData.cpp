#include <cstdio>

#include "BodyData.h"
#include "Error.h"

BodyData::BodyData()
{
	time		 = 0.0;
	h			 = 0.0;

	id			 = 0;
	type		 = 0;
	migType		 = 0;
	migStopAt    = 0;
 	indexOfNN	 = 0;
	distanceOfNN = 0;

	mass		 = 0;
	radius		 = 0;
	density		 = 0;
    cD           = 0;
	gammaStokes	 = 0;
	gammaEpstein = 0;

	y0			 = 0;
	y			 = 0;
	yBetterEst	 = 0;
	yscale		 = 0;
	accel		 = 0;
	error		 = 0;

	bc[0] = bc[1] = bc[2] = bc[3] = bc[4] = bc[5] = 0.0;
	phaseOfBC.position.x = phaseOfBC.position.y = phaseOfBC.position.z = 0.0;
	phaseOfBC.velocity.x = phaseOfBC.velocity.y = phaseOfBC.velocity.z = 0.0;

	angularMomentum[0].x = angularMomentum[0].y = angularMomentum[0].z = 0.0;
	angularMomentum[1].x = angularMomentum[1].y = angularMomentum[1].z = 0.0;
	kineticEnergy[0] = kineticEnergy[1] = 0.0;
	potentialEnergy[0] = potentialEnergy[1] = 0.0;
	totalEnergy[0] = totalEnergy[1] = 0.0;
}

BodyData::~BodyData()
{
	delete[] id;
	delete[] type;
	delete[] migType;
	delete[] migStopAt;
	delete[] indexOfNN;
	delete[] distanceOfNN;

	delete[] mass;
	delete[] radius;
	delete[] density;
    delete[] cD;
	delete[] gammaStokes;
	delete[] gammaEpstein;
	delete[] y0;
	delete[] y;
	delete[] yBetterEst;
	delete[] yscale;
	delete[] accel;
	delete[] error;
}

int BodyData::Allocate()
{
	if (nBodies.total == 0) {
		Error::_errMsg = "host memory allocation";
		Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
		return 1;
	}

	id			= new int[nBodies.total];
	if (id == 0) {
		Error::_errMsg = "host memory allocation";
		Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
		return 1;
	}
	type		= new int[nBodies.total];
	if (type == 0) {
		Error::_errMsg = "host memory allocation";
		Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
		return 1;
	}
	migType		= new int[nBodies.total];
	if (migType == 0) {
		Error::_errMsg = "host memory allocation";
		Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
		return 1;
	}
	migStopAt	= new double[nBodies.total];
	if (migStopAt == 0) {
		Error::_errMsg = "host memory allocation";
		Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
		return 1;
	}
	indexOfNN	= new int[nBodies.total];
	if (indexOfNN == 0) {
		Error::_errMsg = "host memory allocation";
		Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
		return 1;
	}
	distanceOfNN	= new double[nBodies.total];
	if (distanceOfNN == 0) {
		Error::_errMsg = "host memory allocation";
		Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
		return 1;
	}

	mass		= new double[nBodies.total];
	if (mass == 0) {
		Error::_errMsg = "host memory allocation";
		Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
		return 1;
	}
	radius		= new double[nBodies.total];
	if (radius == 0) {
		Error::_errMsg = "host memory allocation";
		Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
		return 1;
	}
	density		= new double[nBodies.total];
	if (density == 0) {
		Error::_errMsg = "host memory allocation";
		Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
		return 1;
	}
	cD  		= new double[nBodies.total];
	if (cD == 0) {
		Error::_errMsg = "host memory allocation";
		Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
		return 1;
	}
	gammaStokes	= new double[nBodies.total];
	if (gammaStokes == 0) {
		Error::_errMsg = "host memory allocation";
		Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
		return 1;
	}
	gammaEpstein= new double[nBodies.total];
	if (gammaEpstein == 0) {
		Error::_errMsg = "host memory allocation";
		Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
		return 1;
	}

	int nVar = nBodies.NOfVar();
	y0			= new double[nVar];
	if (y0 == 0) {
		Error::_errMsg = "host memory allocation";
		Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
		return 1;
	}
	y			= new double[nVar];
	if (y == 0) {
		Error::_errMsg = "host memory allocation";
		Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
		return 1;
	}
	yBetterEst	= new double[nVar];
	if (yBetterEst == 0) {
		Error::_errMsg = "host memory allocation";
		Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
		return 1;
	}
	yscale		= new double[nVar];
	if (y == 0) {
		Error::_errMsg = "host memory allocation";
		Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
		return 1;
	}
	accel		= new double[nVar];
	if (y == 0) {
		Error::_errMsg = "host memory allocation";
		Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
		return 1;
	}
	error		= new double[nVar];
	if (y == 0) {
		Error::_errMsg = "host memory allocation";
		Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
		return 1;
	}

	return 0;
}

void BodyData::Free()
{
	delete[] id;
	delete[] type;
	delete[] migType;
	delete[] migStopAt;
	delete[] indexOfNN;
	delete[] distanceOfNN;

	delete[] mass;
	delete[] radius;
	delete[] density;
    delete[] cD;
	delete[] gammaStokes;
	delete[] gammaEpstein;
	delete[] y0;
	delete[] y;
	delete[] yBetterEst;
	delete[] yscale;
	delete[] accel;
	delete[] error;
}
