#include <cmath>

#include "BodyData.h"
#include "Calculate.h"
#include "Constants.h"
#include "Phase.h"
#include "Vector.h"

#define SQR(a)		((a)*(a))

// TODO: implement
void Calculate::TransformToAC(Phase *phase, Phase *phaseOfBC)
{
}

// TODO: implement
void Calculate::TransformToAC(double *y, double *bc)
{
}

// TODO: Check
void Calculate::TransformToBC(Phase *phase, Phase *phaseOfBC)
{
	phase->position.x -= phaseOfBC->position.x;
	phase->position.y -= phaseOfBC->position.y;
	phase->position.z -= phaseOfBC->position.z;
	phase->velocity.x -= phaseOfBC->velocity.x;
	phase->velocity.y -= phaseOfBC->velocity.y;
	phase->velocity.z -= phaseOfBC->velocity.z;
}

// TODO: Check
void Calculate::TransformToBC(double *y, double *bc)
{
	for (int j=0; j<6; j++)
		y[j] -= bc[j];
}

/**
 * Computes the total mass, the position and the velocity of the barycenter, the length of the position and velocity 
 * vectors, the angular momentum and its length, the kinetic, the potential and the total energy.
 */
int Calculate::Integrals(BodyData *bodyData)
{
	double m;
	TotalMass(bodyData, m);
	bodyData->integrals[ 0] = m;

	PhaseOfBC(bodyData, &bodyData->integrals[1]);
	bodyData->integrals[ 7] = sqrt(SQR(bodyData->integrals[1]) + SQR(bodyData->integrals[2]) + SQR(bodyData->integrals[3]));
	bodyData->integrals[ 8] = sqrt(SQR(bodyData->integrals[4]) + SQR(bodyData->integrals[5]) + SQR(bodyData->integrals[6]));
	Vector c;
	AngularMomentum(bodyData, &c);
	bodyData->integrals[ 9] = c.x;
	bodyData->integrals[10] = c.y;
	bodyData->integrals[11] = c.z;
	bodyData->integrals[12] = c.Length();
	KineticEnergy(bodyData, bodyData->integrals[13]);
	PotentialEnergy(bodyData, bodyData->integrals[14]);
	bodyData->integrals[15] = bodyData->integrals[13] - bodyData->integrals[14];

	return 1;
}

int	Calculate::TotalMass(BodyData *bodyData, double& result)
{
	result = 0.0;
	for (int i = 0; i < bodyData->nBodies.NOfMassive(); i++)
		result += bodyData->mass[i];

	return 0;
}

int	Calculate::PhaseOfBC(BodyData *bodyData, double *bc)
{
	for (int j=0; j<6; j++) {
		bc[j] = 0.0;
	}

	double M = 0.0;
	Calculate::TotalMass(bodyData, M);
	for (int i = 0; i < bodyData->nBodies.total; i++) {
		int i0 = 6*i;
		for (int j=0; j<6; j++) {
			bc[j] += bodyData->mass[i] * bodyData->y0[i0 + j];
		}
	}
	for (int j=0; j<6; j++) {
		bc[j] /= M;
	}

	return 0;
}

int	Calculate::PhaseWithRespectToBC(BodyData *bodyData, double* bc)
{
	for (int i = 0; i < bodyData->nBodies.total; i++) {
		int i0 = 6*i;
		for (int j=0; j<6; j++) {
			bodyData->y0[i0 + j] -= bc[j];
		}
	}

	return 0;
}

int	Calculate::AngularMomentum(BodyData *bodyData, Vector* result)
{
	result->x = 0.0;
	result->y = 0.0;
	result->z = 0.0;

	for (int i = 0; i < bodyData->nBodies.total; i++) {
		int i0 = 6*i;
		Vector r = Vector(bodyData->y0[i0 + 0], bodyData->y0[i0 + 1], bodyData->y0[i0 + 2]);
		Vector v = Vector(bodyData->y0[i0 + 3], bodyData->y0[i0 + 4], bodyData->y0[i0 + 5]);
		Vector c = bodyData->mass[i]*Vector::CrossProduct(r, v);
		result->x += c.x;
		result->y += c.y;
		result->z += c.z;
	}

	return 0;
}

int	Calculate::Energy(BodyData *bodyData, double& result)
{
	double kinetic = 0.0;
	double potential = 0.0;

	KineticEnergy(bodyData, kinetic);
	PotentialEnergy(bodyData, potential);

	result = kinetic - potential;

	return 0;
}

int	Calculate::PotentialEnergy(BodyData *bodyData, double& result)
{
	result = 0.0;

	for (int i = 0; i < bodyData->nBodies.total; i++) {
		int i0 = 6*i;
		for (int j = 0; j < bodyData->nBodies.total; j++) {
			if (i == j)
				continue;
			int j0 = 6*j;
			double dx = bodyData->y0[j0 + 0] - bodyData->y0[i0 + 0];
			double dy = bodyData->y0[j0 + 1] - bodyData->y0[i0 + 1];
			double dz = bodyData->y0[j0 + 2] - bodyData->y0[i0 + 2];
			double rij = sqrt(dx*dx + dy*dy + dz*dz);
			result += bodyData->mass[i]*bodyData->mass[j]/rij;
		}
	}
	result *= 0.5*Constants::Gauss2;

	return 0;
}

int	Calculate::KineticEnergy(BodyData *bodyData, double& result)
{
	result = 0.0;

	for (int i = 0; i < bodyData->nBodies.total; i++) {
		int i0 = 6*i;
		double v2 = SQR(bodyData->y0[i0 + 3]) + SQR(bodyData->y0[i0 + 4]) + SQR(bodyData->y0[i0 + 5]);
		result += 0.5*bodyData->mass[i] * v2;
	}

	return 0;
}
