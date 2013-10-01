#ifndef CALCULATE_H_
#define CALCULATE_H_

class BodyData;
class Phase;
class Vector;

class Calculate
{
public:
	static void TransformToAC(Phase *phase, Phase *phaseOfBC);
	static void TransformToAC(double *y, double *bc);
	static void TransformToBC(Phase *phase, Phase *phaseOfBC);
	static void TransformToBC(double *y, double *bc);

	static int	Integrals(				BodyData *bodyData);
	static int	TotalMass(				BodyData *bodyData, double& result);
	static int	PhaseOfBC(				BodyData *bodyData, double* bc);
	static int	PhaseWithRespectToBC(	BodyData *bodyData, double* bc);
	static int	AngularMomentum(		BodyData *bodyData, Vector* result);
	static int	Energy(					BodyData *bodyData, double& result);
	static int	PotentialEnergy(		BodyData *bodyData, double& result);
	static int	KineticEnergy(			BodyData *bodyData, double& result);
};

#endif