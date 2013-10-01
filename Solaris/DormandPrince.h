#ifndef DORMAND_PRINCE_H_
#define DORMAND_PRINCE_H_

#include <string>

class Acceleration;
class BodyData;
class TimeLine;

class DormandPrince {
public:
	DormandPrince();

	int			Driver(BodyData *bodyData, Acceleration *acceleration, TimeLine *timeLine);
	int 		Step(  BodyData *bodyData, Acceleration *acceleration);
	double		GetErrorMax(const int n, const double *yerr);

	std::string	name;
	std::string	reference;
	double		accuracy;
	double		epsilon;

private:
	int		maxIter;
	int		sizeHeightRKD;
	double	D1[9];
	double	D2[9];
	double	D3[9];
	double	D4[9];
	double	D[9][8];
};

#endif
