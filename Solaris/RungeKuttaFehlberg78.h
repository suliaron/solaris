/*
 * RungeKuttaFehlberg.h
 *
 *  Created on: 2013.04.03.
 *      Author: Dr. Áron Süli (suli)
 */

#ifndef RUNGEKUTTAFEHLBERG78_H_
#define RUNGEKUTTAFEHLBERG78_H_

#include <string>

class Acceleration;
class BodyData;
class TimeLine;

class RungeKuttaFehlberg78 {
public:

	RungeKuttaFehlberg78();

	int			Driver(BodyData *bodyData, Acceleration *acceleration, TimeLine *timeLine);
	//int 		Step(  BodyData *bodyData, Acceleration *acceleration, double *accel, double t, double h, double *yout, double *yerr);
	int 		Step(  BodyData *bodyData, Acceleration *acceleration);
	double		GetErrorMax(const int n, const double *yerr, const double *yscale);

	std::string	name;
	std::string	reference;
	double		accuracy;
	double		epsilon;

private:
	double 	D1_0, D1_1, D1_2, D1_3, D1_4, D1_5,	D1_6, D1_7, D1_8, D1_9, D1_10, D1_11, D1_12;
	double 	D_1_0, D_2_0, D_3_0, D_4_0, D_5_0, D_6_0, D_7_0, D_8_0, D_9_0, D_10_0, D_11_0, D_12_0;
	double	D_2_1;
	double	D_3_2, D_4_2;
	double	D_4_3, D_5_3, D_6_3,   D_8_3,  D_9_3,  D_10_3, D_12_3;
	double	D_5_4, D_6_4, D_7_4,   D_8_4,  D_9_4,  D_10_4, D_12_4;
	double	D_6_5, D_7_5, D_8_5,   D_9_5,  D_10_5, D_11_5, D_12_5;
	double	D_7_6, D_8_6, D_9_6,   D_10_6, D_11_6, D_12_6;
	double	D_8_7, D_9_7, D_10_7,  D_11_7, D_12_7;
	double	D_9_8, D_10_8, D_11_8, D_12_8;
	double	D_10_9,D_11_9, D_12_9;
	double	D_12_11;
};

#endif



