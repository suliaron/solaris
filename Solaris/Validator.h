#ifndef VALIDATOR_H_
#define VALIDATOR_H_

class Validator
{
public:
	static bool GreaterThan(double limit, double value);
	static bool GreaterThanOrEqualTo(double limit, double value);
	static bool LessThan(double limit, double value);
	static bool LessThanOrEqualTo(double limit, double value);
	static bool ElementOfAndContainsEndPoints(double lower, double upper, double value);
	static bool ElementOf(double lower, double upper, double value);
	static bool ElementOfAndContainsLower(double lower, double upper, double value);
	static bool ElementOfAndContainsUpper(double lower, double upper, double value);
};

#endif