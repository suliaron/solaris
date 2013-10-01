#include "Validator.h"

bool Validator::GreaterThan(double limit, double value)
{
    if (value > limit)
    {
        return true;
    }
	return false;
}

bool Validator::GreaterThanOrEqualTo(double limit, double value)
{
    if (value >= limit)
    {
        return true;
    }
	return false;
}

bool Validator::LessThan(double limit, double value)
{
    if (value < limit)
    {
        return true;
    }
	return false;
}

bool Validator::LessThanOrEqualTo(double limit, double value)
{
    if (value <= limit)
    {
        return true;
    }
	return false;
}

bool Validator::ElementOfAndContainsEndPoints(double lower, double upper, double value)
{
    if (value >= lower && value <= upper)
    {
        return true;
    }
	return false;
}

bool Validator::ElementOf(double lower, double upper, double value)
{
    if (value > lower && value < upper)
    {
        return true;
    }
	return false;
}

bool Validator::ElementOfAndContainsLower(double lower, double upper, double value)
{
    if (value >= lower && value < upper)
    {
        return true;
    }
	return false;
}

bool Validator::ElementOfAndContainsUpper(double lower, double upper, double value)
{
    if (value > lower && value <= upper)
    {
        return true;
    }
	return false;
}
