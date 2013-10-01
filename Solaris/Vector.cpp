#include "Vector.h"

Vector::Vector() { 
	x = 0.0;
	y = 0.0;
	z = 0.0;
}

Vector::Vector(double a, double b, double c) {
	x = a;
	y = b;
	z = c;
}

// Copy constructor
Vector::Vector(const Vector &vector) {
	x = vector.x;
	y = vector.y;
	z = vector.z;
}

double Vector::DotProduct(const Vector &v, const Vector &u) {
    return v.x*u.x + v.y*u.y + v.z*u.z;
}

Vector Vector::CrossProduct(const Vector &v, const Vector &u) {
    Vector result;

    result.x = v.y*u.z - v.z*u.y;
    result.y = v.z*u.x - v.x*u.z;
    result.z = v.x*u.y - v.y*u.x;

    return result;
}

Vector operator+(const Vector &v1, const Vector &v2) {
	return Vector(v1.x + v2.x, v1.y + v2.y, v1.z + v2.z);
}

Vector operator*(const Vector &v, double c) {
	return Vector(c * v.x, c * v.y, c * v.z);
}

Vector operator*(double c, const Vector &v) {
	return Vector(c * v.x, c * v.y, c * v.z);
}

//------------------------------------------------------------------
// IO streams
//------------------------------------------------------------------

// Read input Vector format: "(%f)", "(%f, %f)", or "(%f, %f, %f)"
std::istream& operator>>(std::istream& input, Vector& P)
{
	char c;
	input >> c;                // skip '('
	input >> P.x;
	input >> c;                
	if (c == ')') {
		return input;
	}
	// else                    // skip ','
	input >> P.y;
	input >> c;
	if (c == ')') {
		return input;
	}
	// else                    // skip ','
	input >> P.z;
	input >> c;                // skip ')'
	return input;
}

// Write output Vector in format: "(%f)", "(%f, %f)", or "(%f, %f, %f)"
std::ostream& operator<<(std::ostream& output, Vector P)
{
	output << "(" << P.x << ", " << P.y << ", " << P.z << ")";
	return output;
}
