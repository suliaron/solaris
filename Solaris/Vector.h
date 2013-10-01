#ifndef VECTOR_H_
#define VECTOR_H_

#include <iostream>
#include <cmath>

class Vector {
public:

	double	x;
	double	y;
	double	z;

	Vector();
	Vector(double a, double b, double c);
	// Copy constructor
	Vector(const Vector &vector);

	// Computes the length of the vector
	double Length() {
		return std::sqrt(x*x + y*y + z*z);
	}

	double LengthSQR() {
		return x*x + y*y + z*z;
	}

	static double DotProduct(const Vector &v, const Vector &u);
    static Vector CrossProduct(const Vector &v, const Vector &u);

	friend Vector operator+(const Vector &v1, const Vector &v2);
	friend Vector operator*(const Vector &v, double c);
	friend Vector operator*(double c, const Vector &v);

	// Input/Output streams
	friend std::istream& operator>>( std::istream&, Vector&);
	friend std::ostream& operator<<( std::ostream&, Vector);
};

#endif
