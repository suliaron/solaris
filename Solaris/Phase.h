#ifndef PHASE_H_
#define PHASE_H_

#include "Vector.h"

class Phase {
public:
	Phase();
	Phase(int id);
	Phase(int id, double rx, double ry, double rz, double vx, double vy, double vz);
	Phase(const Phase &phase);

	int	bodyId;
	Vector position;
	Vector velocity;

	// Input/Output streams
	friend std::ostream& operator<<( std::ostream&, Phase);

};

#endif
