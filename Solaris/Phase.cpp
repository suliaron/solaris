#include "Phase.h"

Phase::Phase() {
	bodyId = -1;

}

Phase::Phase(int id) {
	bodyId = id;
}

Phase::Phase(int id, double rx, double ry, double rz, double vx, double vy, double vz)
{
	bodyId = id;

	position.x = rx;
	position.y = ry;
	position.z = rz;
	
	velocity.x = vx;
	velocity.y = vy;
	velocity.z = vz;
}

Phase::Phase(const Phase &phase) {
	bodyId	 = phase.bodyId;
	position = Vector(phase.position);
	velocity = Vector(phase.velocity);
}

// Write output Phase in format: "(%f, %f, %f)
//								  (%f, %f, %f)"
std::ostream& operator<<(std::ostream& output, Phase p)
{
	output << p.position << std::endl;
	output << p.velocity << std::endl;

	return output;
}
