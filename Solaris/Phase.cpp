#include "Phase.h"

Phase::Phase() {
	bodyId = -1;

	this->position.x = 0.0;
	this->position.y = 0.0;
	this->position.z = 0.0;
	this->velocity.x = 0.0;
	this->velocity.y = 0.0;
	this->velocity.z = 0.0;
}

Phase::Phase(int id) {
	bodyId = id;

	this->position.x = 0.0;
	this->position.y = 0.0;
	this->position.z = 0.0;
	this->velocity.x = 0.0;
	this->velocity.y = 0.0;
	this->velocity.z = 0.0;
}

Phase::Phase(int id, double rx, double ry, double rz, double vx, double vy, double vz)
{
	bodyId = id;
	this->position.x = rx;
	this->position.y = ry;
	this->position.z = rz;
	this->velocity.x = vx;
	this->velocity.y = vy;
	this->velocity.z = vz;
}

Phase::Phase(const Phase &phase) {
	bodyId = phase.bodyId;
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
