#ifndef COMPONENT_H_
#define COMPONENT_H_

#include <string>

class Component
{
public:
	Component();
	Component(std::string name, double ratio);

	std::string name;
	double		ratio;
};

#endif
