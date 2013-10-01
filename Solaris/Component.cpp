#include "Component.h"

Component::Component() 
{
	ratio = 0.0;
}

Component::Component(std::string n, double r)
{
	name = n;
	ratio = r;
}
