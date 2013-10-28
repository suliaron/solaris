#ifndef NBODIES_H_
#define NBODIES_H_

#include <iostream>
#include <list>

#include "Body.h"

class NBodies {
public:

	NBodies();

	int Count(std::list<Body *> &list);
	int UpdateAfterRemove(BodyType bodyType);

	int NOfMassive();
	int NOfPlAndSpl();
	int NOfVar();

	int centralBody;
	int giantPlanet;
	int rockyPlanet;
	int protoPlanet;
	int superPlanetsimal;
	int planetsimal;
	int testParticle;
	int	total;

	int removed;
	int threshold;

	// Output streams
	friend std::ostream& operator<<( std::ostream&, NBodies nBodies);
};

#endif
