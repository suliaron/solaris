#include <stdio.h>  //included for Eclipse

#include "NBodies.h"
#include "Error.h"

NBodies::NBodies()
{
	centralBody		 = 0;
	giantPlanet		 = 0;
	rockyPlanet		 = 0;
	protoPlanet		 = 0;
	superPlanetsimal = 0;
	planetsimal		 = 0;
	testParticle	 = 0;
	total			 = 0;

	removed			 = 0;
	threshold		 = 10;
}


int NBodies::NOfMassive()
{
	return centralBody + giantPlanet + rockyPlanet + protoPlanet;
}

int NBodies::NOfPlAndSpl()
{
	return superPlanetsimal + planetsimal;
}

int NBodies::NOfVar()
{
	return 6*total;
}

int NBodies::Count(std::list<Body *> &list)
{
	centralBody		 = 0;
	giantPlanet		 = 0;
	rockyPlanet		 = 0;
	protoPlanet		 = 0;
	superPlanetsimal = 0;
	planetsimal		 = 0;
	testParticle	 = 0;
	total			 = 0;
	for (std::list<Body *>::iterator it = list.begin(); it != list.end(); it++) {
		total++;
		switch ((*it)->type) {
			case CentralBody:
				centralBody++;
				break;
			case GiantPlanet:
				giantPlanet++;
				break;
			case RockyPlanet:
				rockyPlanet++;
				break;
			case ProtoPlanet:
				protoPlanet++;
				break;
			case SuperPlanetesimal:
				superPlanetsimal++;
				break;
			case Planetesimal:
				planetsimal++;
				break;
			case TestParticle:
				testParticle++;
				break;
			default:
				Error::_errMsg = "Unknown or undefined Body Type!";
				Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
				return 1;
		}
	}

	return 0;
}

int NBodies::UpdateAfterRemove(BodyType bodyType)
{
	switch (bodyType) {
		case CentralBody:
			centralBody--;
			break;
		case GiantPlanet:
			giantPlanet--;
			break;
		case RockyPlanet:
			rockyPlanet--;
			break;
		case ProtoPlanet:
			protoPlanet--;
			break;
		case SuperPlanetesimal:
			superPlanetsimal--;
			break;
		case Planetesimal:
			planetsimal--;
			break;
		case TestParticle:
			testParticle--;
			break;
		default:
			Error::_errMsg = "Unknown or undefined Body Type!";
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			return 1;
	}
	total--;
	removed++;

	return 0;
}
