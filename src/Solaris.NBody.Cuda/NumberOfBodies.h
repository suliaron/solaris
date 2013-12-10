#pragma once

#include "InteractionBound.h"

struct	NumberOfBodies {

	NumberOfBodies(int s, int gP, int rP, int pP, int sPl, int pl, int tP) : 
		star(s), 
		giantPlanet(gP), 
		rockyPlanet(rP), 
		protoPlanet(pP), 
		superPlanetesimal(sPl), 
		planetesimal(pl), 
		testParticle(tP) 
	{
		total = star + giantPlanet + rockyPlanet + protoPlanet + superPlanetesimal + planetesimal + testParticle;
	}

	int		n_massive() {
		return star + giantPlanet + rockyPlanet + protoPlanet + superPlanetesimal + planetesimal;
	}
	int		n_self_interacting() {
		return star + giantPlanet + rockyPlanet + protoPlanet;
	}

	inline InteractionBound get_self_interacting() {
		sink.x		= 0;
		sink.y		= n_self_interacting();
		source.x	= 0;
		source.y	= n_massive();
		InteractionBound iBound(sink, source);
		return iBound;
	}

	inline InteractionBound get_nonself_interacting() {
		sink.x			= n_self_interacting();
		sink.y			= n_massive();
		source.x		= 0;
		source.y		= n_self_interacting();
		InteractionBound iBound(sink, source);
		return iBound;
	}

	inline InteractionBound get_non_interacting() {
		sink.x			= n_massive();
		sink.y			= total;
		source.x		= 0;
		source.y		= n_massive();
		InteractionBound iBound(sink, source);
		return iBound;
	}

	inline InteractionBound get_bodies_gasdrag() {
		sink.x			= n_self_interacting();
		sink.y			= n_massive();
		source.x		= 0;
		source.y		= 0;
		InteractionBound iBound(sink, source);
		return iBound;
	}

	int		star;
	int		giantPlanet;
	int		rockyPlanet;
	int		protoPlanet;
	int		superPlanetesimal;
	int		planetesimal;
	int		testParticle;
	int		total;

private:
	int2_t		sink;
	int2_t		source;

};
