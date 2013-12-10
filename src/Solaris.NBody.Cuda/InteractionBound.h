#pragma once

#include "config.h"

struct	InteractionBound {
	int2_t	sink;
	int2_t	source;

	InteractionBound(int2_t _sink, int2_t _source) : 
		sink(_sink),
		source(_source) 
	{ }
};
