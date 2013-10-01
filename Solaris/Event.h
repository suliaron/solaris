#ifndef EVENT_H_
#define EVENT_H_

#include <list>
#include "TwoBodyAffair.h"

class Event {
public:
	Event();
	bool Contains(int body2Id);

	std::list<TwoBodyAffair> items;
	// Cumulative counter of the specific event
	int		N;
};

#endif