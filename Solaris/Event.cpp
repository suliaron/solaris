#include "Event.h"

Event::Event()
{
	this->N = 0;
}

bool Event::Contains(int body2Id)
{
	for (std::list<TwoBodyAffair>::iterator it = items.begin(); it != items.end(); it++) {
		if (it->body2Id == body2Id) {
			return true;
		}
	}
	return false;
}

