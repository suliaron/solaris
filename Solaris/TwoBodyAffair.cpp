#include <cstring>
#include <iostream>

#include "Constants.h"
#include "TwoBodyAffair.h"

int TwoBodyAffair::_eventId = 0;

TwoBodyAffair::TwoBodyAffair(EventType type, double time, int idx1, int idx2, int body1Id, int body2Id, double *body1Phase, double *body2Phase)
{
	id = TwoBodyAffair::_eventId++;

	this->type    = type;
	this->time    = time;
	this->idx1    = idx1;
	this->idx2    = idx2;
	this->body1Id = body1Id;
	this->body2Id = body2Id;
	memcpy(this->body1Phase, body1Phase, 6*sizeof(double));
	memcpy(this->body2Phase, body2Phase, 6*sizeof(double));
}

std::ostream& operator<<(std::ostream& output, const TwoBodyAffair& affair)
{
	switch (affair.type) {
		case Ejection:
			output << "Ejection: At " << affair.time * Constants::DayToYear << " [yr] body with id: " << affair.body2Id;
			break;
		case HitCentrum:
			output << "HitCentrum: At " << affair.time * Constants::DayToYear << " [yr] body with id: " << affair.body2Id;
			break;
		default:
			break;
	}

	return output;
}
