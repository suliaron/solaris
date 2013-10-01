#ifndef TWOBODYAFFAIR_H_
#define TWOBODYAFFAIR_H_

#include <iostream> // to avoid the "ISO C++ forbids declaration of ‘ostream’ with no type" error

enum EventType
{
    Ejection		= 0,
    HitCentrum		= 1,
    CloseEncounter	= 2,
    Collision		= 3,
    WeakCapture		= 4
};

class TwoBodyAffair
{
public:
	TwoBodyAffair(EventType type, double time, int idx1, int idx2, int body1Id, int body2Id, double *body1Phase, double *body2Phase);

	int			id;
	EventType	type;
	double		time;
	int			idx1;
	int			idx2;
	int			body1Id;
	int			body2Id;
	double		body1Phase[6];
	double		body2Phase[6];

	friend std::ostream& operator<<(std::ostream& output, const TwoBodyAffair& affair);

private:
	static int _eventId;
};

#endif
