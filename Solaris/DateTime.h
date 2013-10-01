#ifndef DATETIME_H_
#define DATETIME_H_

class DateTime
{
public:
	DateTime();
	DateTime(int year, int month, int day, int hour, int min, int sec, int millisec);

	int Year;
	int Month;
	int Day;
	int Hour;
	int Minute;
	int Second;
	int Millisec;
};

#endif