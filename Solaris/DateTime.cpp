#include "DateTime.h"

DateTime::DateTime()
{
	Year		= 0;
	Month		= 0;
	Day			= 0;
	Hour		= 0;
	Minute		= 0;
	Second		= 0;
	Millisec	= 0;
}

DateTime::DateTime(int year, int month, int day, int hour, int min, int sec, int millisec)
{
	Year		= year;
	Month		= month;
	Day			= day;
	Hour		= hour;
	Minute		= min;
	Second		= sec;
	Millisec	= millisec;
}
