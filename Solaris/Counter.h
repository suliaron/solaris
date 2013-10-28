#ifndef COUNTER_H_
#define COUNTER_H_

#include <iostream>

class Counter
{
public:
	Counter();

	unsigned long long int succededStep;
	unsigned long long int failedStep;

	unsigned long long int ejection;
	unsigned long long int hitCentrum;
	unsigned long long int collision;

	// Output streams
	friend std::ostream& operator<<( std::ostream& output, Counter counter);
};

#endif
