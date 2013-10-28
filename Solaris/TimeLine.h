#ifndef TIMELINE_H_
#define TIMELINE_H_

#include <iostream>

class TimeLine {
public:
	TimeLine();

	bool Forward() { return length >= 0 ? true : false; }

    /// Start time of the integration
    double	start;
    /// The total length of the integration.
    double	length;
    /// Time interval between two subsequent save, where save means to
    /// store the intermediate result of the integration.
    double	output;

	/// Time of the next save operation.
    double	save;
    /// Actual time of the integration.
    double	time;
	/// The number of milleniums, i.e. the number of thousends years
	int		millenium;

    /// The next step size to be tried by the stepper function.
    double	hNext;
    /// The prevoius step size which was successfully accomplished by the stepper function.
    double	hDid;

	/// The sum of all the steps since the start of the simulation (including the synchronization, pre- and main phase.
	double	elapsedTime;
	/// The elapsed (simulation) time after the last save operation.
	double	lastSave;
	/// The sum of the last NSTEP steps, only for debugging purposes.
	double	lastNSteps;

	/// The start time was defined by the user
	bool	startTimeDefined;

	// Output
	friend std::ostream& operator<<(std::ostream& output, TimeLine timeLine);
};

#endif