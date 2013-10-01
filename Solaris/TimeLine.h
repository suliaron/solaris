#ifndef TIMELINE_H_
#define TIMELINE_H_

class TimeLine {
public:
	TimeLine();

	bool Forward() { return length >= 0 ? true : false; }

	/// <summary>
    /// The total length of the integration.
    /// </summary>
    double length;
    /// <summary>
    /// Time interval between two subsequent save, where save means to
    /// store the intermediate result of the integration. If it is not
    /// defined, i.e. it is null, only the start and end state will be stored.
    /// This must be greater than zero in the xml, but will be a signed number
	/// during the execution of the code.
    /// </summary>       
    double output;
    /// <summary>
    /// Time of the next save event, where save means to
    /// store the Phases of all the bodies which are part of the integration
    /// at Save time.
    /// </summary>
    double save;
    /// <summary>
    /// Actual time of the integration.
    /// </summary>
    double time;
    /// <summary>
    /// The next step size to be tried by the stepper function.
    /// </summary>
    double hNext;
    /// <summary>
    /// The prevoius step size which was successfully accomplished by the stepper function.
    /// </summary>
    double hDid;

    /// <summary>
    /// Start time of the integration
    /// </summary>
    double	start;
	/// The start time was defined by the user
	bool	startTimeDefined;
};

#endif