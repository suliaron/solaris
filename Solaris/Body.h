#ifndef BODY_H_
#define BODY_H_

#include <string>
#include "Constants.h"
#include "Phase.h"
#include "OrbitalElement.h"
#include "Characteristics.h"

class Phase;
class OrbitalElement;
class Characteristics;

enum BodyType {
	UndefinedBodyType,
	CentralBody,
    GiantPlanet,
    RockyPlanet,
    ProtoPlanet,
    SuperPlanetesimal,
    Planetesimal,
    TestParticle,
	NOfBodyType
};

enum Ln {
	UndefinedLn,
	L1,
	L2,
	L3,
	L4,
	L5
};

enum MigrationType {
	No,
	TypeI,
	TypeII
};

enum MPCOrbitType {
    UndefinedMPCOrbitType = 0,
    Aten = 2,
    Apollo = 3,
    Amor = 4,
    ObjectWithqLt1_665 = 5,
    Hungaria = 6,
    Phocaea = 7,
    Hilda = 8,
    JupiterTrojan = 9,
    Centaur = 10,
    Plutino = 14,
    OtherResonantTNO = 15,
    Cubewano = 16,
    ScatteredDisk = 17,
    ObjectIsNEO = 2048,
    ObjectIs1kmOrLargerNEO = 4096,
    OneOppositionObjectSeenAtEarlierOpposition = 8192,
    CriticalListNumberedObject = 16384,
    ObjectIsPHA = 32768
};

class Body {
public:

	Body();
	Body(int id);
	Body(BodyType type);
	Body(BodyType type, std::string name);

	inline int GetId() { return _id; }

	double GetGm();
	double CalculateOrbitalPeriod(double mu);

	void Print();

	MPCOrbitType	mPCOrbitType;
	std::string		mPCOrbitTypestr;
	std::string		name;
	std::string		designation;
	std::string		provisionalDesignation;
	std::string		reference;
	std::string		opposition;
	std::string		guid;
	BodyType		type;
	Ln				ln;
	MigrationType	migrationType;
	double			migrationStopAt;

	Phase			*phase;
	OrbitalElement	*orbitalElement;
	Characteristics	*characteristics;

private:
	// id of the body instance
	int			_id;

	// id of the class Body, its value will be assigned to the next instance of Body
	// guaranteeing the uniqueness of the _id field of each Body object.
	static int _bodyId;
};

#endif
