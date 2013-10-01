#ifndef XMLFILEADAPTER_H_
#define XMLFILEADAPTER_H_

#include <sstream>

class BodyGroup;
class BodyGroupList;
class Characteristics;
class Component;
class DragCoefficient;
class EventCondition;
class GasComponent;
class Integrator;
class Nebula;
class OrbitalElement;
class Output;
class Phase;
class PowerLaw;
class Settings;
class Simulation;
class SolidsComponent;
class TimeLine;
class Vector;

#include <list>
#include "Body.h"
#include "tinyxml.h"

class XmlFileAdapter
{
public:
	XmlFileAdapter(void);
	XmlFileAdapter(char *inputPath);

	~XmlFileAdapter(void);

	static int Load(char* path, TiXmlDocument &doc);

	static int DeserializeSimulation(TiXmlDocument &doc, Simulation &simulation);
	static int DeserializeSimulationAttributes(TiXmlAttribute *attribute, Simulation &simulation);

	static int DeserializeSettings(TiXmlElement *xmlElement, Settings *settings);
	static int DeserializeSettingsAttributes(TiXmlAttribute *attribute, Settings *settings);

	static int DeserializeOutput(TiXmlElement *xmlElement, Output *output);

	static int DeserializeIntegrator(TiXmlElement *xmlElement, Integrator *integrator);
	static int DeserializeIntegratorAttributes(TiXmlAttribute *attribute, Integrator *integrator);

	static int DeserializeAccuracy(TiXmlElement *xmlElement, double *accuracy);
	static int DeserializeAccuracyAttributes(TiXmlAttribute *attribute, double *accuracy);

	static int DeserializeTimeLine(TiXmlElement *xmlElement, TimeLine *timeLine);

	static int DeserializeEventCondition(TiXmlElement *xmlElement, EventCondition *eventCondition);
	static int DeserializeEventConditionAttributes(TiXmlAttribute *attribute, EventCondition *eventCondition);

	static int DeserializeDimensionalQuantity(TiXmlElement *xmlElement, double* value, std::string* unit);

	static int DeserializeBodyGroupList(TiXmlElement *xmlElement, BodyGroupList *bodyGroupList);
	static int DeserializeBodyGroup(TiXmlElement *xmlElement, BodyGroup *bodyGroup);
	static int DeserializeBodyGroupAttributes(TiXmlAttribute *attribute, BodyGroup *bodyGroup);

	static int DeserializeBody(TiXmlElement *xmlElement, Body *body);
	static int DeserializeBodyAttributes(TiXmlAttribute *attribute, Body *body);
	static int DeserializeBodyType(std::string type, Body *body);
	static int DeserializeLn(std::string ln, Body *body);
	static int DeserializeMigrationType(TiXmlElement *xmlElement, Body *body);
	static int DeserializeMigrationType(std::string migrationType, Body *body);
	static int DeserializeMigrationTypeAttributes(TiXmlAttribute *attribute, Body *body);
	static int DeserializePhase(TiXmlElement *xmlElement, Phase *phase);
	static int DeserializePosition(TiXmlElement *xmlElement, Vector& r);
	static int DeserializeVelocity(TiXmlElement *xmlElement, Vector& v);
	static int DeserializeVector(TiXmlElement *xmlElement, Vector& vector);
	static int DeserializeOrbitalElement(TiXmlElement *xmlElement, OrbitalElement *oe);
	static int DeserializeCharacteristics(TiXmlElement *xmlElement, Characteristics *characteristics, BodyType bodyType);
	static int DeserializeCharacteristicsAttributes(TiXmlAttribute *attribute, Characteristics *characteristics);
	static int DeserializeDragCoefficient(TiXmlElement *xmlElement, DragCoefficient *dragCoefficient);
	static int DeserializeDragCoefficientAttributes(TiXmlAttribute *attribute, DragCoefficient *dragCoefficient);
	static int DeserializeComponentList(TiXmlElement *xmlElement, std::list<Component> *componentList);
	static int DeserializeComponent(TiXmlElement *xmlElement, Component *component);
	static int DeserializeComponentAttributes(TiXmlAttribute *attribute, Component *component);

	static int DeserializeNebula(TiXmlElement *xmlElement, Nebula *nebula);
	static int DeserializeNebulaAttributes(TiXmlAttribute *attribute, Nebula *nebula);
	static int DeserializeSolidsComponent(TiXmlElement *xmlElement, SolidsComponent *solidsComponent);
	static int DeserializeSolidsDensityFunction(TiXmlElement *xmlElement, SolidsComponent *solidsComponent);
	static int DeserializeSolidsDensityFunctionAttributes(TiXmlAttribute *attribute, SolidsComponent *solidsComponent);

	static int DeserializeGasComponent(TiXmlElement *xmlElement, GasComponent *gasComponent);
	static int DeserializeGasComponentAttributes(TiXmlAttribute *attribute, GasComponent *gasComponent, std::string unitString);
	static int DeserializeGasDecreaseType(std::string type, GasComponent *gasComponent);
	static int DeserializeGasDensityFunction(TiXmlElement *xmlElement, GasComponent *gasComponent);
	static int DeserializeGasDensityFunctionAttributes(TiXmlAttribute *attribute, GasComponent *gasComponent);

	static int DeserializePowerLaw(TiXmlElement *xmlElement, PowerLaw *powerLaw);

	static bool ValidDistanceUnitAttribute(TiXmlElement *xmlElement, std::string& unitValue, std::string unitName);
	static bool ValidAngleUnitAttribute(TiXmlElement *xmlElement, std::string& unitValue, std::string unitName);
	static bool ValidMassUnitAttribute(TiXmlElement *xmlElement, std::string& unitValue, std::string unitName);
	static bool ValidTimeUnitAttribute(TiXmlElement *xmlElement, std::string& unitValue, std::string unitName);
	static bool ValidVelocityUnitAttribute(TiXmlElement *xmlElement, std::string& unitValue, std::string unitName);
	static bool ValidSurfaceDensityUnitAttribute(TiXmlElement *xmlElement, std::string& unitValue, std::string unitName);
	static bool ValidVolumeDensityUnitAttribute(TiXmlElement *xmlElement, std::string& unitValue, std::string unitName);

	TiXmlDocument doc;

private:
	static std::ostringstream _stream; /**< stream used for the conversion */

	std::string _inputPath;
};

#endif
