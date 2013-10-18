#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <sstream>

#include "XmlFileAdapter.h"

#include "Body.h"
#include "BodyGroupList.h"
#include "Constants.h"
#include "DragCoefficient.h"
#include "Ephemeris.h"
#include "Error.h"
#include "EventCondition.h"
#include "GasDecreaseType.h"
#include "Integrator.h"
#include "Nebula.h"
#include "Output.h"
#include "Settings.h"
#include "Simulation.h"
#include "TimeLine.h"
#include "Tools.h"
#include "Units.h"
#include "Validator.h"

std::ostringstream XmlFileAdapter::_stream;

XmlFileAdapter::XmlFileAdapter(void) { }

XmlFileAdapter::XmlFileAdapter(char *inputPath)
{
	_inputPath = inputPath;
}

XmlFileAdapter::~XmlFileAdapter(void) { }

int XmlFileAdapter::Load(char* path, TiXmlDocument &doc) {

	bool loadOkay = doc.LoadFile(path);

	if (!loadOkay) {
		Error::_errMsg = doc.ErrorDesc();
		if (doc.ErrorRow() > 0) {
			// insert the textual representation of row and col in the characters in the stream
			_stream << "! Row, col: [" << doc.ErrorRow() << ", " << doc.ErrorCol() << "]";
			Error::_errMsg += _stream .str() + "!";
		}
		Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
		return 1;
	}
	else {
		std::cerr << "The '" << path << "' file was successfully read from disk." << std::endl;
	}
	return 0;
}

int XmlFileAdapter::DeserializeSimulation(TiXmlDocument &doc, Simulation &simulation)
{
	TiXmlNode *root = doc.FirstChild("Simulation");
	if (root == 0) {
		Error::_errMsg = "Missing Simulation tag";
		Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
		return 1;
	}

	TiXmlElement *xmlElement = root->ToElement();
	if (xmlElement == 0) {
		_stream << "Invalid xml element at row: " << xmlElement->Row() << ", col: " << xmlElement->Column();
		Error::_errMsg = _stream.str();
		Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
		return 1;
	}
	// Iterate over the attributes
	for (TiXmlAttribute *attribute = xmlElement->FirstAttribute(); attribute; attribute = attribute->Next() ) {
		if (DeserializeSimulationAttributes(attribute, simulation) == 1)
			return 1;
	}

	TiXmlNode *node = root->FirstChild("Settings");
	if (node == 0) {
		Error::_errMsg = "Missing Settings tag";
		Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
		return 1;
	}
	xmlElement = node->ToElement();
	if (xmlElement == 0) {
		_stream << "Invalid xml element at row: " << xmlElement->Row() << ", col: " << xmlElement->Column();
		Error::_errMsg = _stream.str();
		Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
		return 1;
	}
	if (DeserializeSettings(xmlElement, &simulation.settings) == 1) {
		return 1;
	}
	
	node = root->FirstChild("BodyGroupList");
	if (node == 0) {
		Error::_errMsg = "Missing BodyGroupList tag";
		Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
		return 1;
	}
	xmlElement = node->ToElement();
	if (xmlElement == 0) {
		_stream << "Invalid xml element at row: " << xmlElement->Row() << ", col: " << xmlElement->Column();
		Error::_errMsg = _stream.str();
		Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
		return 1;
	}
	if (DeserializeBodyGroupList(xmlElement, &simulation.bodyGroupList) == 1) {
		return 1;
	}
	
	node = root->FirstChild("Nebula");
	if (node != 0) {
		xmlElement = node->ToElement();
		if (xmlElement == 0) {
			_stream << "Invalid xml element at row: " << xmlElement->Row() << ", col: " << xmlElement->Column();
			Error::_errMsg = _stream.str();
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			return 1;
		}
		simulation.nebula = new Nebula();
		if (DeserializeNebula(xmlElement, simulation.nebula) == 1) {
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			return 1;
		}		
	}

	return 0;
}

int XmlFileAdapter::DeserializeSimulationAttributes(TiXmlAttribute *attribute, Simulation &simulation)
{
	std::string attributeName = attribute->Name();
	std::string attributeValue = attribute->Value();
	Tools::Trim(attributeValue);
	// Transform the name to lowercase
	std::transform(attributeName.begin(), attributeName.end(), attributeName.begin(), ::tolower);
	if (     attributeName == "name") {
		simulation.name = attributeValue;
	}
	else if (attributeName == "description") {
		simulation.description = attributeValue;
	}
	else {
		_stream << "Unknown attribute: '" << attribute->Name() << "' at row: " << attribute->Row() << ", col: " << attribute->Column();
		Error::_errMsg = _stream.str();
		Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
		return 1;
	}

	return 0;
}

int XmlFileAdapter::DeserializeSettings(TiXmlElement *xmlElement, Settings *settings)
{
	// Iterate over the attributes
	for (TiXmlAttribute *attribute = xmlElement->FirstAttribute(); attribute; attribute = attribute->Next() ) {
		if (DeserializeSettingsAttributes(attribute, settings) == 1)
			return 1;
	}

	TiXmlNode *child = xmlElement->FirstChild("Output");
	Output output;
	// If the Output tag is defined than overwrite the default values
	if (child != 0) {
		if (child->ToElement() == 0) {
			_stream << "Invalid xml element at row: " << xmlElement->Row() << ", col: " << xmlElement->Column();
			Error::_errMsg = _stream.str();
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			return 1;
		}
		if (DeserializeOutput(child->ToElement(), &output) == 1) {
			return 1;
		}
	}
    settings->output = new Output(output);

	child = xmlElement->FirstChild("Integrator");
	Integrator integrator;
	// If the Integrator tag is defined than overwrite the default values
	if (child != 0) {
		if (child->ToElement() == 0) {
			_stream << "Invalid xml element at row: " << xmlElement->Row() << ", col: " << xmlElement->Column();
			Error::_errMsg = _stream.str();
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			return 1;
		}
		if (DeserializeIntegrator(child->ToElement(), &integrator) == 1) {
			return 1;
		}
	}
    settings->integrator = new Integrator(integrator);

	child = xmlElement->FirstChild("TimeLine");
	if (child == 0) {
		Error::_errMsg = "Missing TimeLine tag";
		Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
		return 1;
	}
	if (child->ToElement() == 0) {
		_stream << "Invalid xml element at row: " << xmlElement->Row() << ", col: " << xmlElement->Column();
		Error::_errMsg = _stream.str();
		Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
		return 1;
	}
	TimeLine timeLine;
	if (DeserializeTimeLine(child->ToElement(), &timeLine) == 1) {
		Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
		return 1;
	}
    settings->timeLine = new TimeLine(timeLine);

	double value = 0.0;
	std::string unit;
	child = xmlElement->FirstChild("Ejection");
	if (child != 0) {
		if (child->ToElement() == 0) {
			_stream << "Invalid xml element at row: " << xmlElement->Row() << ", col: " << xmlElement->Column();
			Error::_errMsg = _stream.str();
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			return 1;
		}
		if (DeserializeDimensionalQuantity(child->ToElement(), &value, &unit) == 1) {
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			return 1;
		}
		if (!ValidDistanceUnitAttribute(child->ToElement(), unit, "unit")) {
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			return 1;
		}
		// unit conversion
		UnitTool::DistanceToAu(unit, value);
		settings->ejection = value;
	}

	child = xmlElement->FirstChild("HitCentrum");
	if (child != 0) {
		if (child->ToElement() == 0) {
			_stream << "Invalid xml element at row: " << xmlElement->Row() << ", col: " << xmlElement->Column();
			Error::_errMsg = _stream.str();
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			return 1;
		}
		if (DeserializeDimensionalQuantity(child->ToElement(), &value, &unit) == 1) {
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			return 1;
		}
		if (!ValidDistanceUnitAttribute(child->ToElement(), unit, "unit")) {
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			return 1;
		}
		// unit conversion
		UnitTool::DistanceToAu(unit, value);
		settings->hitCentrum = value;
	}

	EventCondition e;
	child = xmlElement->FirstChild("Collision");
	if (child != 0) {
		if (child->ToElement() == 0) {
			_stream << "Invalid xml element at row: " << xmlElement->Row() << ", col: " << xmlElement->Column();
			Error::_errMsg = _stream.str();
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			return 1;
		}
		if (DeserializeEventCondition(child->ToElement(), &e) == 1) {
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			return 1;
		}
		settings->collision = new EventCondition(e);
	}

	e = EventCondition();
	child = xmlElement->FirstChild("CloseEncounter");
	if (child != 0) {
		if (child->ToElement() == 0) {
			_stream << "Invalid xml element at row: " << xmlElement->Row() << ", col: " << xmlElement->Column();
			Error::_errMsg = _stream.str();
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			return 1;
		}
		if (DeserializeEventCondition(child->ToElement(), &e) == 1) {
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			return 1;
		}
		settings->closeEncounter = new EventCondition(e);
	}

	e = EventCondition();
	child = xmlElement->FirstChild("WeakCapture");
	if (child != 0) {
		if (child->ToElement() == 0) {
			_stream << "Invalid xml element at row: " << xmlElement->Row() << ", col: " << xmlElement->Column();
			Error::_errMsg = _stream.str();
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			return 1;
		}
		if (DeserializeEventCondition(child->ToElement(), &e) == 1) {
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			return 1;
		}
		settings->weakCapture = new EventCondition(e);
	}

	return 0;
}

int XmlFileAdapter::DeserializeSettingsAttributes(TiXmlAttribute *attribute, Settings *settings)
{
	std::string attributeName = attribute->Name();
	std::string attributeValue = attribute->Value();
	Tools::Trim(attributeValue);
	// Transform the name to lowercase
	std::transform(attributeName.begin(), attributeName.end(), attributeName.begin(), ::tolower);
	std::transform(attributeValue.begin(), attributeValue.end(), attributeValue.begin(), ::tolower);
	if (     attributeName == "enabledistinctstarttimes") {
		if (	attributeValue == "true")  {
			settings->enableDistinctStartTimes = true;
		}
		else if(attributeValue == "false") {
			settings->enableDistinctStartTimes = false;
		}
		else {
			_stream << "Invalid value: " << attribute->Name() << "=" << attribute->Value() << " at row: " << attribute->Row() << ", col: " << attribute->Column();
			Error::_errMsg = _stream.str();
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			return 1;
		}
	}
	else if (attributeName == "barycentric") {
		if (	attributeValue == "true")  {
			settings->baryCentric = true;
		}
		else if(attributeValue == "false") {
			settings->baryCentric = false;
		}
		else {
			_stream << "Invalid value: " << attribute->Name() << "=" << attribute->Value() << " at row: " << attribute->Row() << ", col: " << attribute->Column();
			Error::_errMsg = _stream.str();
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			return 1;
		}
	}
	else {
		_stream << "Unknown attribute: '" << attribute->Name() << "' at row: " << attribute->Row() << ", col: " << attribute->Column();
		Error::_errMsg = _stream.str();
		Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
		return 1;
	}

	return 0;
}

int XmlFileAdapter::DeserializeOutput(TiXmlElement *xmlElement, Output *output)
{
	std::string text;

	TiXmlNode *child = xmlElement->FirstChild("Phases");
	if (child != 0) {
		if (child->ToElement() == 0) {
			_stream << "Invalid xml element at row: " << xmlElement->Row() << ", col: " << xmlElement->Column();
			Error::_errMsg = _stream.str();
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			return 1;
		}
		text = child->ToElement()->GetText();
		Tools::Trim(text);
		if (text.length() > 0) output->phases = text;
	}

	child = xmlElement->FirstChild("ConstantProperties");
	if (child != 0) {
		if (child->ToElement() == 0) {
			_stream << "Invalid xml element at row: " << xmlElement->Row() << ", col: " << xmlElement->Column();
			Error::_errMsg = _stream.str();
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			return 1;
		}
		text = child->ToElement()->GetText();
		Tools::Trim(text);
		if (text.length() > 0) output->constantProperties = text;
	}

	child = xmlElement->FirstChild("Integrals");
	if (child != 0) {
		if (child->ToElement() == 0) {
			_stream << "Invalid xml element at row: " << xmlElement->Row() << ", col: " << xmlElement->Column();
			Error::_errMsg = _stream.str();
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			return 1;
		}
		text = child->ToElement()->GetText();
		Tools::Trim(text);
		if (text.length() > 0) output->integrals = text;
	}

	child = xmlElement->FirstChild("VariableProperties");
	if (child != 0) {
		if (child->ToElement() == 0) {
			_stream << "Invalid xml element at row: " << xmlElement->Row() << ", col: " << xmlElement->Column();
			Error::_errMsg = _stream.str();
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			return 1;
		}
		text = child->ToElement()->GetText();
		Tools::Trim(text);
		if (text.length() > 0) output->variableProperties = text;
	}

	child = xmlElement->FirstChild("CompositionProperties");
	if (child != 0) {
		if (child->ToElement() == 0) {
			_stream << "Invalid xml element at row: " << xmlElement->Row() << ", col: " << xmlElement->Column();
			Error::_errMsg = _stream.str();
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			return 1;
		}
		text = child->ToElement()->GetText();
		Tools::Trim(text);
		if (text.length() > 0) output->compositionProperties = text;
	}

	child = xmlElement->FirstChild("TwoBodyAffair");
	if (child != 0) {
		if (child->ToElement() == 0) {
			_stream << "Invalid xml element at row: " << xmlElement->Row() << ", col: " << xmlElement->Column();
			Error::_errMsg = _stream.str();
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			return 1;
		}
		text = child->ToElement()->GetText();
		Tools::Trim(text);
		if (text.length() > 0) output->twoBodyAffair = text;
	}

	child = xmlElement->FirstChild("Log");
	if (child != 0) {
		if (child->ToElement() == 0) {
			_stream << "Invalid xml element at row: " << xmlElement->Row() << ", col: " << xmlElement->Column();
			Error::_errMsg = _stream.str();
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			return 1;
		}
		text = child->ToElement()->GetText();
		Tools::Trim(text);
		if (text.length() > 0) output->log = text;
	}

	return 0;
}

int XmlFileAdapter::DeserializeIntegrator(TiXmlElement *xmlElement, Integrator *integrator)
{
	// Iterate over the attributes
	for (TiXmlAttribute *attribute = xmlElement->FirstAttribute(); attribute; attribute = attribute->Next() ) {
		if (DeserializeIntegratorAttributes(attribute, integrator) == 1) {
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			return 1;
		}
	}

	TiXmlNode *child = xmlElement->FirstChild("Accuracy");
	if (child != 0 ) {
		if (child->ToElement() == 0) {
			_stream << "Invalid xml element at row: " << xmlElement->Row() << ", col: " << xmlElement->Column();
			Error::_errMsg = _stream.str();
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			return 1;
		}
		double accuracy = 0.0;
		if (DeserializeAccuracy(child->ToElement(), &accuracy) == 1) {
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			return 1;
		}
		integrator->accuracy = accuracy;
	}

	return 0;
}

int XmlFileAdapter::DeserializeIntegratorAttributes(TiXmlAttribute *attribute, Integrator *integrator)
{
	std::string attributeName = attribute->Name();
	std::string attributeValue = attribute->Value();
	Tools::Trim(attributeValue);
	// Transform the name and value to lowercase
	std::transform(attributeName.begin(),  attributeName.end(),  attributeName.begin(),  ::tolower);
	std::transform(attributeValue.begin(), attributeValue.end(), attributeValue.begin(), ::tolower);
	if (     attributeName == "name" || attributeName == "xsi:type") {
		integrator->name = attributeValue;
	}
	else if (attributeName == "xmlns:xsi" || attributeName == "xmlns:xsd") {
		;
	}
	else {
		_stream << "Unknown attribute: '" << attribute->Name() << "' at row: " << attribute->Row() << ", col: " << attribute->Column();
		Error::_errMsg = _stream.str();
		Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
		return 1;
	}

	return 0;
}

int XmlFileAdapter::DeserializeAccuracy(TiXmlElement *xmlElement, double *accuracy)
{
	// Iterate over the attributes
	for (TiXmlAttribute *attribute = xmlElement->FirstAttribute(); attribute; attribute = attribute->Next() ) {
		if (DeserializeAccuracyAttributes(attribute, accuracy) == 1) {
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			return 1;
		}
	}

	return 0;
}

int XmlFileAdapter::DeserializeAccuracyAttributes(TiXmlAttribute *attribute, double *accuracy)
{
	std::string attributeName = attribute->Name();
	std::string attributeValue = attribute->Value();
	Tools::Trim(attributeValue);
	// Transform the name to lowercase
	std::transform(attributeName.begin(), attributeName.end(), attributeName.begin(), ::tolower);
	if (     attributeName == "value") {
		if (attribute->QueryDoubleValue(accuracy) != TIXML_SUCCESS) {
			_stream << "Invalid value: " << attribute->Name() << "=" << attribute->Value() << " at row: " << attribute->Row() << ", col: " << attribute->Column();
			Error::_errMsg = _stream.str();
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			return 1;
		}
		if (!Validator::ElementOfAndContainsEndPoints(-16.0, 0.0, *accuracy)) {
			_stream << "Value out of range: " << attribute->Name() << "=" << attribute->Value() << " at row: " << attribute->Row() << ", col: " << attribute->Column();
			Error::_errMsg = _stream.str();
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			return 1;
		}
	}
	else {
		_stream << "Unknown attribute: '" << attribute->Name() << "' at row: " << attribute->Row() << ", col: " << attribute->Column();
		Error::_errMsg = _stream.str();
		Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
		return 1;
	}

	return 0;
}

int XmlFileAdapter::DeserializeTimeLine(TiXmlElement *xmlElement, TimeLine *timeLine)
{
	std::string unitString;
	if (!ValidTimeUnitAttribute(xmlElement, unitString, "unit")) {
		Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
		return 1;
	}
	double length = 0.0;
	double output = 0.0;
	double start = 0.0;
	// Iterate over the attributes
	for (TiXmlAttribute *attribute = xmlElement->FirstAttribute(); attribute; attribute = attribute->Next() ) {
		std::string attributeName = attribute->Name();
		std::string attributeValue = attribute->Value();
		Tools::Trim(attributeValue);
		// Transform the name to lowercase
		std::transform(attributeName.begin(), attributeName.end(), attributeName.begin(), ::tolower);
		if (     attributeName == "start") {
			if (attribute->QueryDoubleValue(&start) != TIXML_SUCCESS) {
				_stream << "Invalid value: " << attribute->Name() << "=" << attribute->Value() << " at row: " << attribute->Row() << ", col: " << attribute->Column();
				Error::_errMsg = _stream.str();
				Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
				return 1;
			}
			// Convert to day
			UnitTool::TimeToDay(unitString, start);
			timeLine->start = start;
			timeLine->startTimeDefined = true;
		}
		else if (attributeName == "length") {
			if (attribute->QueryDoubleValue(&length) != TIXML_SUCCESS) {
				_stream << "Invalid value: " << attribute->Name() << "=" << attribute->Value() << " at row: " << attribute->Row() << ", col: " << attribute->Column();
				Error::_errMsg = _stream.str();
				Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
				return 1;
			}
			UnitTool::TimeToDay(unitString, length);
			timeLine->length = length;
		}
		else if (attributeName == "output") {
			if (attribute->QueryDoubleValue(&output) != TIXML_SUCCESS) {
				_stream << "Invalid value: " << attribute->Name() << "=" << attribute->Value() << " at row: " << attribute->Row() << ", col: " << attribute->Column();
				Error::_errMsg = _stream.str();
				Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
				return 1;
			}
			if (!Validator::GreaterThan(0.0, output)) {
				_stream << "Value out of range: " << attribute->Name() << "=" << attribute->Value() << " at row: " << attribute->Row() << ", col: " << attribute->Column();
				Error::_errMsg = _stream.str();
				Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
				return 1;
			}
			UnitTool::TimeToDay(unitString, output);
			timeLine->output = output;
		}
		else if (attributeName == "unit") {
			;
		}
		else {
			_stream << "Unknown attribute: '" << attribute->Name() << "' at row: " << attribute->Row() << ", col: " << attribute->Column();
			Error::_errMsg = _stream.str();
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			return 1;
		}
	}

	return 0;
}

int XmlFileAdapter::DeserializeEventCondition(TiXmlElement *xmlElement, EventCondition *eventCondition)
{
	// Iterate over the attributes
	for (TiXmlAttribute *attribute = xmlElement->FirstAttribute(); attribute; attribute = attribute->Next() ) {
		if (DeserializeEventConditionAttributes(attribute, eventCondition) == 1) {
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			return 1;
		}
	}

	return 0;
}

int XmlFileAdapter::DeserializeEventConditionAttributes(TiXmlAttribute *attribute, EventCondition *eventCondition)
{
	std::string attributeName = attribute->Name();
	std::string attributeValue = attribute->Value();
	Tools::Trim(attributeValue);
	// Transform the name and value to lowercase
	std::transform(attributeName.begin(), attributeName.end(), attributeName.begin(), ::tolower);
	std::transform(attributeValue.begin(), attributeValue.end(), attributeValue.begin(), ::tolower);
	if (     attributeName == "factor") {
		double factor = 0.0;
		if (attribute->QueryDoubleValue(&factor) != TIXML_SUCCESS) {
			_stream << "Invalid value: " << attribute->Name() << "=" << attribute->Value() << " at row: " << attribute->Row() << ", col: " << attribute->Column();
			Error::_errMsg = _stream.str();
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			return 1;
		}
		if (!Validator::GreaterThanOrEqualTo(0.0, factor)) {
			_stream << "Value out of range: " << attribute->Name() << "=" << attribute->Value() << " at row: " << attribute->Row() << ", col: " << attribute->Column();
			Error::_errMsg = _stream.str();
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			return 1;
		}
		eventCondition->factor = factor;
	}
	else if (attributeName == "stop") {
		if (	attributeValue == "true")  {
			eventCondition->stop = true;
		}
		else if(attributeValue == "false") {
			eventCondition->stop = false;
		}
		else {
			_stream << "Invalid value: " << attribute->Name() << "=" << attribute->Value() << " at row: " << attribute->Row() << ", col: " << attribute->Column();
			Error::_errMsg = _stream.str();
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			return 1;
		}
	}
	else {
		_stream << "Unknown attribute: '" << attribute->Name() << "' at row: " << attribute->Row() << ", col: " << attribute->Column();
		Error::_errMsg = _stream.str();
		Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
		return 1;
	}

	return 0;
}

int XmlFileAdapter::DeserializeBodyGroupList(TiXmlElement *xmlElement, BodyGroupList *list)
{
	if (xmlElement->FirstChild("BodyGroup") == 0) {
		Error::_errMsg = "Missing BodyGroup tag";
		Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
		return 1;
	}
	// Iterate over the BodyGroup tags
	for (TiXmlNode *child = xmlElement->FirstChild("BodyGroup"); child; child = child->NextSibling() ) {
		BodyGroup bodyGroup;
		if (child->ToElement() == 0) {
			_stream << "Invalid xml element at row: " << xmlElement->Row() << ", col: " << xmlElement->Column();
			Error::_errMsg = _stream.str();
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			return 1;
		}
		if (DeserializeBodyGroup(child->ToElement(), &bodyGroup) == 1) {
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			return 1;
		}
		list->items.push_back(bodyGroup);
	}

	return 0;
}

int XmlFileAdapter::DeserializeBodyGroup(TiXmlElement *xmlElement, BodyGroup *bodyGroup)
{
	// Iterate over the attributes
	for (TiXmlAttribute *attribute = xmlElement->FirstAttribute(); attribute; attribute = attribute->Next() ) {
		if (DeserializeBodyGroupAttributes(attribute, bodyGroup) == 1)
			return 1;
	}

	TiXmlNode *child = xmlElement->FirstChild("Items");
	if (child == 0) {
		Error::_errMsg = "Missing Items tag";
		Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
		return 1;
	}

	if (child->FirstChild("Body") == 0)  {
		Error::_errMsg = "Missing Body tag";
		Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
		return 1;
	}
	// Iterate over the Body tags
	for (child = child->FirstChild("Body"); child; child = child->NextSibling() ) {
		Body body;
		if (child->ToElement() == 0) {
			_stream << "Invalid xml element at row: " << xmlElement->Row() << ", col: " << xmlElement->Column();
			Error::_errMsg = _stream.str();
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			return 1;
		}
		if (DeserializeBody(child->ToElement(), &body) == 1) {
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			return 1;
		}
		bodyGroup->items.push_back(body);
	}

	return 0;
}

int XmlFileAdapter::DeserializeBodyGroupAttributes(TiXmlAttribute *attribute, BodyGroup *bodyGroup)
{
	std::string attributeName = attribute->Name();
	std::string attributeValue = attribute->Value();
	Tools::Trim(attributeValue);
	// Transform the name to lowercase
	std::transform(attributeName.begin(), attributeName.end(), attributeName.begin(), ::tolower);
	if (     attributeName == "description") {
		bodyGroup->description = attributeValue;
	}
	else if (attributeName == "epoch") {
		bodyGroup->epoch = attributeValue;
	}
	else if (attributeName == "offset") {
		double offset = 0.0;
		if (attribute->QueryDoubleValue(&offset) != TIXML_SUCCESS) {
			_stream << "Invalid value: " << attribute->Name() << "=" << attribute->Value() << " at row: " << attribute->Row() << ", col: " << attribute->Column();
			Error::_errMsg = _stream.str();
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			return 1;
		}
		bodyGroup->offset = offset;
	}
	else if (attributeName == "referenceframe") {
		bodyGroup->referenceFrame = attributeValue;
	}
	else if (attributeName == "guid") {
		bodyGroup->guid = attributeValue;
	}
	else {
		_stream << "Unknown attribute: '" << attribute->Name() << "' at row: " << attribute->Row() << ", col: " << attribute->Column();
		Error::_errMsg = _stream.str();
		Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
		return 1;
	}

	return 0;
}

int XmlFileAdapter::DeserializeBody(TiXmlElement *xmlElement, Body *body)
{
	// Iterate over the attributes
	for (TiXmlAttribute *attribute = xmlElement->FirstAttribute(); attribute; attribute = attribute->Next() ) {
		if (DeserializeBodyAttributes(attribute, body) == 1)
			return 1;
	}
	
	TiXmlNode *phaseChild	= xmlElement->FirstChild("Phase");
	TiXmlNode *oeChild		= xmlElement->FirstChild("OrbitalElement");
	if (body->type == CentralBody && phaseChild == 0 && oeChild == 0) {
		body->phase = new Phase(body->GetId());
	}
	else {
		if (phaseChild != 0 && oeChild != 0 ) {
			_stream << "You cannot define both Phase and OrbitalElement tags for the same Body! Row: " << xmlElement->Row() << ", col: " << xmlElement->Column();
			Error::_errMsg = _stream.str();
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			return 1;
		}
		else if (phaseChild == 0 && oeChild == 0) {
			_stream << "Either Phase or OrbitalElement tag must be defined for a Body! Row: " << xmlElement->Row() << ", col: " << xmlElement->Column();
			Error::_errMsg = _stream.str();
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			return 1;
		}
		else if (phaseChild != 0) {
			Phase phase;
			if (phaseChild->ToElement() == 0) {
				_stream << "Invalid xml element at row: " << xmlElement->Row() << ", col: " << xmlElement->Column();
				Error::_errMsg = _stream.str();
				Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
				return 1;
			}
			if (DeserializePhase(phaseChild->ToElement(), &phase) == 1) {
				Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
				return 1;
			}
			phase.bodyId = body->GetId();
			body->phase = new Phase(phase);
		}
		else {
			OrbitalElement orbElem;
			if (oeChild->ToElement() == 0) {
				_stream << "Invalid xml element at row: " << xmlElement->Row() << ", col: " << xmlElement->Column();
				Error::_errMsg = _stream.str();
				Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
				return 1;
			}
			if (DeserializeOrbitalElement(oeChild->ToElement(), &orbElem) == 1) {
				Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
				return 1;
			}
			body->orbitalElement = new OrbitalElement(orbElem);
		}
	}

	if (body->type != TestParticle) {
		TiXmlNode *child = xmlElement->FirstChild("Characteristics");
		if (child == 0) {
			Error::_errMsg = "Missing Characteristics tag";
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			return 1;
		}
		if (child->ToElement() == 0) {
			_stream << "Invalid xml element at row: " << xmlElement->Row() << ", col: " << xmlElement->Column();
			Error::_errMsg = _stream.str();
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			return 1;
		}
		Characteristics characteristics;
		if (DeserializeCharacteristics(child->ToElement(), &characteristics, body->type) == 1) {
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			return 1;
		}
		body->characteristics = new Characteristics(characteristics);

		//DragCoefficient dragCoefficient;
		//child = xmlElement->FirstChild("DragCoefficient");
		//if (child != 0) {
		//	if (child->ToElement() == 0) {
		//		_stream << "Invalid xml element at row: " << xmlElement->Row() << ", col: " << xmlElement->Column();
		//		Error::_errMsg = _stream.str();
		//		Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
		//		return 1;
		//	}
		//	if (DeserializeDragCoefficient(child->ToElement(), &dragCoefficient) == 1) {
		//		Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
		//		return 1;
		//	}
		//	body->characteristics->stokes = dragCoefficient.stokes;
		//}

		child = xmlElement->FirstChild("Migration");
		if (child != 0) {
			if (child->ToElement() == 0) {
				_stream << "Invalid xml element at row: " << xmlElement->Row() << ", col: " << xmlElement->Column();
				Error::_errMsg = _stream.str();
				Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
				return 1;
			}
			if (DeserializeMigrationType(child->ToElement(), body) == 1) {
				Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
				return 1;
			}
		}
	}

	return 0;
}

int XmlFileAdapter::DeserializeBodyAttributes(TiXmlAttribute *attribute, Body *body)
{
	std::string attributeName = attribute->Name();
	std::string attributeValue = attribute->Value();
	Tools::Trim(attributeValue);
	// Transform the name to lowercase
	std::transform(attributeName.begin(), attributeName.end(), attributeName.begin(), ::tolower);
	if (     attributeName == "name") {
		body->name = attributeValue;
	}
	else if (attributeName == "type") {
		std::string bodyType = attributeValue;
		std::transform(bodyType.begin(), bodyType.end(), bodyType.begin(), ::tolower);
		if (SetBodyType(bodyType, body) == 1) {
			_stream << "Unknown type: " << attribute->Name() << "=" << attribute->Value() << " at row: " << attribute->Row() << ", col: " << attribute->Column();
			Error::_errMsg = _stream.str();
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			return 1;
		}
	}
	else if (attributeName == "mpcorbittype") {
		std::string mpcOrbitType = attributeValue;
		std::transform(mpcOrbitType.begin(), mpcOrbitType.end(), mpcOrbitType.begin(), ::tolower);
		if (SetMPCOrbitType(mpcOrbitType, body) == 1) {
			_stream << "Unknown type: " << attribute->Name() << "=" << attribute->Value() << " at row: " << attribute->Row() << ", col: " << attribute->Column();
			Error::_errMsg = _stream.str();
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			return 1;
		}
		body->mPCOrbitTypestr = attributeValue;
	}
	else if (attributeName == "designation" || attributeName == "des") {
		body->designation = attributeValue;
	}
	else if (attributeName == "provisionaldesignation" || attributeName == "provdes") {
		body->provisionalDesignation = attributeValue;
	}
	else if (attributeName == "reference" || attributeName == "ref") {
		body->reference = attributeValue;
	}
	else if (attributeName == "opposition" || attributeName == "opp") {
		body->opposition = attributeValue;
	}
	else if (attributeName == "ln") {
		std::string bodyLn = attributeValue;
		std::transform(bodyLn.begin(), bodyLn.end(), bodyLn.begin(), ::tolower);
		if (SetLn(bodyLn, body)) {
			_stream << "Unknown ln: " << attribute->Name() << "=" << attribute->Value() << " at row: " << attribute->Row() << ", col: " << attribute->Column();
			Error::_errMsg = _stream.str();
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			return 1;
		}
	}
	else if (attributeName == "guid") {
		body->guid = attributeValue;
	}
	else {
		_stream << "Unknown attribute: '" << attribute->Name() << "' at row: " << attribute->Row() << ", col: " << attribute->Column();
		Error::_errMsg = _stream.str();
		Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
		return 1;
	}

	return 0;
}

int XmlFileAdapter::SetBodyType(std::string type, Body *body)
{
	if (     type == "centralbody")			{ body->type = CentralBody;			}
	else if (type == "giantplanet")			{ body->type = GiantPlanet;			}
	else if (type == "rockyplanet")			{ body->type = RockyPlanet;			}
	else if (type == "protoplanet")			{ body->type = ProtoPlanet;			}
	else if (type == "superplanetesimal")	{ body->type = SuperPlanetesimal;	}
	else if (type == "planetesimal")		{ body->type = Planetesimal;		}
	else if (type == "testparticle")		{ body->type = TestParticle;		}
	else {
		return 1;
	}

	return 0;
}

int XmlFileAdapter::SetMPCOrbitType(std::string type, Body *body)
{
	if (     type == "aten")											{ body->mPCOrbitType = Aten;										}
	else if (type == "apollo")											{ body->mPCOrbitType = Apollo;										}
	else if (type == "amor")											{ body->mPCOrbitType = Amor;										}
	else if (type == "objectwithqlt1_665")								{ body->mPCOrbitType = ObjectWithqLt1_665;							}
	else if (type == "hungaria")										{ body->mPCOrbitType = Hungaria;									}
	else if (type == "phocaea")											{ body->mPCOrbitType = Phocaea;										}
	else if (type == "hilda")											{ body->mPCOrbitType = Hilda;										}
	else if (type == "JupiterTrojan")									{ body->mPCOrbitType = JupiterTrojan;								}
	else if (type == "Centaur")											{ body->mPCOrbitType = Centaur;										}
	else if (type == "Plutino")											{ body->mPCOrbitType = Plutino;										}
	else if (type == "OtherResonantTNO")								{ body->mPCOrbitType = OtherResonantTNO;							}
	else if (type == "Cubewano")										{ body->mPCOrbitType = Cubewano;									}
	else if (type == "ScatteredDisk")									{ body->mPCOrbitType = ScatteredDisk;								}
	else if (type == "ObjectIsNEO")										{ body->mPCOrbitType = ObjectIsNEO;									}
	else if (type == "ObjectIs1kmOrLargerNEO")							{ body->mPCOrbitType = ObjectIs1kmOrLargerNEO;						}
	else if (type == "OneOppositionObjectSeenAtEarlierOpposition")		{ body->mPCOrbitType = OneOppositionObjectSeenAtEarlierOpposition;	}
	else if (type == "CriticalListNumberedObject")						{ body->mPCOrbitType = CriticalListNumberedObject;					}
	else if (type == "ObjectIsPHA")										{ body->mPCOrbitType = ObjectIsPHA;									}
	else {
		return 1;
	}

	return 0;
}

int XmlFileAdapter::SetLn(std::string ln, Body *body)
{
	if (     ln == "l1") { body->ln = L1; }
	else if (ln == "l2") { body->ln = L2; }
	else if (ln == "l3") { body->ln = L3; }
	else if (ln == "l4") { body->ln = L4; }
	else if (ln == "l5") { body->ln = L5; }
	else {
		return 1;
	}

	return 0; 
}

int XmlFileAdapter::DeserializeMigration(TiXmlElement *xmlElement, Body *body)
{
	// Iterate over the attributes
	for (TiXmlAttribute *attribute = xmlElement->FirstAttribute(); attribute; attribute = attribute->Next() ) {
		if (DeserializeMigrationAttributes(attribute, body) == 1) {
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			return 1;
		}
	}

	return 0;
}

int XmlFileAdapter::DeserializeMigrationAttributes(TiXmlAttribute *attribute, Body *body)
{
	double x = 0.0;

	std::string attributeName = attribute->Name();
	std::string attributeValue = attribute->Value();
	Tools::Trim(attributeValue);
	// Transform the name to lowercase
	std::transform(attributeName.begin(), attributeName.end(), attributeName.begin(), ::tolower);
	if (     attributeName == "type") {
		std::string migType = attributeValue;
		std::transform(migType.begin(), migType.end(), migType.begin(), ::tolower);
		if (SetMigrationType(migType, body) == 1) {
			_stream << "Unknown type: " << attribute->Name() << "=" << attribute->Value() << " at row: " << attribute->Row() << ", col: " << attribute->Column();
			Error::_errMsg = _stream.str();
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			return 1;
		}
	}
	else if (attributeName == "stopat") {
		if (attribute->QueryDoubleValue(&x) != TIXML_SUCCESS) {
			_stream << "Invalid value: " << attribute->Name() << "=" << attribute->Value() << " at row: " << attribute->Row() << ", col: " << attribute->Column();
			Error::_errMsg = _stream.str();
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			return 1;
		}
		if (!Validator::GreaterThan(0.0, x)) {
			_stream << "Value out of range: " << attribute->Name() << "=" << attribute->Value() << " at row: " << attribute->Row() << ", col: " << attribute->Column();
			Error::_errMsg = _stream.str();
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			return 1;
		}
		body->migrationStopAt = x;
	}
	else {
		_stream << "Unknown attribute: '" << attribute->Name() << "' at row: " << attribute->Row() << ", col: " << attribute->Column();
		Error::_errMsg = _stream.str();
		Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
		return 1;
	}

	return 0;
}

int XmlFileAdapter::SetMigrationType(std::string migrationType, Body *body)
{
	if (     migrationType == "i")  { body->migrationType = TypeI;  }
	else if (migrationType == "ii") { body->migrationType = TypeII; }
	else {
		return 1;
	}

	return 0;
}

int XmlFileAdapter::DeserializePhase(TiXmlElement *xmlElement, Phase *phase) 
{
	TiXmlNode *child = xmlElement->FirstChild("Position");
	if (child == 0) {
		Error::_errMsg = "Missing Position tag";
		Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
		return 1;
	}
	if (child->ToElement() == 0) {
		_stream << "Invalid xml element at row: " << xmlElement->Row() << ", col: " << xmlElement->Column();
		Error::_errMsg = _stream.str();
		Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
		return 1;
	}
	if (DeserializePosition(child->ToElement(), phase->position) == 1) {
		Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
		return 1;
	}

	child = xmlElement->FirstChild("Velocity");
	if (child == 0) {
		Error::_errMsg = "Missing Velocity tag";
		Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
		return 1;
	}
	if (child->ToElement() == 0) {
		_stream << "Invalid xml element at row: " << xmlElement->Row() << ", col: " << xmlElement->Column();
		Error::_errMsg = _stream.str();
		Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
		return 1;
	}
	if (DeserializeVelocity(child->ToElement(), phase->velocity) == 1) {
		Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
		return 1;
	}
	
	return 0;
}

int XmlFileAdapter::DeserializePosition(TiXmlElement *xmlElement, Vector& position)
{
	std::string unitString;
	if (!ValidDistanceUnitAttribute(xmlElement, unitString, "unit")) {
		Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
		return 1;
	}

	if (DeserializeVector(xmlElement, position) == 1) {
		Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
		return 1;
	}
	// length unit conversion
	UnitTool::DistanceToAu(unitString, position.x);
	UnitTool::DistanceToAu(unitString, position.y);
	UnitTool::DistanceToAu(unitString, position.z);

	return 0;
}

int XmlFileAdapter::DeserializeVelocity(TiXmlElement *xmlElement, Vector& velocity)
{
	std::string unitString;
	if (!ValidVelocityUnitAttribute(xmlElement, unitString, "unit")) {
		Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
		return 1;
	}
	if (DeserializeVector(xmlElement, velocity) == 1) {
		Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
		return 1;
	}
	// velocity unit conversion
	UnitTool::VelocityToCM(unitString, velocity.x);
	UnitTool::VelocityToCM(unitString, velocity.y);
	UnitTool::VelocityToCM(unitString, velocity.z);

	return 0;
}

int XmlFileAdapter::DeserializeVector(TiXmlElement *xmlElement, Vector& vector)
{
	double x, y, z;
	x = y = z = 0.0;
	// Iterate over the attributes
	for (TiXmlAttribute *attribute = xmlElement->FirstAttribute(); attribute; attribute = attribute->Next() ) {
		std::string attributeName = attribute->Name();
		std::string attributeValue = attribute->Value();
		Tools::Trim(attributeValue);
		// Transform the name to lowercase
		std::transform(attributeName.begin(), attributeName.end(), attributeName.begin(), ::tolower);

		if (     attributeName == "x") {
			if (attribute->QueryDoubleValue(&x) != TIXML_SUCCESS) {
				_stream << "Invalid value: " << attribute->Name() << "=" << attribute->Value() << " at row: " << attribute->Row() << ", col: " << attribute->Column();
				Error::_errMsg = _stream.str();
				Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
				return 1;
			}
		}
		else if (attributeName == "y") {
			if (attribute->QueryDoubleValue(&y) != TIXML_SUCCESS) {
				_stream << "Invalid value: " << attribute->Name() << "=" << attribute->Value() << " at row: " << attribute->Row() << ", col: " << attribute->Column();
				Error::_errMsg = _stream.str();
				Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
				return 1;
			}
		}
		else if (attributeName == "z") {
			if (attribute->QueryDoubleValue(&z) != TIXML_SUCCESS) {
				_stream << "Invalid value: " << attribute->Name() << "=" << attribute->Value() << " at row: " << attribute->Row() << ", col: " << attribute->Column();
				Error::_errMsg = _stream.str();
				Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
				return 1;
			}
		}
		else if (attributeName == "unit") {
			;
		}
		else {
			_stream << "Unknown attribute: '" << attribute->Name() << "' at row: " << attribute->Row() << ", col: " << attribute->Column();
			Error::_errMsg = _stream.str();
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			return 1;
		}
	}

	vector.x = x;
	vector.y = y;
	vector.z = z;

	return 0;
}

int XmlFileAdapter::DeserializeOrbitalElement(TiXmlElement *xmlElement, OrbitalElement *oe)
{
	std::string distanceUnitString;
	std::string angleUnitString;

	if (!ValidDistanceUnitAttribute(xmlElement, distanceUnitString, "distanceunit")) {
		Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
		return 1;
	}
	if (!ValidAngleUnitAttribute(xmlElement, angleUnitString, "angleunit")) {
		Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
		return 1;
	}

	double sma, ecc, inc, peri, node, mean;
	sma = ecc = inc = peri = node = mean = 0.0;
	// Iterate over the attributes
	for (TiXmlAttribute *attribute = xmlElement->FirstAttribute(); attribute; attribute = attribute->Next() ) {
		std::string attributeName = attribute->Name();
		std::string attributeValue = attribute->Value();
		// Transform the name to lowercase
		std::transform(attributeName.begin(), attributeName.end(), attributeName.begin(), ::tolower);
		if (     attributeName == "semimajoraxis" || attributeName == "a") {
			if (attribute->QueryDoubleValue(&sma) != TIXML_SUCCESS) {
				_stream << "Invalid value: " << attribute->Name() << "=" << attribute->Value() << " at row: " << attribute->Row() << ", col: " << attribute->Column();
				Error::_errMsg = _stream.str();
				Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
				return 1;
			}
			if (!Validator::GreaterThan(0, sma)) {
				_stream << "Value out of range: " << attribute->Name() << "=" << attribute->Value() << " at row: " << attribute->Row() << ", col: " << attribute->Column();
				Error::_errMsg = _stream.str();
				Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
				return 1;
			}
			UnitTool::DistanceToAu(distanceUnitString, sma);
		}
		else if (attributeName == "eccentricity" || attributeName == "e") {
			if (attribute->QueryDoubleValue(&ecc) != TIXML_SUCCESS) {
				_stream << "Invalid value: " << attribute->Name() << "=" << attribute->Value() << " at row: " << attribute->Row() << ", col: " << attribute->Column();
				Error::_errMsg = _stream.str();
				Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
				return 1;
			}
			if (!Validator::ElementOfAndContainsLower(0.0, 1.0, ecc)) {
				_stream << "Value out of range: " << attribute->Name() << "=" << attribute->Value() << " at row: " << attribute->Row() << ", col: " << attribute->Column();
				Error::_errMsg = _stream.str();
				Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
				return 1;
			}
		}
		else if (attributeName == "inclination" || attributeName == "inc") {
			if (attribute->QueryDoubleValue(&inc) != TIXML_SUCCESS) {
				_stream << "Invalid value: " << attribute->Name() << "=" << attribute->Value() << " at row: " << attribute->Row() << ", col: " << attribute->Column();
				Error::_errMsg = _stream.str();
				Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
				return 1;
			}
			UnitTool::AngleToRadian(angleUnitString, inc);
			if (!Validator::ElementOfAndContainsEndPoints(0, Constants::Pi, inc)) {
				_stream << "Value out of range: " << attribute->Name() << "=" << attribute->Value() << " at row: " << attribute->Row() << ", col: " << attribute->Column();
				Error::_errMsg = _stream.str();
				Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
				return 1;
			}
		}
		else if (attributeName == "argumentofpericenter" || attributeName == "peri") {
			if (attribute->QueryDoubleValue(&peri) != TIXML_SUCCESS) {
				_stream << "Invalid value: " << attribute->Name() << "=" << attribute->Value() << " at row: " << attribute->Row() << ", col: " << attribute->Column();
				Error::_errMsg = _stream.str();
				Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
				return 1;
			}
			UnitTool::AngleToRadian(angleUnitString, peri);
			Ephemeris::ShiftIntoRange(0.0, 2.0*Constants::Pi, peri);
		}
		else if (attributeName == "longitudeofnode" || attributeName == "node") {
			if (attribute->QueryDoubleValue(&node) != TIXML_SUCCESS) {
				_stream << "Invalid value: " << attribute->Name() << "=" << attribute->Value() << " at row: " << attribute->Row() << ", col: " << attribute->Column();
				Error::_errMsg = _stream.str();
				Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
				return 1;
			}
			UnitTool::AngleToRadian(angleUnitString, node);
			Ephemeris::ShiftIntoRange(0.0, 2.0*Constants::Pi, node);
		}
		else if (attributeName == "meananomaly" || attributeName == "m") {
			if (attribute->QueryDoubleValue(&mean) != TIXML_SUCCESS) {
				_stream << "Invalid value: " << attribute->Name() << "=" << attribute->Value() << " at row: " << attribute->Row() << ", col: " << attribute->Column();
				Error::_errMsg = _stream.str();
				Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
				return 1;
			}
			UnitTool::AngleToRadian(angleUnitString, mean);
			Ephemeris::ShiftIntoRange(0.0, 2.0*Constants::Pi, mean);
		}
		else if (attributeName == "angleunit" || attributeName == "distanceunit") {
			;
		}
		else {
			_stream << "Unknown attribute: '" << attribute->Name() << "' at row: " << attribute->Row() << ", col: " << attribute->Column();
			Error::_errMsg = _stream.str();
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			return 1;
		}
	}

	oe->semiMajorAxis			= sma;
	oe->eccentricity			= ecc;
	oe->inclination				= inc;
	oe->argumentOfPericenter	= peri;
	oe->longitudeOfNode			= node;
	oe->meanAnomaly				= mean;

	return 0;
}

// The BodyType is needed since in the case of super-planetesimals both the radius and the density can be separately defined
int XmlFileAdapter::DeserializeCharacteristics(TiXmlElement *xmlElement, Characteristics *characteristics, BodyType bodyType)
{
	for (TiXmlAttribute *attribute = xmlElement->FirstAttribute(); attribute; attribute = attribute->Next() ) {
		if (DeserializeCharacteristicsAttributes(attribute, characteristics) == 1) {
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			return 1;
		}
	}

	double value = 0.0;
	std::string unit;

	TiXmlNode *child = xmlElement->FirstChild("Mass");
	if (child == 0) {
		Error::_errMsg = "Missing Mass tag";
		Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
		return 1;
	}
	else {
		if (child->ToElement() == 0) {
			_stream << "Invalid xml element at row: " << xmlElement->Row() << ", col: " << xmlElement->Column();
			Error::_errMsg = _stream.str();
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			return 1;
		}
		if (DeserializeDimensionalQuantity(child->ToElement(), &value, &unit) == 1) {
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			return 1;
		}
		if (!ValidMassUnitAttribute(child->ToElement(), unit, "unit")) {
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			return 1;
		}
		// unit conversion
		UnitTool::MassToSolar(unit, value);
		characteristics->mass = value;
	}

	child = xmlElement->FirstChild("Radius");
	if (child != 0) {
		if (child->ToElement() == 0) {
			_stream << "Invalid xml element at row: " << xmlElement->Row() << ", col: " << xmlElement->Column();
			Error::_errMsg = _stream.str();
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			return 1;
		}
		if (DeserializeDimensionalQuantity(child->ToElement(), &value, &unit) == 1) {
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			return 1;
		}
		if (!ValidDistanceUnitAttribute(child->ToElement(), unit, "unit")) {
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			return 1;
		}
		// unit conversion
		UnitTool::DistanceToAu(unit, value);
		characteristics->radius = value;
	}

	child = xmlElement->FirstChild("Density");
	// Both the radius and the density cannot be defined simultaneously for bodies other than super-planetesimals
	if (bodyType != SuperPlanetesimal && child != 0 && characteristics->radius > 0) {
		_stream << "The radius and density of a body cannot be defined simultaneously! Row: " << child->Row() << ", col: " << child->Column();
		Error::_errMsg = _stream.str();
		Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
		return 1;
	}
	if (child != 0) {
		if (child->ToElement() == 0) {
			_stream << "Invalid xml element at row: " << xmlElement->Row() << ", col: " << xmlElement->Column();
			Error::_errMsg = _stream.str();
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			return 1;
		}
		if (DeserializeDimensionalQuantity(child->ToElement(), &value, &unit) == 1) {
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			return 1;
		}
		if (!ValidVolumeDensityUnitAttribute(child->ToElement(), unit, "unit")) {
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			return 1;
		}
		// unit conversion
		UnitTool::VolumeDensityToCM(unit, value);
		characteristics->density = value;
	}

	child = xmlElement->FirstChild("ComponentList");
	if (child != 0) {
		if (child->ToElement() == 0) {
			_stream << "Invalid xml element at row: " << xmlElement->Row() << ", col: " << xmlElement->Column();
			Error::_errMsg = _stream.str();
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			return 1;
		}
		if (DeserializeComponentList(child->ToElement(), &(characteristics->componentList)) == 1) {
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			return 1;
		}
		if (!characteristics->componentList.empty()) {
			double sum = 100.0;
			for (std::list<Component>::iterator it = characteristics->componentList.begin(); it != characteristics->componentList.end(); it++) {
					sum -= it->ratio;
			}
			if (fabs(sum) > 1.0e-4 ) {
				_stream << "The sum of the ComponentList tag is not 100%! Row: " << child->Row() << ", col: " << child->Column();
				Error::_errMsg = _stream.str();
				Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			}
		}
		else {
			_stream << "ComponentList tag is empty! Row: " << child->Row() << ", col: " << child->Column();
			Error::_errMsg = _stream.str();
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			return 1;
		}
	}

	return 0;
}

int XmlFileAdapter::DeserializeCharacteristicsAttributes(TiXmlAttribute *attribute, Characteristics *characteristics)
{
	std::string attributeName = attribute->Name();
	std::string attributeValue = attribute->Value();
	// Transform the name to lowercase
	std::transform(attributeName.begin(), attributeName.end(), attributeName.begin(), ::tolower);
	if (     attributeName == "absvismag") {
		double value;
		if (attribute->QueryDoubleValue(&value) != TIXML_SUCCESS) {
			_stream << "Invalid value: " << attribute->Name() << "=" << attribute->Value() << " at row: " << attribute->Row() << ", col: " << attribute->Column();
			Error::_errMsg = _stream.str();
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			return 1;
		}
		characteristics->absVisMag = value;
	}
	else if (attributeName == "cd") {
		double value;
		if (attribute->QueryDoubleValue(&value) != TIXML_SUCCESS) {
			_stream << "Invalid value: " << attribute->Name() << "=" << attribute->Value() << " at row: " << attribute->Row() << ", col: " << attribute->Column();
			Error::_errMsg = _stream.str();
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			return 1;
		}
		if (!Validator::GreaterThanOrEqualTo(0.0, value)) {
			_stream << "Value out of range: " << attribute->Name() << "=" << attribute->Value() << " at row: " << attribute->Row() << ", col: " << attribute->Column();
			Error::_errMsg = _stream.str();
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			return 1;
		}
		characteristics->stokes = cd;
	}
	else {
		_stream << "Unknown attribute: '" << attribute->Name() << "' at row: " << attribute->Row() << ", col: " << attribute->Column();
		Error::_errMsg = _stream.str();
		Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
		return 1;
	}

	return 0;
}

//int XmlFileAdapter::DeserializeDragCoefficient(TiXmlElement *xmlElement, DragCoefficient *dragCoefficient)
//{
//	for (TiXmlAttribute *attribute = xmlElement->FirstAttribute(); attribute; attribute = attribute->Next() ) {
//		if (DeserializeDragCoefficientAttributes(attribute, dragCoefficient) == 1) {
//			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
//			return 1;
//		}
//	}
//
//	return 0;
//}
//
//int XmlFileAdapter::DeserializeDragCoefficientAttributes(TiXmlAttribute *attribute, DragCoefficient *dragCoefficient)
//{
//	std::string attributeName = attribute->Name();
//	std::string attributeValue = attribute->Value();
//	// Transform the name to lowercase
//	std::transform(attributeName.begin(), attributeName.end(), attributeName.begin(), ::tolower);
//	if (     attributeName == "stokes") {
//		double value;
//		if (attribute->QueryDoubleValue(&value) != TIXML_SUCCESS) {
//			_stream << "Invalid value: " << attribute->Name() << "=" << attribute->Value() << " at row: " << attribute->Row() << ", col: " << attribute->Column();
//			Error::_errMsg = _stream.str();
//			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
//			return 1;
//		}
//		if (!Validator::GreaterThanOrEqualTo(0.0, value)) {
//			_stream << "Value out of range: " << attribute->Name() << "=" << attribute->Value() << " at row: " << attribute->Row() << ", col: " << attribute->Column();
//			Error::_errMsg = _stream.str();
//			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
//			return 1;
//		}
//		dragCoefficient->stokes = value;
//	}
//	else {
//		_stream << "Unknown attribute: '" << attribute->Name() << "' at row: " << attribute->Row() << ", col: " << attribute->Column();
//		Error::_errMsg = _stream.str();
//		Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
//		return 1;
//	}
//
//	return 0;
//}

int XmlFileAdapter::DeserializeComponentList(TiXmlElement *xmlElement, std::list<Component> *componentList)
{
	// Iterate over the Component tags
	for (TiXmlNode *child = xmlElement->FirstChild("Component"); child; child = child->NextSibling() ) {
		Component component;
		if (child->ToElement() == 0) {
			_stream << "Invalid xml element at row: " << xmlElement->Row() << ", col: " << xmlElement->Column();
			Error::_errMsg = _stream.str();
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			return 1;
		}
		if (DeserializeComponent(child->ToElement(), &component) == 1) {
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			return 1;
		}
		componentList->push_back(component);
	}

	return 0;
}

int XmlFileAdapter::DeserializeComponent(TiXmlElement *xmlElement, Component *component)
{
	// Iterate over the attributes
	for (TiXmlAttribute *attribute = xmlElement->FirstAttribute(); attribute; attribute = attribute->Next() ) {
		if (DeserializeComponentAttributes(attribute, component) == 1) {
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			return 1;
		}
	}

	return 0;
}

int XmlFileAdapter::DeserializeComponentAttributes(TiXmlAttribute *attribute, Component *component)
{
	std::string attributeName = attribute->Name();
	std::string attributeValue = attribute->Value();
	// Transform the name to lowercase
	std::transform(attributeName.begin(), attributeName.end(), attributeName.begin(), ::tolower);
	if (     attributeName == "name") {
		component->name = attributeValue;
	}
	else if (attributeName == "ratio") {
		double value;
		if (attribute->QueryDoubleValue(&value) != TIXML_SUCCESS) {
			_stream << "Invalid value: " << attribute->Name() << "=" << attribute->Value() << " at row: " << attribute->Row() << ", col: " << attribute->Column();
			Error::_errMsg = _stream.str();
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			return 1;
		}
		if (!Validator::ElementOfAndContainsEndPoints(0.0, 100.0, value)) {
			_stream << "Value out of range: " << attribute->Name() << "=" << attribute->Value() << " at row: " << attribute->Row() << ", col: " << attribute->Column();
			Error::_errMsg = _stream.str();
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			return 1;
		}
		component->ratio = value;
	}
	else {
		_stream << "Unknown attribute: '" << attribute->Name() << "' at row: " << attribute->Row() << ", col: " << attribute->Column();
		Error::_errMsg = _stream.str();
		Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
		return 1;
	}

	return 0;
}

int XmlFileAdapter::DeserializeNebula(TiXmlElement *xmlElement, Nebula *nebula)
{
	// Iterate over the attributes
	for (TiXmlAttribute *attribute = xmlElement->FirstAttribute(); attribute; attribute = attribute->Next() ) {
		if (DeserializeNebulaAttributes(attribute, nebula) == 1) {
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			return 1;
		}
	}

	std::string text;
	double value = 0.0;
	TiXmlNode *child = xmlElement->FirstChild("MassFactor");
	if (child != 0) {
		TiXmlNode *node = child->ToElement();
		if (node == 0) {
			_stream << "Invalid xml element at row: " << xmlElement->Row() << ", col: " << xmlElement->Column();
			Error::_errMsg = _stream.str();
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			return 1;
		}
		// TODO: howto check whether the tag contains data or not. If it is empty the next line will cause on error!!
		text = child->ToElement()->GetText();
		if (text.length() > 0) {
			// TODO: implement check for number
			value = atof(text.c_str());
			if (!Validator::GreaterThan(0.0, value)) {
				_stream << "Value out of range at row: " << child->Row() << ", col: " << child->Column();
				Error::_errMsg = _stream.str();
				Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
				return 1;
			}
			nebula->massFactor = value;
		} else {
			_stream << "Missing value at row: " << child->Row() << ", col: " << child->Column();
			Error::_errMsg = _stream.str();
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			return 1;
		}
	}

	child = xmlElement->FirstChild("GasToDustRatio");
	if (child != 0) {
		if (child->ToElement() == 0) {
			_stream << "Invalid xml element at row: " << xmlElement->Row() << ", col: " << xmlElement->Column();
			Error::_errMsg = _stream.str();
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			return 1;
		}
		text = child->ToElement()->GetText();
		if (text.length() > 0) {
			// TODO: implement check for number
			value = atof(text.c_str());
			if (!Validator::GreaterThan(0.0, value)) {
				_stream << "Value out of range at row: " << child->Row() << ", col: " << child->Column();
				Error::_errMsg = _stream.str();
				Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
				return 1;
			}
			nebula->gasToDustRatio = value;
		} else {
			_stream << "Missing value at row: " << child->Row() << ", col: " << child->Column();
			Error::_errMsg = _stream.str();
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			return 1;
		}
	}

	child = xmlElement->FirstChild("SnowLine");
	if (child != 0) {
		if (child->ToElement() == 0) {
			_stream << "Invalid xml element at row: " << xmlElement->Row() << ", col: " << xmlElement->Column();
			Error::_errMsg = _stream.str();
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			return 1;
		}
		text = child->ToElement()->GetText();
		if (text.length() > 0) {
			// TODO: implement check for number
			value = atof(text.c_str());
			if (!Validator::GreaterThan(0.0, value)) {
				_stream << "Value out of range at row: " << child->Row() << ", col: " << child->Column();
				Error::_errMsg = _stream.str();
				Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
				return 1;
			}
			nebula->snowLine = value;
		} else {
			_stream << "Missing value at row: " << child->Row() << ", col: " << child->Column();
			Error::_errMsg = _stream.str();
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			return 1;
		}
	}

	child = xmlElement->FirstChild("SolidsComponent");
	if (child != 0) {
		if (child->ToElement() == 0) {
			_stream << "Invalid xml element at row: " << xmlElement->Row() << ", col: " << xmlElement->Column();
			Error::_errMsg = _stream.str();
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			return 1;
		}
		if (DeserializeSolidsComponent(child->ToElement(), &(nebula->solidsComponent)) == 1)
			return 1;
	}

	child = xmlElement->FirstChild("GasComponent");
	if (child != 0) {
		if (child->ToElement() == 0) {
			_stream << "Invalid xml element at row: " << xmlElement->Row() << ", col: " << xmlElement->Column();
			Error::_errMsg = _stream.str();
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			return 1;
		}
		if (DeserializeGasComponent(child->ToElement(), &(nebula->gasComponent)) == 1)
			return 1;
	}

	return 0;
}

int XmlFileAdapter::DeserializeNebulaAttributes(TiXmlAttribute *attribute, Nebula *nebula)
{
	std::string attributeName = attribute->Name();
	std::string attributeValue = attribute->Value();
	// Transform the name to lowercase
	std::transform(attributeName.begin(), attributeName.end(), attributeName.begin(), ::tolower);
	if (     attributeName == "name") {
		nebula->name = attributeValue;
	}
	else if (attributeName == "description") {
		nebula->description = attributeValue;
	}
	else if (attributeName == "path") {
		nebula->path = attributeValue;
	}
	else {
		_stream << "Unknown attribute: '" << attribute->Name() << "' at row: " << attribute->Row() << ", col: " << attribute->Column();
		Error::_errMsg = _stream.str();
		Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
		return 1;
	}

	return 0;
}

int XmlFileAdapter::DeserializeSolidsComponent(TiXmlElement *xmlElement, SolidsComponent *solidsComponent)
{
	std::string text;
	double value = 0.0;
	TiXmlNode *child = xmlElement->FirstChild("IceCondensationFactor");
	if (child != 0) {
		if (child->ToElement() == 0) {
			_stream << "Invalid xml element at row: " << xmlElement->Row() << ", col: " << xmlElement->Column();
			Error::_errMsg = _stream.str();
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			return 1;
		}
		text = child->ToElement()->GetText();
		if (text.length() > 0) {
			// TODO: implement check for number
			value = atof(text.c_str());
			if (!Validator::GreaterThan(0.0, value)) {
				_stream << "Value out of range at row: " << child->Row() << ", col: " << child->Column();
				Error::_errMsg = _stream.str();
				Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
				return 1;
			}
			solidsComponent->SetIceCondensationFactor(value);
		} else {
			_stream << "Missing value at row: " << child->Row() << ", col: " << child->Column();
			Error::_errMsg = _stream.str();
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			return 1;
		}
	}

	child = xmlElement->FirstChild("SolidsDensityFunction");
	if (child != 0) {
		if (child->ToElement() == 0) {
			_stream << "Invalid xml element at row: " << xmlElement->Row() << ", col: " << xmlElement->Column();
			Error::_errMsg = _stream.str();
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			return 1;
		}
		if (DeserializeSolidsDensityFunction(child->ToElement(), solidsComponent) == 1) {
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			return 1;
		}
	}

	return 0;
}

int XmlFileAdapter::DeserializeSolidsDensityFunction(TiXmlElement *xmlElement, SolidsComponent *solidsComponent)
{
	// Iterate over the attributes
	for (TiXmlAttribute *attribute = xmlElement->FirstAttribute(); attribute; attribute = attribute->Next() ) {
		if (DeserializeSolidsDensityFunctionAttributes(attribute, solidsComponent) == 1) {
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			return 1;
		}
	}

	TiXmlNode *child = xmlElement->FirstChild("Density");
	if (child != 0) {
		if (child->ToElement() == 0) {
			_stream << "Invalid xml element at row: " << xmlElement->Row() << ", col: " << xmlElement->Column();
			Error::_errMsg = _stream.str();
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			return 1;
		}
		double value = 0.0;
		std::string unit;
		if (DeserializeDimensionalQuantity(child->ToElement(), &value, &unit) == 1) {
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			return 1;
		}
		if (!ValidSurfaceDensityUnitAttribute(child->ToElement(), unit, "unit")) {
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			return 1;
		}
		// unit conversion
		UnitTool::SurfaceDensityToCM(unit, value);
		solidsComponent->GetSolidsDensityFunction()->c = value;
	}

	return 0;
}

int XmlFileAdapter::DeserializeSolidsDensityFunctionAttributes(TiXmlAttribute *attribute, SolidsComponent *solidsComponent)
{
	std::string attributeName = attribute->Name();
	std::string attributeValue = attribute->Value();
	// Transform the name to lowercase
	std::transform(attributeName.begin(), attributeName.end(), attributeName.begin(), ::tolower);
	if (attributeName == "index") {
		double value;
		if (attribute->QueryDoubleValue(&value) != TIXML_SUCCESS) {
			_stream << "Invalid value: " << attribute->Name() << "=" << attribute->Value() << " at row: " << attribute->Row() << ", col: " << attribute->Column();
			Error::_errMsg = _stream.str();
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			return 1;
		}
		solidsComponent->GetSolidsDensityFunction()->index = value;
	}
	else {
		_stream << "Unknown attribute: '" << attribute->Name() << "' at row: " << attribute->Row() << ", col: " << attribute->Column();
		Error::_errMsg = _stream.str();
		Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
		return 1;
	}

	return 0;
}

int XmlFileAdapter::DeserializeGasComponent(TiXmlElement *xmlElement, GasComponent *gasComponent)
{
	std::string unitString;
	if (!ValidTimeUnitAttribute(xmlElement, unitString, "unit")) {
		Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
		return 1;
	}

	// Iterate over the attributes
	for (TiXmlAttribute *attribute = xmlElement->FirstAttribute(); attribute; attribute = attribute->Next() ) {
		if (DeserializeGasComponentAttributes(attribute, gasComponent, unitString) == 1) {
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			return 1;
		}
	}

	TiXmlNode *child = xmlElement->FirstChild("Eta");
	if (child != 0) {
		if (child->ToElement() == 0) {
			_stream << "Invalid xml element at row: " << xmlElement->Row() << ", col: " << xmlElement->Column();
			Error::_errMsg = _stream.str();
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			return 1;
		}
		if (DeserializePowerLaw(child->ToElement(), &(gasComponent->eta)) == 1) {
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			return 1;
		}
	}

	child = xmlElement->FirstChild("ScaleHeight");
	if (child != 0) {
		if (child->ToElement() == 0) {
			_stream << "Invalid xml element at row: " << xmlElement->Row() << ", col: " << xmlElement->Column();
			Error::_errMsg = _stream.str();
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			return 1;
		}
		if (DeserializePowerLaw(child->ToElement(), &(gasComponent->scaleHeight)) == 1) {
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			return 1;
		}
	}

	child = xmlElement->FirstChild("Tau");
	if (child != 0) {
		if (child->ToElement() == 0) {
			_stream << "Invalid xml element at row: " << xmlElement->Row() << ", col: " << xmlElement->Column();
			Error::_errMsg = _stream.str();
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			return 1;
		}
		if (DeserializePowerLaw(child->ToElement(), &(gasComponent->tau)) == 1) {
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			return 1;
		}
	}

	child = xmlElement->FirstChild("GasDensityFunction");
	if (child != 0) {
		if (child->ToElement() == 0) {
			_stream << "Invalid xml element at row: " << xmlElement->Row() << ", col: " << xmlElement->Column();
			Error::_errMsg = _stream.str();
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			return 1;
		}
		if (DeserializeGasDensityFunction(child->ToElement(), gasComponent) == 1) {
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			return 1;
		}
	}

	return 0;
}

int XmlFileAdapter::DeserializeGasComponentAttributes(TiXmlAttribute *attribute, GasComponent *gasComponent, std::string unitString)
{
	double value = 0.0;

	std::string attributeName = attribute->Name();
	std::string attributeValue = attribute->Value();
	// Transform the name to lowercase
	std::transform(attributeName.begin(), attributeName.end(), attributeName.begin(), ::tolower);
	if (     attributeName == "alpha") {
		if (attribute->QueryDoubleValue(&value) != TIXML_SUCCESS) {
			_stream << "Invalid value: " << attribute->Name() << "=" << attribute->Value() << " at row: " << attribute->Row() << ", col: " << attribute->Column();
			Error::_errMsg = _stream.str();
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			return 1;
		}
		if (!Validator::GreaterThan(0.0, value)) {
			_stream << "Value out of range: " << attribute->Name() << "=" << attribute->Value() << " at row: " << attribute->Row() << ", col: " << attribute->Column();
			Error::_errMsg = _stream.str();
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			return 1;
		}
		gasComponent->alpha = value;
	}
	else if (attributeName == "type") {
		std::string type = attributeValue;
		std::transform(type.begin(), type.end(), type.begin(), ::tolower);
		if (SetGasDecreaseType(type, gasComponent) == 1) {
			_stream << "Unknown type: " << attribute->Name() << "=" << attribute->Value() << " at row: " << attribute->Row() << ", col: " << attribute->Column();
			Error::_errMsg = _stream.str();
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			return 1;
		}
	}
	else if (attributeName == "timescale") {
		if (attribute->QueryDoubleValue(&value) != TIXML_SUCCESS) {
			_stream << "Invalid value: " << attribute->Name() << "=" << attribute->Value() << " at row: " << attribute->Row() << ", col: " << attribute->Column();
			Error::_errMsg = _stream.str();
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			return 1;
		}
		if (!Validator::GreaterThanOrEqualTo(0.0, value)) {
			_stream << "Value out of range: " << attribute->Name() << "=" << attribute->Value() << " at row: " << attribute->Row() << ", col: " << attribute->Column();
			Error::_errMsg = _stream.str();
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			return 1;
		}
		// Convert to day
		UnitTool::TimeToDay(unitString, value);
		gasComponent->timeScale = value;
	}
	else if (attributeName == "t0") {
		if (attribute->QueryDoubleValue(&value) != TIXML_SUCCESS) {
			_stream << "Invalid value: " << attribute->Name() << "=" << attribute->Value() << " at row: " << attribute->Row() << ", col: " << attribute->Column();
			Error::_errMsg = _stream.str();
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			return 1;
		}
		// Convert to day
		UnitTool::TimeToDay(unitString, value);
		gasComponent->t0 = value;
	}
	else if (attributeName == "t1") {
		if (attribute->QueryDoubleValue(&value) != TIXML_SUCCESS) {
			_stream << "Invalid value: " << attribute->Name() << "=" << attribute->Value() << " at row: " << attribute->Row() << ", col: " << attribute->Column();
			Error::_errMsg = _stream.str();
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			return 1;
		}
		// Convert to day
		UnitTool::TimeToDay(unitString, value);
		gasComponent->t1 = value;
	}
	else if (attributeName == "unit") {
		;
	}
	else {
		_stream << "Unknown attribute: '" << attribute->Name() << "' at row: " << attribute->Row() << ", col: " << attribute->Column();
		Error::_errMsg = _stream.str();
		Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
		return 1;
	}

	return 0;
}

int XmlFileAdapter::SetGasDecreaseType(std::string type, GasComponent *gasComponent)
{
	if (     type == "constant")	{ gasComponent->type = Constant;	}
	else if (type == "linear")		{ gasComponent->type = Linear;		}
	else if (type == "exponential")	{ gasComponent->type = Exponential; }
	else {
		return 1;
	}

	return 0;
}

int XmlFileAdapter::DeserializeGasDensityFunction(TiXmlElement *xmlElement, GasComponent *gasComponent)
{
	// Iterate over the attributes
	for (TiXmlAttribute *attribute = xmlElement->FirstAttribute(); attribute; attribute = attribute->Next() ) {
		if (DeserializeGasDensityFunctionAttributes(attribute, gasComponent) == 1) {
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			return 1;
		}
	}

	TiXmlNode *child = xmlElement->FirstChild("Density");
	if (child != 0) {
		if (child->ToElement() == 0) {
			_stream << "Invalid xml element at row: " << xmlElement->Row() << ", col: " << xmlElement->Column();
			Error::_errMsg = _stream.str();
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			return 1;
		}
		double value = 0.0;
		std::string unit;
		if (DeserializeDimensionalQuantity(child->ToElement(), &value, &unit) == 1) {
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			return 1;
		}
		if (!ValidVolumeDensityUnitAttribute(child->ToElement(), unit, "unit")) {
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			return 1;
		}
		// unit conversion
		UnitTool::VolumeDensityToCM(unit, value);
		gasComponent->gasDensityFunction.c = value;
	}

	return 0;
}

int XmlFileAdapter::DeserializeGasDensityFunctionAttributes(TiXmlAttribute *attribute, GasComponent *gasComponent)
{
	std::string attributeName = attribute->Name();
	std::string attributeValue = attribute->Value();
	// Transform the name to lowercase
	std::transform(attributeName.begin(), attributeName.end(), attributeName.begin(), ::tolower);
	if (attributeName == "index") {
		double value;
		if (attribute->QueryDoubleValue(&value) != TIXML_SUCCESS) {
			_stream << "Invalid value: " << attribute->Name() << "=" << attribute->Value() << " at row: " << attribute->Row() << ", col: " << attribute->Column();
			Error::_errMsg = _stream.str();
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			return 1;
		}
		gasComponent->gasDensityFunction.index = value;
	}
	else {
		_stream << "Unknown attribute: '" << attribute->Name() << "' at row: " << attribute->Row() << ", col: " << attribute->Column();
		Error::_errMsg = _stream.str();
		Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
		return 1;
	}

	return 0;
}

int XmlFileAdapter::DeserializePowerLaw(TiXmlElement *xmlElement, PowerLaw *powerLaw)
{
	double value = 0.0;
	// Iterate over the attributes
	for (TiXmlAttribute *attribute = xmlElement->FirstAttribute(); attribute; attribute = attribute->Next() ) {
		std::string attributeName = attribute->Name();
		std::string attributeValue = attribute->Value();
		// Transform the name to lowercase
		std::transform(attributeName.begin(), attributeName.end(), attributeName.begin(), ::tolower);
		if (     attributeName == "c") {
			if (attribute->QueryDoubleValue(&value) != TIXML_SUCCESS) {
				_stream << "Invalid value: " << attribute->Name() << "=" << attribute->Value() << " at row: " << attribute->Row() << ", col: " << attribute->Column();
				Error::_errMsg = _stream.str();
				Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
				return 1;
			}
			powerLaw->c = value;
		}
		else if (attributeName == "index") {
			if (attribute->QueryDoubleValue(&value) != TIXML_SUCCESS) {
				_stream << "Invalid value: " << attribute->Name() << "=" << attribute->Value() << " at row: " << attribute->Row() << ", col: " << attribute->Column();
				Error::_errMsg = _stream.str();
				Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
				return 1;
			}
			powerLaw->index = value;
		}
		else {
			_stream << "Unknown attribute: '" << attribute->Name() << "' at row: " << attribute->Row() << ", col: " << attribute->Column();
			Error::_errMsg = _stream.str();
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			return 1;
		}
	}

	return 0;
}

int XmlFileAdapter::DeserializeDimensionalQuantity(TiXmlElement *xmlElement, double* value, std::string* unit)
{
	// Iterate over the attributes
	for (TiXmlAttribute *attribute = xmlElement->FirstAttribute(); attribute; attribute = attribute->Next() ) {
		std::string attributeName = attribute->Name();
		std::string attributeValue = attribute->Value();
		// Transform the name to lowercase
		std::transform(attributeName.begin(), attributeName.end(), attributeName.begin(), ::tolower);
		if (     attributeName == "value") {
			if (attribute->QueryDoubleValue(value) != TIXML_SUCCESS) {
				_stream << "Invalid value: " << attribute->Name() << "=" << attribute->Value() << " at row: " << attribute->Row() << ", col: " << attribute->Column();
				Error::_errMsg = _stream.str();
				Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
				return 1;
			}
			if (!Validator::GreaterThanOrEqualTo(0.0, *value)) {
				_stream << "Value out of range: " << attribute->Name() << "=" << attribute->Value() << " at row: " << attribute->Row() << ", col: " << attribute->Column();
				Error::_errMsg = _stream.str();
				Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
				return 1;
			}
		}
		else if (attributeName == "unit") {
			*unit = attributeValue;
		}
		else {
			_stream << "Unknown attribute: '" << attribute->Name() << "' at row: " << attribute->Row() << ", col: " << attribute->Column();
			Error::_errMsg = _stream.str();
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			return 1;
		}
	}

	return 0;
}

bool XmlFileAdapter::ValidDistanceUnitAttribute(TiXmlElement *xmlElement, std::string& unitValue, std::string unitName)
{
	for (TiXmlAttribute *attribute = xmlElement->FirstAttribute(); attribute; attribute = attribute->Next() ) {
		double dummy = 0.0;
		std::string attributeName = attribute->Name();
		std::string attributeValue = attribute->Value();
		// Transform the name and value to lowercase
		std::transform(attributeName.begin(), attributeName.end(), attributeName.begin(), ::tolower);
		std::transform(attributeValue.begin(), attributeValue.end(), attributeValue.begin(), ::tolower);
		if (attributeName == unitName) {
			if (UnitTool::DistanceToAu(attributeValue, dummy) == 1) {
				_stream << "Unrecognized dimension for attribute '" << attribute->Name() << "' at row: " << attribute->Row() << ", col: " << attribute->Column();
				Error::_errMsg = _stream.str();
				return false;
			}
			unitValue = attributeValue;
			break;
		}
	}
	// If unit was not defined use the default
	if (unitValue.length() == 0 ) unitValue = "au";

	return true;
}

bool XmlFileAdapter::ValidAngleUnitAttribute(TiXmlElement *xmlElement, std::string& unitValue, std::string unitName)
{
	for (TiXmlAttribute *attribute = xmlElement->FirstAttribute(); attribute; attribute = attribute->Next() ) {
		double dummy = 0.0;
		std::string attributeName = attribute->Name();
		std::string attributeValue = attribute->Value();
		// Transform the name and value to lowercase
		std::transform(attributeName.begin(), attributeName.end(), attributeName.begin(), ::tolower);
		std::transform(attributeValue.begin(), attributeValue.end(), attributeValue.begin(), ::tolower);
		if (attributeName == unitName) {
			if (UnitTool::AngleToRadian(attributeValue, dummy) == 1) {
				_stream << "Unrecognized dimension for attribute '" << attribute->Name() << "' at row: " << attribute->Row() << ", col: " << attribute->Column();
				Error::_errMsg = _stream.str();
				return false;
			}
			unitValue = attributeValue;
			break;
		}
	}
	// If unit was not defined use the default
	if (unitValue.length() == 0 ) unitValue = "radian";

	return true;
}

bool XmlFileAdapter::ValidTimeUnitAttribute(TiXmlElement *xmlElement, std::string& unitValue, std::string unitName)
{
	for (TiXmlAttribute *attribute = xmlElement->FirstAttribute(); attribute; attribute = attribute->Next() ) {
		double dummy = 0.0;
		std::string attributeName = attribute->Name();
		std::string attributeValue = attribute->Value();
		// Transform the name and value to lowercase
		std::transform(attributeName.begin(), attributeName.end(), attributeName.begin(), ::tolower);
		std::transform(attributeValue.begin(), attributeValue.end(), attributeValue.begin(), ::tolower);
		if (attributeName == unitName) {
			if (UnitTool::TimeToDay(attributeValue, dummy) == 1) {
				_stream << "Unrecognized dimension for attribute '" << attribute->Name() << "' at row: " << attribute->Row() << ", col: " << attribute->Column();
				Error::_errMsg = _stream.str();
				return false;
			}
			unitValue = attributeValue;
			break;
		}
	}
	// If unit was not defined use the default
	if (unitValue.length() == 0 ) unitValue = "day";

	return true;
}

bool XmlFileAdapter::ValidVelocityUnitAttribute(TiXmlElement *xmlElement, std::string& unitValue, std::string unitName)
{
	for (TiXmlAttribute *attribute = xmlElement->FirstAttribute(); attribute; attribute = attribute->Next() ) {
		double dummy = 0.0;
		std::string attributeName = attribute->Name();
		std::string attributeValue = attribute->Value();
		// Transform the name and value to lowercase
		std::transform(attributeName.begin(), attributeName.end(), attributeName.begin(), ::tolower);
		std::transform(attributeValue.begin(), attributeValue.end(), attributeValue.begin(), ::tolower);
		if (attributeName == unitName) {
			if (UnitTool::VelocityToCM(attributeValue, dummy) == 1) {
				_stream << "Unrecognized dimension for attribute '" << attribute->Name() << "' at row: " << attribute->Row() << ", col: " << attribute->Column();
				Error::_errMsg = _stream.str();
				return false;
			}
			unitValue = attributeValue;
			break;
		}
	}
	// If unit was not defined use the default
	if (unitValue.length() == 0 ) unitValue = "auday";

	return true;
}

bool XmlFileAdapter::ValidMassUnitAttribute(TiXmlElement *xmlElement, std::string& unitValue, std::string unitName)
{
	for (TiXmlAttribute *attribute = xmlElement->FirstAttribute(); attribute; attribute = attribute->Next() ) {
		double dummy = 0.0;
		std::string attributeName = attribute->Name();
		std::string attributeValue = attribute->Value();
		// Transform the name and value to lowercase
		std::transform(attributeName.begin(), attributeName.end(), attributeName.begin(), ::tolower);
		std::transform(attributeValue.begin(), attributeValue.end(), attributeValue.begin(), ::tolower);
		if (attributeName == unitName) {
			if (UnitTool::MassToSolar(attributeValue, dummy) == 1) {
				_stream << "Unrecognized dimension for attribute '" << attribute->Name() << "' at row: " << attribute->Row() << ", col: " << attribute->Column();
				Error::_errMsg = _stream.str();
				return false;
			}
			unitValue = attributeValue;
			break;
		}
	}
	// If unit was not defined use the default
	if (unitValue.length() == 0 ) unitValue = "solar";

	return true;
}

bool XmlFileAdapter::ValidVolumeDensityUnitAttribute(TiXmlElement *xmlElement, std::string& unitValue, std::string unitName)
{
	for (TiXmlAttribute *attribute = xmlElement->FirstAttribute(); attribute; attribute = attribute->Next() ) {
		double dummy = 0.0;
		std::string attributeName = attribute->Name();
		std::string attributeValue = attribute->Value();
		// Transform the name and value to lowercase
		std::transform(attributeName.begin(), attributeName.end(), attributeName.begin(), ::tolower);
		std::transform(attributeValue.begin(), attributeValue.end(), attributeValue.begin(), ::tolower);
		if (attributeName == unitName) {
			if (UnitTool::VolumeDensityToCM(attributeValue, dummy) == 1) {
				_stream << "Unrecognized dimension for attribute '" << attribute->Name() << "' at row: " << attribute->Row() << ", col: " << attribute->Column();
				Error::_errMsg = _stream.str();
				return false;
			}
			unitValue = attributeValue;
			break;
		}
	}
	// If unit was not defined use the default
	if (unitValue.length() == 0 ) unitValue = "solarau3";

	return true;
}

bool XmlFileAdapter::ValidSurfaceDensityUnitAttribute(TiXmlElement *xmlElement, std::string& unitValue, std::string unitName)
{
	for (TiXmlAttribute *attribute = xmlElement->FirstAttribute(); attribute; attribute = attribute->Next() ) {
		double dummy = 0.0;
		std::string attributeName = attribute->Name();
		std::string attributeValue = attribute->Value();
		// Transform the name and value to lowercase
		std::transform(attributeName.begin(), attributeName.end(), attributeName.begin(), ::tolower);
		std::transform(attributeValue.begin(), attributeValue.end(), attributeValue.begin(), ::tolower);
		if (attributeName == unitName) {
			if (UnitTool::SurfaceDensityToCM(attributeValue, dummy) == 1) {
				_stream << "Unrecognized dimension for attribute '" << attribute->Name() << "' at row: " << attribute->Row() << ", col: " << attribute->Column();
				Error::_errMsg = _stream.str();
				return false;
			}
			unitValue = attributeValue;
			break;
		}
	}
	// If unit was not defined use the default
	if (unitValue.length() == 0 ) unitValue = "solarau2";

	return true;
}
