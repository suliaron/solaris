#include <cmath>
#include <iostream>
#include <fstream>

#include "Acceleration.h"
#include "BinaryFileAdapter.h"
#include "Body.h"
#include "BodyData.h"
#include "BodyGroup.h"
#include "BodyGroupList.h"
#include "Calculate.h"
#include "Characteristics.h"
#include "Component.h"
#include "Constants.h"
#include "Counter.h"
#include "DateTime.h"
#include "DormandPrince.h"
#include "DragCoefficient.h"
#include "Ephemeris.h"
#include "Error.h"
#include "Event.h"
#include "EventCondition.h"
#include "FargoParameters.h"
#include "GasComponent.h"
#include "GasDecreaseType.h"
#include "Integrator.h"
#include "IntegratorType.h"
#include "NBodies.h"
#include "Nebula.h"
#include "OrbitalElement.h"
#include "Output.h"
#include "Phase.h"
#include "PowerLaw.h"
#include "RungeKutta4.h"
#include "RungeKutta56.h"
#include "RungeKuttaFehlberg78.h"
#include "Settings.h"
#include "Simulation.h"
#include "Simulator.h"
#include "SolidsComponent.h"
#include "TimeLine.h"
#include "tinystr.h"
#include "tinyxml.h"
#include "Tokenizer.h"
#include "Tools.h"
#include "TwoBodyAffair.h"
#include "Units.h"
#include "Validator.h"
#include "Vector.h"
#include "XmlFileAdapter.h"

#define SQR(a)		((a)*(a))
#define CUBE(a)		((a)*(a)*(a))
#define FORTH(a)	((a)*(a)*(a)*(a))
#define FIFTH(a)	((a)*(a)*(a)*(a)*(a))

//#define EPSTEIN_STOKES_VERSUS_DISTANCE

//#define TEST_CALCULATE
//#define TEST_CONSTANTS
//#define TEST_TYPEI_MIGRATION_TIME
//#define TEST_MEAN_FREE_PATH
//#define TEST_BODYGROUP
//#define TEST_BODYGROUPLIST
#define TEST_TOOLS

//#define TEST_DUSTPARTICLE

bool TestConstants()
{
	bool	failed = false;
	double	value = 0.0;

	// GramPerCm3ToSolarPerAu3
	{
		value = Constants::GramPerCm3ToSolarPerAu3;
		if (fabs(value - 1.6831291259964250163336161614994e6) > 1.0e-10 )
		{
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			failed = true;
		}
	}

	return failed;		
}

bool TestTypeIMigrationTime()
{
	bool	failed = false;
	double	value = 0.0;

	// 
	{
		GasComponent gasComponent;

		double	Mc = 1.0;		// Solar mass
		double	m = 0.5*Constants::EarthToSolar;	// Solar mass
		double	beta = 1.5;	
		double	Sigma_gas = 150.0*Constants::GramPerCm2ToSolarPerAu2;	// Solar mass / AU^2
		double	a = 1.0; // AU
		double	e = 0.0;
		double	r = 2.0;
		double	hr = 0.056;
		double	frac = 1.0;//(1.0 + FIFTH(e*r/(1.3*gasComponent.scaleHeight(r))) )/(1.0 - FORTH(e*r/(1.1*gasComponent.scaleHeight(r))) );
		double	T = 2.0*Constants::Pi/Constants::Gauss*sqrt(CUBE(a)/(Mc+m));
		double	O = 2.0*Constants::Pi/T;
		double	O1 = Constants::Gauss * sqrt((Mc + m)/CUBE(a));

		double a1 = 2.0/(2.7+beta);
		double a2 = (Mc/m);
		double a3 = Mc/(Sigma_gas*a*a);
		double a4 = SQR(hr);

		double tm = a1*a2*a3*a4*frac/O;
		tm *= Constants::DayToYear;
		tm = 1.0*tm;

	}

	return failed;		
}

bool TestCalculate()
{
	bool	failed = false;
	double	value = 0.0;

	// Test Calculate::TotalMass()
	{
		BodyData bodyData;

		bodyData.nBodies.giantPlanet = 3;
		bodyData.nBodies.total = 3;
		bodyData.Allocate();

		for (int i=0; i<bodyData.nBodies.total; i++)
			bodyData.mass[i] = 1.0/bodyData.nBodies.total;

		double result = 0.0;
		Calculate::TotalMass(&bodyData, result);

		if (fabs(result - 1.0) > 1.0e-10 )
		{
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			failed = true;
		}
	}

	// Test Calculate::PhaseOfBC()
	{
		BodyData bodyData;

		bodyData.nBodies.giantPlanet = 3;
		bodyData.nBodies.total = 3;
		bodyData.Allocate();

		for (int i=0; i<bodyData.nBodies.total; i++) {
			bodyData.mass[i] = 1.0/bodyData.nBodies.total;
		}
		bodyData.y0[0] = 0.0;
		bodyData.y0[1] = 0.0;
		bodyData.y0[2] = 0.0;
		bodyData.y0[3] = 0.0;
		bodyData.y0[4] = 0.0;
		bodyData.y0[5] = 0.0;

		bodyData.y0[6] = 1.0;
		bodyData.y0[7] = 0.0;
		bodyData.y0[8] = 0.0;
		bodyData.y0[9] = 0.0;
		bodyData.y0[10] = 0.0;
		bodyData.y0[11] = 0.0;

		bodyData.y0[12] = 0.5;
		bodyData.y0[13] = sqrt(3.0)/2.0;
		bodyData.y0[14] = 0.0;
		bodyData.y0[15] = 0.0;
		bodyData.y0[16] = 0.0;
		bodyData.y0[17] = 0.0;

		double bc[6] = { 1.0/2.0, 1.0/(2.0*sqrt(3.0)), 0.0, 0.0, 0.0, 0.0 };
		double result[6] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
		Calculate::PhaseOfBC(&bodyData, result);

		for (int i=0; i<6; i++) {
			if (fabs(result[i] - bc[i]) > 1.0e-14 )
			{
				Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
				failed = true;
			}
		}
	}

	// Test Calculate::PhaseWithRespectToBC()
	{
		BodyData bodyData;

		bodyData.nBodies.giantPlanet = 3;
		bodyData.nBodies.total = 3;
		bodyData.Allocate();

		for (int i=0; i<bodyData.nBodies.total; i++) {
			bodyData.mass[i] = 1.0/bodyData.nBodies.total;
		}
		bodyData.y0[0] = 0.0;
		bodyData.y0[1] = 0.0;
		bodyData.y0[2] = 0.0;
		bodyData.y0[3] = 0.0;
		bodyData.y0[4] = 0.0;
		bodyData.y0[5] = 0.0;

		bodyData.y0[6] = 1.0;
		bodyData.y0[7] = 0.0;
		bodyData.y0[8] = 0.0;
		bodyData.y0[9] = 0.0;
		bodyData.y0[10] = 0.0;
		bodyData.y0[11] = 0.0;

		bodyData.y0[12] = 0.5;
		bodyData.y0[13] = sqrt(3.0)/2.0;
		bodyData.y0[14] = 0.0;
		bodyData.y0[15] = 0.0;
		bodyData.y0[16] = 0.0;
		bodyData.y0[17] = 0.0;

		bodyData.y[0] = -0.5;
		bodyData.y[1] = -1.0/(2.0*sqrt(3.0));
		bodyData.y[2] = 0.0;
		bodyData.y[3] = 0.0;
		bodyData.y[4] = 0.0;
		bodyData.y[5] = 0.0;

		bodyData.y[6] = 0.5;
		bodyData.y[7] = -1.0/(2.0*sqrt(3.0));
		bodyData.y[8] = 0.0;
		bodyData.y[9] = 0.0;
		bodyData.y[10] = 0.0;
		bodyData.y[11] = 0.0;

		bodyData.y[12] = 0.0;
		bodyData.y[13] = 1.0/sqrt(3.0);
		bodyData.y[14] = 0.0;
		bodyData.y[15] = 0.0;
		bodyData.y[16] = 0.0;
		bodyData.y[17] = 0.0;

		double result[6] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
		Calculate::PhaseOfBC(&bodyData, result);
		Calculate::PhaseWithRespectToBC(&bodyData, result);

		int nvar = bodyData.nBodies.NOfVar();
		for (int i=0; i<nvar; i++) {
			if (fabs(bodyData.y[i] - bodyData.y0[i]) > 1.0e-14 )
			{
				Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
				failed = true;
			}
		}

		Calculate::PhaseOfBC(&bodyData, result);
		for (int i=0; i<6; i++) {
			if (fabs(result[i] - 0.0) > 1.0e-14 )
			{
				Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
				failed = true;
			}
		}
	}

	// Test Calculate::AngularMomentum()
	{
		BodyData bodyData;

		bodyData.nBodies.giantPlanet = 3;
		bodyData.nBodies.total = 3;
		bodyData.Allocate();

		double m = 1.0/bodyData.nBodies.total;
		for (int i=0; i<bodyData.nBodies.total; i++) {
			bodyData.mass[i] = m;
		}
		bodyData.y0[0] = 0.0;
		bodyData.y0[1] = 0.0;
		bodyData.y0[2] = 0.0;
		bodyData.y0[3] = 0.0;
		bodyData.y0[4] = 0.0;
		bodyData.y0[5] = 0.0;

		bodyData.y0[6] = 1.0;
		bodyData.y0[7] = 0.0;
		bodyData.y0[8] = 0.0;
		bodyData.y0[9] = -bodyData.y0[7];
		bodyData.y0[10] = bodyData.y0[6];
		bodyData.y0[11] = 0.0;

		bodyData.y0[12] = 0.5;
		bodyData.y0[13] = sqrt(3.0)/2.0;
		bodyData.y0[14] = 0.0;
		bodyData.y0[15] = -bodyData.y0[13];
		bodyData.y0[16] = bodyData.y0[12];
		bodyData.y0[17] = 0.0;

		Vector result;
		Vector c(0.0, 0.0, 2*m);
		Calculate::AngularMomentum(&bodyData, &result);

		if (fabs(result.x - c.x) > 1.0e-14 )
		{
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			failed = true;
		}

		if (fabs(result.y - c.y) > 1.0e-14 )
		{
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			failed = true;
		}

		if (fabs(result.z - c.z) > 1.0e-14 )
		{
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			failed = true;
		}
	}

	// Test Calculate::PotentialEnergy()
	{
		BodyData bodyData;

		bodyData.nBodies.giantPlanet = 3;
		bodyData.nBodies.total = 3;
		bodyData.Allocate();

		double m = 1.0/bodyData.nBodies.total;
		for (int i=0; i<bodyData.nBodies.total; i++) {
			bodyData.mass[i] = m;
		}
		bodyData.y0[0] = 0.0;
		bodyData.y0[1] = 0.0;
		bodyData.y0[2] = 0.0;
		bodyData.y0[3] = 0.0;
		bodyData.y0[4] = 0.0;
		bodyData.y0[5] = 0.0;

		bodyData.y0[6] = 1.0;
		bodyData.y0[7] = 0.0;
		bodyData.y0[8] = 0.0;
		bodyData.y0[9] = -bodyData.y0[7];
		bodyData.y0[10] = bodyData.y0[6];
		bodyData.y0[11] = 0.0;

		bodyData.y0[12] = 0.5;
		bodyData.y0[13] = sqrt(3.0)/2.0;
		bodyData.y0[14] = 0.0;
		bodyData.y0[15] = -bodyData.y0[13];
		bodyData.y0[16] = bodyData.y0[12];
		bodyData.y0[17] = 0.0;

		double result;
		Calculate::PotentialEnergy(&bodyData, result);

		if (fabs(result - Constants::Gauss2/3.0) > 1.0e-14 )
		{
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			failed = true;
		}
	}

	// Test Calculate::KineticEnergy()
	{
		BodyData bodyData;

		bodyData.nBodies.giantPlanet = 3;
		bodyData.nBodies.total = 3;
		bodyData.Allocate();

		double m = 1.0/bodyData.nBodies.total;
		for (int i=0; i<bodyData.nBodies.total; i++) {
			bodyData.mass[i] = m;
		}
		bodyData.y0[0] = 0.0;
		bodyData.y0[1] = 0.0;
		bodyData.y0[2] = 0.0;
		bodyData.y0[3] = 0.0;
		bodyData.y0[4] = 0.0;
		bodyData.y0[5] = 0.0;

		bodyData.y0[6] = 1.0;
		bodyData.y0[7] = 0.0;
		bodyData.y0[8] = 0.0;
		bodyData.y0[9] = -bodyData.y0[7];
		bodyData.y0[10] = bodyData.y0[6];
		bodyData.y0[11] = 0.0;

		bodyData.y0[12] = 0.5;
		bodyData.y0[13] = sqrt(3.0)/2.0;
		bodyData.y0[14] = 0.0;
		bodyData.y0[15] = -bodyData.y0[13];
		bodyData.y0[16] = bodyData.y0[12];
		bodyData.y0[17] = 0.0;

		double result;
		Calculate::KineticEnergy(&bodyData, result);
		if (fabs(result - m) > 1.0e-14 )
		{
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			failed = true;
		}

		Calculate::Energy(&bodyData, result);
		if (fabs(result - 1.0/3.0*(1.0 - Constants::Gauss2)) > 1.0e-14 )
		{
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			failed = true;
		}
	}

	return failed;
}

bool TestMeanFreePath()
{
	bool	failed = false;
	double	value = 0.0;
    double  expected = 0.0;

	// Test mean free path
	{
        GasComponent gasComponent;

        value = gasComponent.MeanFreePath_SI(1.4e-6);
        expected = 2.3 * Constants::ProtonMass_SI / (sqrt(2.0) * Constants::Pi * SQR(3.0e-10) * 1.4e-6);
        if (fabs(value - expected) > 1.0e-5 )
		{
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			failed = true;
		}

        value = gasComponent.MeanFreePath_CMU(1.4e-9 * Constants::GramPerCm3ToSolarPerAu3);
        double value_SI = value * Constants::AuToMeter;
        expected = 2.3 * Constants::ProtonMass_CMU / (sqrt(2.0) * Constants::Pi * SQR(3.0e-10 * Constants::MeterToAu) * (1.4e-9 * Constants::GramPerCm3ToSolarPerAu3));
        double expected_SI = expected * Constants::AuToMeter;
        double lambda = gasComponent.meanFreePath.Evaluate(1.0);
        lambda *= Constants::AuToMeter;
        if (fabs(value - expected) > 1.0e-5 )
		{
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			failed = true;
		}
    }

    // Test sound speed
    {
        GasComponent gasComponent;
        value = gasComponent.SoundSpeed_SI(100.0);
        expected = 599.07171977237708467220571287692;       // m/s
        if (fabs(value - expected) > 1.0e-5 )
		{
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			failed = true;
		}

        value = gasComponent.SoundSpeed_CMU(100.0);
        value *= (Constants::AuToMeter / Constants::DayToSecond);
        expected = 599.07171977237708467220571287692;       // m/s
        if (fabs(value - expected) > 1.0e-5 )
		{
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			failed = true;
		}
    }

    // Test temperature
    {
        GasComponent gasComponent;
        value = gasComponent.Temperature_SI(1.0 * Constants::SolarToKilogram, 1.0 * Constants::AuToMeter);
        // If ch = 0.045
        //expected = 500.69882237231788725014278527085;       // K
        // If ch = 0.02
        expected = 98.903471085889933;                      // K
        if (fabs(value - expected) > 1.0e-5 )
		{
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			failed = true;
		}
        
        value = gasComponent.Temperature_CMU(1.0, 1.0);
        // If ch = 0.045
        //expected = 500.69882237231788725014278527085;		// K
        // If ch = 0.02
        expected = 98.903471085889933;                      // K
        if (fabs(value - expected) > 1.0 )
		{
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			failed = true;
		}

        value = gasComponent.Temperature_CMU(1.0, 0.1);
        // If ch = 0.045
        //expected = 500.69882237231788725014278527085;       // K
        // If ch = 0.02
        //expected = 989.03471085889933;                      // K
        expected = 312.67170117983335;                        // K
        if (fabs(value - expected) > 1.0 )
		{
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			failed = true;
		}
    }

	// Test mean thermal velocity
    {
        GasComponent gasComponent;
		value = gasComponent.MeanThermalSpeed_SI(1.0 * Constants::SolarToKilogram, 1.0 * Constants::AuToMeter);
        expected = 950.7244052398592; // m / s
        if (fabs(value - expected) > 1.0e-5 )
		{
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			failed = true;
		}

		value = gasComponent.MeanThermalSpeed_CMU(1.0, 1.0);
        expected = 950.7244052398592 * (Constants::MeterToAu)/(Constants::SecondToDay); // m / s
        if (fabs(value - expected) > 1.0e-5 )
		{
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			failed = true;
		}
    }

	return failed;		
}

// Test the functions of the BodyGroupList class
bool TestBodyGroupList()
{
	bool	failed = false;
	double	value = 0.0;

	// Test data structure
	{
		BodyGroupList bodyGroupList;

		BodyGroup bodyGroup0("89A79746-3D47-49A0-9026-157F61420376", "Test body group 0", "",          0.0, "J2000.0");
		BodyGroup bodyGroup1("8C4C23A0-D162-419D-95EB-46F5034572A1", "Test body group 1", "20010101",  0.0, "J2000.0");
		BodyGroup bodyGroup2("A130417E-D47C-43F7-B6D9-A1FD8CCCAAB0", "Test body group 2", "20020101", 10.0, "J2000.0");

		Body centralBody(CentralBody);
		centralBody.characteristics = new Characteristics(1.0);

		Body protoPlanet(ProtoPlanet);
		protoPlanet.characteristics = new Characteristics(1.0 * Constants::EarthMoonToSolar);
	
		Body rockyPlanet(RockyPlanet);
		rockyPlanet.characteristics = new Characteristics(1.0 * Constants::EarthToSolar);
	
		Body giantPlanet(GiantPlanet);
		giantPlanet.characteristics = new Characteristics(1.0 * Constants::JupiterToSolar);
	
		Body superPl(SuperPlanetesimal);
		superPl.characteristics = new Characteristics(1.0e3 * Constants::KilogramToSolar);
	
		Body planetesimal(Planetesimal);
		planetesimal.characteristics = new Characteristics(1.0 * Constants::KilogramToSolar);
	
		Body testParticle(TestParticle);

		bodyGroup0.items.push_back(centralBody);
		bodyGroup0.items.push_back(protoPlanet);
		bodyGroup0.items.push_back(rockyPlanet);

		bodyGroup1.items.push_back(giantPlanet);
		bodyGroup1.items.push_back(superPl);

		bodyGroup2.items.push_back(planetesimal);
		bodyGroup2.items.push_back(testParticle);

		bodyGroupList.items.push_back(bodyGroup0);
		bodyGroupList.items.push_back(bodyGroup1);
		bodyGroupList.items.push_back(bodyGroup2);

		BodyGroup bg = bodyGroupList.items.front();
		Body body = bodyGroupList.items.front().items.front();
	}


    // DistinctEpochs()
    {
        BodyGroupList bodyGroupList;

        std::list<double> epochs;
        bool increasing = true;

        int result = bodyGroupList.DistinctEpochs(epochs, increasing);
        if (result != 0)
        {
			std::cout << "Failed at bodyGroupList.DistinctEpochs(epochs, increasing)" << std::endl;
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			failed = true;
        }

        BodyGroup bg0;
        bg0.epoch = "19900101";
        BodyGroup bg1;
        bg1.epoch = "19910101";
        BodyGroup bg2;
        bg2.epoch = "19930101";
        BodyGroup bg3;
        bg3.epoch = "19930101";

        bodyGroupList.items.push_back(bg0);
        bodyGroupList.items.push_back(bg1);
        bodyGroupList.items.push_back(bg2);
        bodyGroupList.items.push_back(bg3);

        result = bodyGroupList.DistinctEpochs(epochs, increasing);
        if (result != 0)
        {
			std::cout << "Failed at bodyGroupList.DistinctEpochs(epochs, increasing)" << std::endl;
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			failed = true;
        }
        if (epochs.size() != 3)
        {
			std::cout << "Failed at bodyGroupList.DistinctEpochs(epochs, increasing)" << std::endl;
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			failed = true;
        }
    }

    // GetEpoch()
    {
        BodyGroupList bodyGroupList;

        int result = bodyGroupList.GetEpoch(value, First);
        if (result != -1)
        {
			std::cout << "Failed at bodyGroupList.GetEpoch(value, First)" << std::endl;
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			failed = true;
        }
        
        result = bodyGroupList.GetEpoch(value, Last);
        if (result != -1)
        {
			std::cout << "Failed at bodyGroupList.GetEpoch(value, Last)" << std::endl;
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			failed = true;
        }

        BodyGroup bg0;
        bg0.epoch = "19900101";
        BodyGroup bg1;
        bg1.epoch = "19910101";
        BodyGroup bg2;
        bg2.epoch = "19920101";
        BodyGroup bg3;
        bg3.epoch = "19930101";

        bodyGroupList.items.push_back(bg0);
        bodyGroupList.items.push_back(bg1);
        bodyGroupList.items.push_back(bg2);
        bodyGroupList.items.push_back(bg3);

        result = bodyGroupList.GetEpoch(value, First);
        if (result != 0 && value != 2447892.5)
        {
			std::cout << "Failed at bodyGroupList.GetEpoch(value, First)" << std::endl;
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			failed = true;
        }

        result = bodyGroupList.GetEpoch(value, Last);
        if (result != 0 && value != 2448999.5)
        {
			std::cout << "Failed at bodyGroupList.GetEpoch(value, Last)" << std::endl;
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			failed = true;
        }
    }

	// SetStartTime()
	{
        BodyGroupList bodyGroupList;

		int result = bodyGroupList.SetStartTime(0.0);
		if (result != 0)
        {
			std::cout << "Failed at bodyGroupList.SetStartTime(0.0)" << std::endl;
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			failed = true;
        }

        BodyGroup bg0;
        bg0.epoch = "19900101";
        BodyGroup bg1;
        bg1.epoch = "19910101";
        BodyGroup bg2;
        bg2.epoch = "19920101";
        BodyGroup bg3;

        bodyGroupList.items.push_back(bg0);
        bodyGroupList.items.push_back(bg1);
        bodyGroupList.items.push_back(bg2);
        bodyGroupList.items.push_back(bg3);

		result = bodyGroupList.SetStartTime(100.0);
		if (result != 0)
        {
			std::cout << "Failed at bodyGroupList.SetStartTime(100.0)" << std::endl;
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			failed = true;
        }
		std::list<BodyGroup>::iterator it = bodyGroupList.items.begin();
		if (it->startTime != 2447892.5)
        {
			std::cout << "Failed at it->startTime != 2447892.5" << std::endl;
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			failed = true;
        }
		it++;
		if (it->startTime != 2448257.5)
        {
			std::cout << "Failed at it->startTime != 2448257.5" << std::endl;
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			failed = true;
        }
		it++;
		if (it->startTime != 2448622.5)
        {
			std::cout << "Failed at it->startTime != 2448622.5" << std::endl;
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			failed = true;
        }
		it++;
		if (it->startTime != 100.0)
        {
			std::cout << "Failed at it->startTime != 100.0" << std::endl;
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			failed = true;
        }
	}

    // DistinctStartTimes()
    {
        BodyGroupList bodyGroupList;

        std::list<double> startTimes;
        bool increasing = true;

		bodyGroupList.DistinctStartTimes(startTimes, increasing);
		if (startTimes.size() != 0)
        {
			std::cout << "Failed at bodyGroupList.DistinctStartTimes(startTimes, increasing)" << std::endl;
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			failed = true;
        }

        BodyGroup bg0;
        bg0.epoch = "19900101";
        BodyGroup bg1;
        BodyGroup bg2;
        bg2.epoch = "19930101";
        BodyGroup bg3;
        bg3.epoch = "19930101";

        bodyGroupList.items.push_back(bg0);
        bodyGroupList.items.push_back(bg1);
        bodyGroupList.items.push_back(bg2);
        bodyGroupList.items.push_back(bg3);

		bodyGroupList.SetStartTime(100.0);
		bodyGroupList.DistinctStartTimes(startTimes, increasing);

		if (startTimes.size() != 3)
        {
			std::cout << "Failed at bodyGroupList.DistinctStartTimes(startTimes, increasing)" << std::endl;
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			failed = true;
        }
    }

	// GetStartTime()
    {
		double startTimeOfMainIntegrationPhase = 100.0;
        BodyGroupList bodyGroupList;

        BodyGroup bg0;
        bg0.epoch = "19900101";
        BodyGroup bg1;
        BodyGroup bg2;
        bg2.epoch = "19930101";
        BodyGroup bg3;
        bg3.epoch = "19930101";

        bodyGroupList.items.push_back(bg0);
        bodyGroupList.items.push_back(bg1);
        bodyGroupList.items.push_back(bg2);
        bodyGroupList.items.push_back(bg3);

		int result = bodyGroupList.SetStartTime(startTimeOfMainIntegrationPhase);
		if (result != 0)
        {
			std::cout << "Failed at bodyGroupList.SetStartTime(startTimeOfMainIntegrationPhase)" << std::endl;
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			failed = true;
        }

		double startTime = bodyGroupList.GetStartTime(First);
		if (startTime != 100.0)
        {
			std::cout << "Failed at bodyGroupList.GetStartTime(First)" << std::endl;
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			failed = true;
        }

		startTime = bodyGroupList.GetStartTime(Last);
		if (startTime != 2448988.5)
        {
			std::cout << "Failed at bodyGroupList.GetStartTime(Last)" << std::endl;
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			failed = true;
        }
    }

	// DistinctStartTimesOfMassiveBodies()
	{
		double startTimeOfMainIntegrationPhase = 100.0;
        BodyGroupList bodyGroupList;

        BodyGroup bg0;
        bg0.epoch = "19900101";
        BodyGroup bg1;
        BodyGroup bg2;
        bg2.epoch = "19930101";

		Body star(CentralBody);
		star.characteristics = new Characteristics(1.0);

		Body planet(RockyPlanet);
		planet.characteristics = new Characteristics(1.0 * Constants::EarthToSolar);

		Body testP(TestParticle);

		bg0.items.push_back(star);
		bg0.items.push_back(planet);

		bg1.items.push_back(testP);
		bg1.items.push_back(testP);

		bg2.items.push_back(testP);
		bg2.items.push_back(testP);

        bodyGroupList.items.push_back(bg0);
        bodyGroupList.items.push_back(bg1);
        bodyGroupList.items.push_back(bg2);

		int result = bodyGroupList.SetStartTime(startTimeOfMainIntegrationPhase);
		if (result != 0)
        {
			std::cout << "Failed at bodyGroupList.SetStartTime(startTimeOfMainIntegrationPhase)" << std::endl;
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			failed = true;
        }

		int nBG;
		result = bodyGroupList.DistinctStartTimesOfMassiveBodies(nBG);
		if (result != 0)
        {
			std::cout << "Failed at bodyGroupList.DistinctStartTimesOfMassiveBodies(nBG)" << std::endl;
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			failed = true;
        }
		if (nBG != 1)
        {
			std::cout << "Failed at bodyGroupList.DistinctStartTimesOfMassiveBodies(nBG)" << std::endl;
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			failed = true;
        }

		bg0.items.clear();
		bg1.items.clear();
		bg2.items.clear();
		bodyGroupList.items.clear();

		bg0.items.push_back(star);
		bg0.items.push_back(planet);

		bg1.items.push_back(testP);
		bg1.items.push_back(planet);

		bg2.items.push_back(testP);
		bg2.items.push_back(testP);

		bodyGroupList.items.push_back(bg0);
        bodyGroupList.items.push_back(bg1);
        bodyGroupList.items.push_back(bg2);

		result = bodyGroupList.SetStartTime(startTimeOfMainIntegrationPhase);
		if (result != 0)
        {
			std::cout << "Failed at bodyGroupList.SetStartTime(startTimeOfMainIntegrationPhase)" << std::endl;
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			failed = true;
        }
		result = bodyGroupList.DistinctStartTimesOfMassiveBodies(nBG);
		if (result != 0)
        {
			std::cout << "Failed at bodyGroupList.DistinctStartTimesOfMassiveBodies(nBG)" << std::endl;
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			failed = true;
        }
		if (nBG != 2)
        {
			std::cout << "Failed at bodyGroupList.DistinctStartTimesOfMassiveBodies(nBG)" << std::endl;
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			failed = true;
        }
	}

	// CountBy()
	{
        BodyGroupList bodyGroupList;

        BodyGroup bg0;
        bg0.epoch = "19900101";
        BodyGroup bg1;
        BodyGroup bg2;
        bg2.epoch = "19930101";

		Body star(CentralBody);
		star.characteristics = new Characteristics(1.0);

		Body proto(ProtoPlanet);
		proto.characteristics = new Characteristics(1.0 * Constants::EarthMoonToSolar);
	
		Body rockyPlanet(RockyPlanet);
		rockyPlanet.characteristics = new Characteristics(1.0 * Constants::EarthToSolar);
	
		Body giant(GiantPlanet);
		giant.characteristics = new Characteristics(1.0 * Constants::JupiterToSolar);
	
		Body superPl(SuperPlanetesimal);
		superPl.characteristics = new Characteristics(1.0e3 * Constants::KilogramToSolar);
	
		Body Pl(Planetesimal);
		Pl.characteristics = new Characteristics(1.0 * Constants::KilogramToSolar);
	
		Body testP(TestParticle);

		bg0.items.push_back(star);

		bg0.items.push_back(proto);
		bg0.items.push_back(proto);

		bg0.items.push_back(rockyPlanet);
		bg0.items.push_back(rockyPlanet);
		bg0.items.push_back(rockyPlanet);

		bg2.items.push_back(giant);
		bg2.items.push_back(giant);
		bg2.items.push_back(giant);
		bg2.items.push_back(giant);

		bg1.items.push_back(superPl);
		bg1.items.push_back(superPl);
		bg1.items.push_back(superPl);
		bg1.items.push_back(superPl);
		bg1.items.push_back(superPl);

		bg1.items.push_back(Pl);
		bg1.items.push_back(Pl);
		bg1.items.push_back(Pl);
		bg1.items.push_back(Pl);
		bg1.items.push_back(Pl);
		bg1.items.push_back(Pl);

		bg1.items.push_back(testP);
		bg1.items.push_back(testP);
		bg1.items.push_back(testP);
		bg1.items.push_back(testP);
		bg1.items.push_back(testP);
		bg1.items.push_back(testP);
		bg1.items.push_back(testP);


		bodyGroupList.items.push_back(bg0);
        bodyGroupList.items.push_back(bg1);
        bodyGroupList.items.push_back(bg2);

		int nBody = bodyGroupList.CountBy(CentralBody);
		if (nBody != 1)
        {
			std::cout << "Failed at bodyGroupList.CountBy(CentralBody)" << std::endl;
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			failed = true;
        }
		nBody = bodyGroupList.CountBy(ProtoPlanet);
		if (nBody != 2)
        {
			std::cout << "Failed at bodyGroupList.CountBy(ProtoPlanet)" << std::endl;
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			failed = true;
        }
		nBody = bodyGroupList.CountBy(RockyPlanet);
		if (nBody != 3)
        {
			std::cout << "Failed at bodyGroupList.CountBy(RockyPlanet)" << std::endl;
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			failed = true;
        }
		nBody = bodyGroupList.CountBy(GiantPlanet);
		if (nBody != 4)
        {
			std::cout << "Failed at bodyGroupList.CountBy(GiantPlanet)" << std::endl;
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			failed = true;
        }
		nBody = bodyGroupList.CountBy(SuperPlanetesimal);
		if (nBody != 5)
        {
			std::cout << "Failed at bodyGroupList.CountBy(SuperPlanetesimal)" << std::endl;
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			failed = true;
        }
		nBody = bodyGroupList.CountBy(Planetesimal);
		if (nBody != 6)
        {
			std::cout << "Failed at bodyGroupList.CountBy(Planetesimal)" << std::endl;
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			failed = true;
        }
		nBody = bodyGroupList.CountBy(TestParticle);
		if (nBody != 7)
        {
			std::cout << "Failed at bodyGroupList.CountBy(TestParticle)" << std::endl;
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			failed = true;
        }
	}

	// FindBy(double startTime, std::list<BodyGroup *> &result)
	{
		double startTimeOfMainIntegrationPhase = 100.0;
        BodyGroupList bodyGroupList;

        BodyGroup bg0("bg0");
        bg0.epoch = "19900101";
        BodyGroup bg1("bg1");
        BodyGroup bg2("bg2");
        bg2.epoch = "19930101";
        BodyGroup bg3("bg3");
        bg3.epoch = "19900101";

		bodyGroupList.items.push_back(bg0);
		bodyGroupList.items.push_back(bg1);
		bodyGroupList.items.push_back(bg2);
		bodyGroupList.items.push_back(bg3);

		bodyGroupList.SetStartTime(startTimeOfMainIntegrationPhase);

		std::list<BodyGroup *> result;
		bodyGroupList.FindBy(100.0, result);
		std::list<BodyGroup *>::iterator it = result.begin();
		if ((*it)->guid != "bg1") 
        {
			std::cout << "Failed at bodyGroupList.FindBy(100.0, result)" << std::endl;
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			failed = true;
        }

		result.clear();
		bodyGroupList.FindBy(2447892.5, result);
		it = result.begin();
		if ((*it)->guid != "bg0") 
        {
			std::cout << "Failed at bodyGroupList.FindBy(2447892.5, result)" << std::endl;
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			failed = true;
        }
		it++;
		if ((*it)->guid != "bg3") 
        {
			std::cout << "Failed at (*it)->guid != \"bg3\"" << std::endl;
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			failed = true;
        }
	}

	// FindBy(BodyType type, std::list<Body *> &result)
	{
        BodyGroupList bodyGroupList;

        BodyGroup bg0;
        bg0.epoch = "19900101";
        BodyGroup bg1;
        BodyGroup bg2;
        bg2.epoch = "19930101";

		Body star(CentralBody);
		Body proto(ProtoPlanet);
		Body rockyPlanet(RockyPlanet);
		Body giant(GiantPlanet);
		Body superPl(SuperPlanetesimal);
		Body Pl(Planetesimal);
		Body testP(TestParticle);

		bg0.items.push_back(star);

		bg0.items.push_back(proto);
		bg0.items.push_back(proto);

		bg0.items.push_back(rockyPlanet);
		bg0.items.push_back(rockyPlanet);
		bg0.items.push_back(rockyPlanet);

		bg2.items.push_back(giant);
		bg2.items.push_back(giant);
		bg2.items.push_back(giant);
		bg2.items.push_back(giant);

		bg1.items.push_back(superPl);
		bg1.items.push_back(superPl);
		bg1.items.push_back(superPl);
		bg1.items.push_back(superPl);
		bg1.items.push_back(superPl);

		bg1.items.push_back(Pl);
		bg1.items.push_back(Pl);
		bg1.items.push_back(Pl);
		bg1.items.push_back(Pl);
		bg1.items.push_back(Pl);
		bg1.items.push_back(Pl);

		bg1.items.push_back(testP);
		bg1.items.push_back(testP);
		bg1.items.push_back(testP);
		bg1.items.push_back(testP);
		bg1.items.push_back(testP);
		bg1.items.push_back(testP);
		bg1.items.push_back(testP);

		bodyGroupList.items.push_back(bg0);
        bodyGroupList.items.push_back(bg1);
        bodyGroupList.items.push_back(bg2);

		std::list<Body *> result;
		bodyGroupList.FindBy(CentralBody, result);
		if (result.size() != 1)
        {
			std::cout << "Failed at bodyGroupList.FindBy(CentralBody, result)" << std::endl;
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			failed = true;
        }
		std::list<Body *>::iterator it = result.begin();
		for ( ; it != result.end(); it++) 
		{
			if ((*it)->type != CentralBody)
			{
				std::cout << "Failed at (*it)->type != CentralBody" << std::endl;
				Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
				failed = true;
			}
		}
		result.clear();

		bodyGroupList.FindBy(ProtoPlanet, result);
		if (result.size() != 2)
        {
			std::cout << "Failed at bodyGroupList.FindBy(ProtoPlanet, result)" << std::endl;
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			failed = true;
        }
		it = result.begin();
		for ( ; it != result.end(); it++) 
		{
			if ((*it)->type != ProtoPlanet)
			{
				std::cout << "Failed at bodyGroupList.FindBy(ProtoPlanet, result)" << std::endl;
				Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
				failed = true;
			}
		}
		result.clear();
	
		bodyGroupList.FindBy(GiantPlanet, result);
		if (result.size() != 4)
        {
			std::cout << "Failed at bodyGroupList.FindBy(GiantPlanet, result)" << std::endl;
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			failed = true;
        }
		it = result.begin();
		for ( ; it != result.end(); it++) 
		{
			if ((*it)->type != GiantPlanet)
			{
				std::cout << "Failed at bodyGroupList.FindBy(GiantPlanet, result)" << std::endl;
				Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
				failed = true;
			}
		}
	}

	// GetBodyGroupWithMassiveBodies()
	{
        BodyGroupList bodyGroupList;

        BodyGroup bg0("bg0");
        bg0.epoch = "19900101";
        BodyGroup bg1("bg1");
        BodyGroup bg2("bg2");
        bg2.epoch = "19930101";

		Body star(CentralBody);
		star.characteristics = new Characteristics(1.0);
		Body proto(ProtoPlanet);
		proto.characteristics = new Characteristics(1.0);
		Body testP(TestParticle);

		bg0.items.push_back(star);
		bg0.items.push_back(proto);
		bg0.items.push_back(proto);

		bg1.items.push_back(testP);
		bg1.items.push_back(testP);
		bg1.items.push_back(testP);

		bg2.items.push_back(testP);
		bg2.items.push_back(testP);
		bg2.items.push_back(testP);

		bodyGroupList.items.push_back(bg0);
		bodyGroupList.items.push_back(bg1);
		bodyGroupList.items.push_back(bg2);

		std::list<BodyGroup>::iterator it;
		bool contains = bodyGroupList.GetBodyGroupWithMassiveBodies(it);
		if (!contains)
		{
			std::cout << "Failed at bodyGroupList.GetBodyGroupWithMassiveBodies(it)" << std::endl;
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			failed = true;
		}
		if (it->guid != "bg0")
		{
			std::cout << "Failed at it->guid != \"bg0\"" << std::endl;
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			failed = true;
		}
	}

	return failed;
}

bool TestBodyGroup()
{
	bool	failed = false;
	double	value = 0.0;

    // 	int		CountBy(BodyType type);
    {
		BodyGroup bodyGroup;

		Body star(CentralBody);
		star.characteristics = new Characteristics(1.0);

		Body giant(GiantPlanet);
		giant.characteristics = new Characteristics(1.0 * Constants::JupiterToSolar);

		Body rockyPlanet(RockyPlanet);
		rockyPlanet.characteristics = new Characteristics(1.0 * Constants::EarthToSolar);
	
		Body proto(ProtoPlanet);
		proto.characteristics = new Characteristics(1.0 * Constants::EarthMoonToSolar);
	
		Body superPl(SuperPlanetesimal);
		superPl.characteristics = new Characteristics(1.0e3 * Constants::KilogramToSolar);
	
		Body Pl(Planetesimal);
		Pl.characteristics = new Characteristics(1.0 * Constants::KilogramToSolar);
	
		Body testP(TestParticle);

		bodyGroup.items.push_back(star);
		bodyGroup.items.push_back(star);
		bodyGroup.items.push_back(star);

		int count = bodyGroup.CountBy(CentralBody);
		if (3 != count)
		{
			std::cout << "Failed at bodyGroup.CountBy(CentralBody)" << std::endl;
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			failed = true;
		}

		bodyGroup.items.push_back(giant);
		count = bodyGroup.CountBy(GiantPlanet);
		if (1 != count)
		{
			std::cout << "Failed at bodyGroup.CountBy(GiantPlanet)" << std::endl;
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			failed = true;
		}

		count = bodyGroup.CountBy(RockyPlanet);
		if (0 != count)
		{
			std::cout << "Failed at bodyGroup.CountBy(RockyPlanet)" << std::endl;
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			failed = true;
		}
	}

	// int		CountBy(double mass);
	{
		BodyGroup bodyGroup;

		Body star(CentralBody);
		star.characteristics = new Characteristics(1.0);

		Body giant(GiantPlanet);
		giant.characteristics = new Characteristics(1.0 * Constants::JupiterToSolar);

		Body rockyPlanet(RockyPlanet);
		rockyPlanet.characteristics = new Characteristics(1.0 * Constants::EarthToSolar);
	
		Body proto(ProtoPlanet);
		proto.characteristics = new Characteristics(1.0 * Constants::EarthMoonToSolar);
	
		Body superPl(SuperPlanetesimal);
		superPl.characteristics = new Characteristics(1.0e3 * Constants::KilogramToSolar);
	
		Body Pl(Planetesimal);
		Pl.characteristics = new Characteristics(1.0 * Constants::KilogramToSolar);
	
		Body testP(TestParticle);

		bodyGroup.items.push_back(star);
		bodyGroup.items.push_back(star);
		bodyGroup.items.push_back(star);

		int count = bodyGroup.CountBy(1.0);
		if (3 != count)
		{
			std::cout << "Failed at bodyGroup.CountBy(1.0)" << std::endl;
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			failed = true;
		}

		bodyGroup.items.clear();
		bodyGroup.items.push_back(giant);
		bodyGroup.items.push_back(giant);
		count = bodyGroup.CountBy(1.0 * Constants::JupiterToSolar);
		if (2 != count)
		{
			std::cout << "Failed at bodyGroup.CountBy(1.0 * Constants::JupiterToSolar)" << std::endl;
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			failed = true;
		}
	}

	// 	void	FindBy(BodyType type, std::list<Body *> &result);
	{
		BodyGroup bodyGroup;

		Body star(CentralBody);
		star.characteristics = new Characteristics(1.0);

		Body giant(GiantPlanet);
		giant.characteristics = new Characteristics(1.0 * Constants::JupiterToSolar);

		Body rockyPlanet(RockyPlanet);
		rockyPlanet.characteristics = new Characteristics(1.0 * Constants::EarthToSolar);
	
		Body proto(ProtoPlanet);
		proto.characteristics = new Characteristics(1.0 * Constants::EarthMoonToSolar);
	
		Body superPl(SuperPlanetesimal);
		superPl.characteristics = new Characteristics(1.0e3 * Constants::KilogramToSolar);
	
		Body Pl(Planetesimal);
		Pl.characteristics = new Characteristics(1.0 * Constants::KilogramToSolar);
	
		Body testP(TestParticle);

		bodyGroup.items.push_back(star);
		bodyGroup.items.push_back(star);
		bodyGroup.items.push_back(star);
		bodyGroup.items.push_back(rockyPlanet);
		bodyGroup.items.push_back(superPl);

		std::list<Body *> result;
		bodyGroup.FindBy(CentralBody, result);
		if (3 != result.size())
		{
			std::cout << "Failed at bodyGroup.FindBy(CentralBody, result)" << std::endl;
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			failed = true;
		}
		std::list<Body *>::iterator it = result.begin();
		for ( ; it != result.end(); it++) 
		{
			if ((*it)->type != CentralBody)
			{
				std::cout << "Failed at (*it)->type != CentralBody" << std::endl;
				Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
				failed = true;
			}
		}
		result.clear();

		bodyGroup.FindBy(RockyPlanet, result);
		if (1 != result.size())
		{
			std::cout << "Failed at bodyGroup.FindBy(CentralBody, result)" << std::endl;
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			failed = true;
		}
		it = result.begin();
		for ( ; it != result.end(); it++) 
		{
			if ((*it)->type != RockyPlanet)
			{
				std::cout << "Failed at (*it)->type != RockyPlanet" << std::endl;
				Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
				failed = true;
			}
		}
		result.clear();
	}

	// int		ToBodyList(std::list<Body *> &result)
	{
		BodyGroup bodyGroup;

		Body star(CentralBody);
		star.characteristics = new Characteristics(1.0);

		Body giant(GiantPlanet);
		giant.characteristics = new Characteristics(1.0 * Constants::JupiterToSolar);

		Body rockyPlanet(RockyPlanet);
		rockyPlanet.characteristics = new Characteristics(1.0 * Constants::EarthToSolar);
	
		Body proto(ProtoPlanet);
		proto.characteristics = new Characteristics(1.0 * Constants::EarthMoonToSolar);
	
		Body superPl(SuperPlanetesimal);
		superPl.characteristics = new Characteristics(1.0e3 * Constants::KilogramToSolar);
	
		Body Pl(Planetesimal);
		Pl.characteristics = new Characteristics(1.0 * Constants::KilogramToSolar);
	
		Body testP(TestParticle);

		bodyGroup.items.push_back(star);
		bodyGroup.items.push_back(star);
		bodyGroup.items.push_back(star);
		bodyGroup.items.push_back(rockyPlanet);
		bodyGroup.items.push_back(superPl);

		std::list<Body *> result;
		bodyGroup.ToBodyList(result);
		if (5 != result.size())
		{
			std::cout << "Failed at bodyGroup.ToBodyList(result)" << std::endl;
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			failed = true;
		}
		std::list<Body *>::iterator it = result.begin();
		std::list<Body>::iterator _it = bodyGroup.items.begin();
		for ( ; it != result.end() && _it != bodyGroup.items.end(); it++, _it++) 
		{
			if ((*it)->type != _it->type)
			{
				std::cout << "Failed at (*it)->type != _it->type" << std::endl;
				Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
				failed = true;
			}
		}
	}
	
	// int		SetStartTime(double startTimeOfMainSimulation)
	{
		double startTimeOfMainSimulation = 100.0;

		BodyGroup bodyGroup0("89A79746-3D47-49A0-9026-157F61420376", "Test body group 0", "",          0.0, "J2000.0");
		BodyGroup bodyGroup1("8C4C23A0-D162-419D-95EB-46F5034572A1", "Test body group 1", "20010101",  0.0, "J2000.0");
		BodyGroup bodyGroup2("A130417E-D47C-43F7-B6D9-A1FD8CCCAAB0", "Test body group 2", "20020101", 10.0, "J2000.0");
		BodyGroup bodyGroup3("6CDFBA3F-8AFA-485F-8507-17C8C9E673A7", "Test body group 3", "",         20.0, "J2000.0");

		bodyGroup0.SetStartTime(startTimeOfMainSimulation);
		if (startTimeOfMainSimulation != bodyGroup0.startTime)
		{
			std::cout << "Failed at bodyGroup0.SetStartTime(startTimeOfMainSimulation)" << std::endl;
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			failed = true;
		}

		bodyGroup1.SetStartTime(startTimeOfMainSimulation);
		if (2451910.5 != bodyGroup1.startTime)
		{
			std::cout << "Failed at bodyGroup1.SetStartTime(startTimeOfMainSimulation)" << std::endl;
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			failed = true;
		}

		bodyGroup2.SetStartTime(startTimeOfMainSimulation);
		if (2452275.5 + 10.0 != bodyGroup2.startTime)
		{
			std::cout << "Failed at bodyGroup2.SetStartTime(startTimeOfMainSimulation)" << std::endl;
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			failed = true;
		}

		bodyGroup3.SetStartTime(startTimeOfMainSimulation);
		if (startTimeOfMainSimulation + 20.0 != bodyGroup3.startTime)
		{
			std::cout << "Failed at bodyGroup3.SetStartTime(startTimeOfMainSimulation)" << std::endl;
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			failed = true;
		}
	}

	// bool	ContainsMassiveBody()
	{
		BodyGroup bodyGroup;

		Body star(CentralBody);
		star.characteristics = new Characteristics(1.0);

		Body giant(GiantPlanet);
		giant.characteristics = new Characteristics(1.0 * Constants::JupiterToSolar);

		Body rockyPlanet(RockyPlanet);
		rockyPlanet.characteristics = new Characteristics(1.0 * Constants::EarthToSolar);
	
		Body proto(ProtoPlanet);
		proto.characteristics = new Characteristics(1.0 * Constants::EarthMoonToSolar);
	
		Body superPl(SuperPlanetesimal);
		superPl.characteristics = new Characteristics(1.0e3 * Constants::KilogramToSolar);
	
		Body Pl(Planetesimal);
		Pl.characteristics = new Characteristics(1.0 * Constants::KilogramToSolar);
	
		Body testP(TestParticle);

		bodyGroup.items.push_back(testP);
		bodyGroup.items.push_back(testP);
		bodyGroup.items.push_back(testP);
		bodyGroup.items.push_back(testP);

		bool result = bodyGroup.ContainsMassiveBody();
		if (true == result)
		{
			std::cout << "Failed at bodyGroup.ContainsMassiveBody()" << std::endl;
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			failed = true;
		}

		bodyGroup.items.push_back(proto);
		result = bodyGroup.ContainsMassiveBody();
		if (false == result)
		{
			std::cout << "Failed at bodyGroup.ContainsMassiveBody()" << std::endl;
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			failed = true;
		}
	}

	return failed;
}

bool TestTools()
{
	bool	failed = false;
	double	value = 0.0;

    // enum OS	GetOs()
    {
		OS os = Tools::GetOs();
		if (Windows != os)
		{
			std::cout << "Failed at Tools::GetOs()" << std::endl;
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			failed = true;
		}
	}

	// std::string Tools::GetWorkingDirectory()
	{
		std::string wd;
		int result = Tools::GetWorkingDirectory(wd);
		std::cout << "Current working directory: " << wd << std::endl;
	}

	// static std::string GetDirectory(const std::string path, const char directorySeparator)
	{
		std::string path = "E:\\Work\\VSSolutions\\solaris\\Debug\\Test.exe";
		char dirSep = '\\';

		std::string directory = Tools::GetDirectory(path, dirSep);
		if ("E:\\Work\\VSSolutions\\solaris\\Debug" != directory)
		{
			std::cout << "Failed at Tools::GetDirectory()" << std::endl;
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			failed = true;
		}

		path = "Test.exe";
		directory = Tools::GetDirectory(path, dirSep);
		if ("" != directory)
		{
			std::cout << "Failed at Tools::GetDirectory()" << std::endl;
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			failed = true;
		}
	}

	// static std::string GetFileName(const std::string path, const char directorySeparator)
	{
		std::string path = "E:\\Work\\VSSolutions\\solaris\\Debug\\Test.exe";
		char dirSep = '\\';

		std::string fileName = Tools::GetFileName(path, dirSep);
		if ("Test.exe" != fileName)
		{
			std::cout << "Failed at Tools::GetFileName()" << std::endl;
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			failed = true;
		}

		path = "Test.exe";
		fileName = Tools::GetFileName(path, dirSep);
		if ("Test.exe" != fileName)
		{
			std::cout << "Failed at Tools::GetFileName()" << std::endl;
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			failed = true;
		}
	}


	return failed;
}
/*

cd 'G:\Work\VSSolutions\Solaris\TestCases\UnifiedDragForce
set title "Drag laws for different regimes"
set xlabel "R [m]" 
set ylabel "a [m/s^2]"
set logscale
set grid
autom=1.495978707e11
dtos=86400.0
aud2toms2=autom/(dtos*dtos)
lambda=4.6e-014*autom
lambda1=lambda/10.0
lambda2=10.0*lambda
set arrow from lambda, 10 to lambda, 1e11 nohead lc rgb 'red'
set arrow from lambda1, 10 to lambda1, 1e11 nohead lc rgb 'black'
set arrow from lambda2, 10 to lambda2, 1e11 nohead lc rgb 'blue'
set out "DragForceWithTransition.jpeg"
set term jpeg large size 800,600
p "DragForceLin.txt" u ($1*autom):($2*aud2toms2) w l title "Epstein",  "DragForceLin.txt" u ($1*autom):($3*aud2toms2) w l title "Stokes", "DragForceLin.txt" u ($1*autom):($4*aud2toms2) w p title "Transition"
set term windows
*/

void CalculateGraphOfUnifiedDragLaw(const double mC, const double r, const double density, const double Cd, const double Dlambda, const char* path, const bool append)
{
	GasComponent gasComponent;
    // Calculate the mean free path at r
    double lambda   = gasComponent.meanFreePath.Evaluate(r);
    double lambda1  = lambda / Dlambda;
    double lambda2  = lambda * Dlambda;

    double radius0  = 1.0e-8 * Constants::MeterToAu;    // m 1.0e-8 m = 0.01 micro m
    double radius1  = 1.0 * Constants::MeterToAu;       // m
    double radius   = radius0;
    double dR       = radius0;

	// Calculate the relative velocity of the solid with respect to the gas
	double f = sqrt(1.0 - 2.0*gasComponent.eta.Evaluate(r));
	Vector vCircular(Constants::Gauss * sqrt( mC/r ), 0.0, 0.0);
	Vector vGas = f * vCircular;
    Vector vRelative = vCircular + (-1.0 * vGas);
    double v = vRelative.Length();

	// Calculate the mean thermal velocity of the gas molecules
	double vth = gasComponent.MeanThermalSpeed_CMU(mC, r);

    double gammaE = 1.0/(density * lambda1);
    double gammaS = 3.0/8.0 * Cd / (density * lambda2);

    double K = gammaS * v / (gammaE * vth);
    double eta = lambda2 / lambda1;         // = Dlambda^2
    double kappa = log10(K) / log10(eta);   // ha Dlambda = 10, akkor log10(eta) = 2
    double gamma = gammaE * vth * pow(lambda1, -kappa);;

    std::ofstream file(path, append);
    if (!file.is_open())
    {
        std::cout << "Unable to open file";
        return;
    }

    int step = 0;
    do {
        gammaE = 1.0/(density * radius);
        gammaS = 3.0/8.0 * Cd * gammaE;
        double fd = 0.0;
        if (radius > lambda1 && radius < lambda2)
        {
            fd = gamma * pow(radius, kappa);
        }

        // Write the data to the text file
        file << radius << " " << gammaE * vth << " " << gammaS * v << " " << fd << "\n";

        step++;
        if (step % 100 == 0) {
            dR *= 10.0;
        }
        radius += dR;
    } while (radius <= radius1);

    file.close();
}

/*
cd 'G:\Work\VSSolutions\Solaris\TestCases\UnifiedDragForce
set title "The mean free path of gas molecules"
set logscale
set grid
set xlabel "r [AU]"
set ylabel "Mean free path [m]"
set out "MeanFreePath.jpeg"
set term jpeg large size 800,600
p "MeanFreePath.txt" u 1:3 w l title "lambda"
set term windows
*/
void CalculateGraphOfMeanFreePath(const char* path, const bool append)
{
	GasComponent gasComponent;

    std::ofstream file(path, append);
    if (!file.is_open())
    {
        std::cout << "Unable to open file";
        return;
    }

    double r = 0.01; // AU
    do {
        // Calculate the mean free path at r
        double lambda   = gasComponent.meanFreePath.Evaluate(r);
        // Write the data to the text file
        file << r << " " << lambda << " " << lambda * Constants::AuToMeter << "\n";
        r += 0.1;  // AU
    } while (r <= 100.0);

    file.close();
}

/*
cd 'G:\Work\VSSolutions\Solaris\TestCases\UnifiedDragForce
set title "Acceleration caused by the Epstein and Stokes drag force"
set logscale
set grid
autom=1.495978707e11
dtos=86400.0
aud2toms2=autom/(dtos*dtos)
set xlabel "r [AU]"
set ylabel "a [m/s^2]"
set out "EpsteinStokesVersusDistanceR1e-6m.jpeg"
set term jpeg large size 800,600
p "EpsteinStokesVersusDistanceR1e-6.txt" u 1:($2*aud2toms2) w l title "Epstein", "EpsteinStokesVersusDistanceR1e-6.txt" u 1:($3*aud2toms2) w l title "Stokes"

p "EpsteinStokesVersusDistanceR1e-6.txt" u 1:($2*aud2toms2) w l title "Epstein", "EpsteinStokesVersusDistanceR1e-6.txt" u 1:($3*aud2toms2) w l title "Stokes", "EpsteinStokesVersusDistanceR1e-5.txt" u 1:($2*aud2toms2) w l title "Epstein", "EpsteinStokesVersusDistanceR1e-5.txt" u 1:($3*aud2toms2) w l title "Stokes"
set term windows
*/
void CalculateGraphOfEpsteinStokesVersusDistance(const double mC, const double density, const double R, const double Cd, const char* path, const bool append)
{
	GasComponent gasComponent;

    std::ofstream file(path, append);
    if (!file.is_open())
    {
        std::cout << "Unable to open file";
        return;
    }

    int step = 0;
    double r = 0.01;       // AU
    do {
	    // Calculate the relative velocity of the solid with respect to the gas
	    double f = sqrt(1.0 - 2.0*gasComponent.eta.Evaluate(r));
	    Vector vCircular(Constants::Gauss * sqrt( mC/r ), 0.0, 0.0);
	    Vector vGas = f * vCircular;
        Vector vRelative = vCircular + (-1.0 * vGas);
        double v = vRelative.Length();

	    // Calculate the mean thermal velocity of the gas molecules
	    double vth = gasComponent.MeanThermalSpeed_CMU(mC, r);

        double gammaE = 1.0/(density * R);
        double gammaS = 3.0/8.0 * Cd * gammaE;

        // Write the data to the text file
        file << r << " " << gammaE * vth << " " << gammaS * v << "\n";

        step++;
        r += 0.1;
    } while (r <= 100.0);

    file.close();
}

int		main(int argc, char* argv[], char* envp[])
{
	bool failed = false;

#ifdef TEST_DUSTPARTICLE
	int result = TestDustParticle(argc, argv);
#endif

#ifdef TEST_CALCULATE
	failed = TestCalculate();
	if (!failed) {
		std::cout << "TestCalculate() passed." << std::endl;
	}
	else {
		Error::PrintStackTrace();
	}
#endif

#ifdef TEST_CONSTANTS
	failed = TestConstants();
	if (!failed) {
		std::cout << "TestConstants() passed." << std::endl;
	}
	else {
		Error::PrintStackTrace();
	}
#endif

#ifdef TEST_TYPEI_MIGRATION_TIME
	failed = TestTypeIMigrationTime();
	if (!failed) {
		std::cout << "TestTypeIMigrationTime() passed." << std::endl;
	}
	else {
		Error::PrintStackTrace();
	}
#endif

#ifdef TEST_MEAN_FREE_PATH
    failed = TestMeanFreePath();
	if (!failed) {
		std::cout << "TestMeanFreePath() passed." << std::endl;
	}
	else {
		Error::PrintStackTrace();
	}
#endif

#ifdef TEST_BODYGROUPLIST
    failed = TestBodyGroupList();
	if (!failed) {
		std::cout << "TestBodyGroupList() passed." << std::endl;
	}
	else {
		Error::PrintStackTrace();
	}
#endif

#ifdef TEST_BODYGROUP
	failed = TestBodyGroup();
	if (!failed) {
		std::cout << "TestBodyGroup() passed." << std::endl;
	}
	else {
		Error::PrintStackTrace();
	}
#endif

#ifdef TEST_TOOLS
	failed = TestTools();
	if (!failed) {
		std::cout << "TestTools() passed." << std::endl;
	}
	else {
		Error::PrintStackTrace();
	}
#endif

#ifdef EPSTEIN_STOKES_VERSUS_DISTANCE
    const double mC = 1.0;      // mass of the star in solra mass unit
    const double r0 = 1.0;      // AU
    const double density = 1.0 * Constants::GramPerCm3ToSolarPerAu3;
    const double Cd = 1.0;      // drag coefficient
    //CalculateGraphOfUnifiedDragLaw(mC, r0, density, Cd, 10.0, "G:\\Work\\VSSolutions\\Solaris\\TestCases\\UnifiedDragForce\\DragForceLin.txt", false);

    //CalculateGraphOfMeanFreePath("G:\\Work\\VSSolutions\\Solaris\\TestCases\\UnifiedDragForce\\MeanFreePath.txt", false);

    double R = 1.0e-6 * Constants::MeterToAu;  // m
    CalculateGraphOfEpsteinStokesVersusDistance(mC, density, R, Cd, 
		"G:\\Work\\VSSolutions\\Solaris\\TestCases\\UnifiedDragForce\\EpsteinStokesVersusDistanceR1e-6.txt", false);
    
    R = 1.0e-5 * Constants::MeterToAu;  // m
    CalculateGraphOfEpsteinStokesVersusDistance(mC, density, R, Cd, 
		"G:\\Work\\VSSolutions\\Solaris\\TestCases\\UnifiedDragForce\\EpsteinStokesVersusDistanceR1e-5.txt", false);

    R = 1.0e-4 * Constants::MeterToAu;  // m
    CalculateGraphOfEpsteinStokesVersusDistance(mC, density, R, Cd, 
		"G:\\Work\\VSSolutions\\Solaris\\TestCases\\UnifiedDragForce\\EpsteinStokesVersusDistanceR1e-4.txt", false);
#endif

    if (failed)
		Error::PrintStackTrace();
	else
		std::cout << "Test passed." << std::endl;


	return 0;
}