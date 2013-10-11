#include <cmath>
#include <iostream>
#include <fstream>

#include "../Solaris/Acceleration.h"
#include "../Solaris/BinaryFileAdapter.h"
#include "../Solaris/Body.h"
#include "../Solaris/BodyData.h"
#include "../Solaris/BodyGroup.h"
#include "../Solaris/BodyGroupList.h"
#include "../Solaris/Calculate.h"
#include "../Solaris/Characteristics.h"
#include "../Solaris/Component.h"
#include "../Solaris/Constants.h"
#include "../Solaris/DateTime.h"
#include "../Solaris/DormandPrince.h"
#include "../Solaris/DragCoefficient.h"
#include "../Solaris/Ephemeris.h"
#include "../Solaris/Error.h"
#include "../Solaris/Event.h"
#include "../Solaris/EventCondition.h"
#include "../Solaris/FargoParameters.h"
#include "../Solaris/GasComponent.h"
#include "../Solaris/GasDecreaseType.h"
#include "../Solaris/Integrator.h"
#include "../Solaris/IntegratorType.h"
#include "../Solaris/NBodies.h"
#include "../Solaris/Nebula.h"
#include "../Solaris/OrbitalElement.h"
#include "../Solaris/Output.h"
#include "../Solaris/Phase.h"
#include "../Solaris/PowerLaw.h"
#include "../Solaris/RungeKutta4.h"
#include "../Solaris/RungeKutta56.h"
#include "../Solaris/RungeKuttaFehlberg78.h"
#include "../Solaris/Settings.h"
#include "../Solaris/Simulation.h"
#include "../Solaris/Simulator.h"
#include "../Solaris/SolidsComponent.h"
#include "../Solaris/TimeLine.h"
#include "../Solaris/tinystr.h"
#include "../Solaris/tinyxml.h"
#include "../Solaris/Tokenizer.h"
#include "../Solaris/Tools.h"
#include "../Solaris/TwoBodyAffair.h"
#include "../Solaris/Units.h"
#include "../Solaris/Validator.h"
#include "../Solaris/Vector.h"
#include "../Solaris/XmlFileAdapter.h"


#define SQR(a)		((a)*(a))
#define CUBE(a)		((a)*(a)*(a))
#define FORTH(a)	((a)*(a)*(a)*(a))
#define FIFTH(a)	((a)*(a)*(a)*(a)*(a))

//#define EPSTEIN_STOKES_VERSUS_DISTANCE

//#define TEST_CALCULATE
//#define TEST_CONSTANTS
//#define TEST_TYPEI_MIGRATION_TIME
//#define TEST_MEAN_FREE_PATH
//#define TEST_BODYGROUPLIST

#define TEST_DUSTPARTICLE

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

    // DistinctEpochs()
    {
        BodyGroupList bodyGroupList;

        std::list<double> epochs;
        bool increasing = true;

        int result = bodyGroupList.DistinctEpochs(epochs, increasing);
        if (result != 0)
        {
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
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			failed = true;
        }
        if (epochs.size() != 3)
        {
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
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			failed = true;
        }
        
        result = bodyGroupList.GetEpoch(value, Last);
        if (result != -1)
        {
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
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			failed = true;
        }

        result = bodyGroupList.GetEpoch(value, Last);
        if (result != 0 && value != 2448999.5)
        {
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
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			failed = true;
        }
		std::list<BodyGroup>::iterator it = bodyGroupList.items.begin();
		if (it->startTime != 2447892.5)
        {
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			failed = true;
        }
		it++;
		if (it->startTime != 2448257.5)
        {
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			failed = true;
        }
		it++;
		if (it->startTime != 2448622.5)
        {
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			failed = true;
        }
		it++;
		if (it->startTime != 100.0)
        {
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
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			failed = true;
        }

		double startTime = bodyGroupList.GetStartTime(First);
		if (startTime != 100.0)
        {
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			failed = true;
        }

		startTime = bodyGroupList.GetStartTime(Last);
		if (startTime != 2448988.5)
        {
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
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			failed = true;
        }

		int nBG;
		result = bodyGroupList.DistinctStartTimesOfMassiveBodies(nBG);
		if (result != 0)
        {
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			failed = true;
        }
		if (nBG != 1)
        {
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
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			failed = true;
        }
		result = bodyGroupList.DistinctStartTimesOfMassiveBodies(nBG);
		if (result != 0)
        {
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			failed = true;
        }
		if (nBG != 2)
        {
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
		giant.characteristics = new Characteristics(1.0e3 * Constants::KilogramToSolar);
	
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
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			failed = true;
        }
		nBody = bodyGroupList.CountBy(ProtoPlanet);
		if (nBody != 2)
        {
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			failed = true;
        }
		nBody = bodyGroupList.CountBy(RockyPlanet);
		if (nBody != 3)
        {
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			failed = true;
        }
		nBody = bodyGroupList.CountBy(GiantPlanet);
		if (nBody != 4)
        {
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			failed = true;
        }
		nBody = bodyGroupList.CountBy(SuperPlanetesimal);
		if (nBody != 5)
        {
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			failed = true;
        }
		nBody = bodyGroupList.CountBy(Planetesimal);
		if (nBody != 6)
        {
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			failed = true;
        }
		nBody = bodyGroupList.CountBy(TestParticle);
		if (nBody != 7)
        {
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
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			failed = true;
        }

		result.clear();
		bodyGroupList.FindBy(2447892.5, result);
		it = result.begin();
		if ((*it)->guid != "bg0") 
        {
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			failed = true;
        }
		it++;
		if ((*it)->guid != "bg3") 
        {
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
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			failed = true;
        }
		std::list<Body *>::iterator it = result.begin();
		for ( ; it != result.end(); it++) 
		{
			if ((*it)->type != CentralBody)
			{
				Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
				failed = true;
			}
		}
		result.clear();

		bodyGroupList.FindBy(ProtoPlanet, result);
		if (result.size() != 2)
        {
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			failed = true;
        }
		it = result.begin();
		for ( ; it != result.end(); it++) 
		{
			if ((*it)->type != ProtoPlanet)
			{
				Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
				failed = true;
			}
		}
		result.clear();
	
		bodyGroupList.FindBy(GiantPlanet, result);
		if (result.size() != 4)
        {
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			failed = true;
        }
		it = result.begin();
		for ( ; it != result.end(); it++) 
		{
			if ((*it)->type != GiantPlanet)
			{
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
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			failed = true;
		}
		if (it->guid != "bg0")
		{
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

int ProcessArgv(int argc, char* argv[], std::string &directory, std::string &fileName, std::string &runType)
{
	if (argc < 2) { // Check the value of argc.
        Error::_errMsg = Constants::Usage;
        return 1;
    }
    for (int i = 1; i < argc; i++) {
        if (     strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0) {
			std::cout << Constants::CodeName << ":" << std::endl << Constants::Usage;
			exit(0);
		}
        else if (strcmp(argv[i], "-i") == 0) {
			i++;
			std::string input(argv[i]);			
			Tools::SplitPath(input, Output::directorySeparator, directory, fileName);
            runType = "New";
        } 
        else if (strcmp(argv[i], "-c") == 0) {
			i++;
			std::string input(argv[i]);			
			Tools::SplitPath(input, Output::directorySeparator, directory, fileName);
            runType = "Continue";
        } 
		else {
			Error::_errMsg = "Invalid argument.\n" + Constants::Usage;
            return 1;
		}
	}
	
	// If the directory is empty, then use the current/working directory
	if (directory.length() == 0) {
		directory = Tools::GetWorkingDirectory();
	}

	return 0;
}

/**
 * Loads the input data and stores it into the simulation object.
 *
 * @param inputPath the input path of the data file
 * @param simulation the object where the input data will be stored
 * @return 0 on success 1 on error
 */
int LoadInput(char* inputPath, Simulation &simulation)
{
	XmlFileAdapter xml(inputPath);
	if (XmlFileAdapter::Load(inputPath, xml.doc) == 1) {
		Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
		return 1;
	}

	if (XmlFileAdapter::DeserializeSimulation(xml.doc, simulation) == 1) {
		Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
		return 1;
	}

	if (simulation.nebula != 0 && simulation.nebula->path.length() > 0) {
		simulation.fargoParameters = new FargoParameters();
		simulation.fargoParameters->ReadConfigFile(simulation.nebula->path);
		if (simulation.fargoParameters->ParseConfig(true) == 1) {
    		Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			return 1;
        }
	}

	return 0;
}

int TestDustParticle(int argc, char* argv[])
{
	{
		if (Tools::GetDirectorySeparator(&Output::directorySeparator) == 1) {
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			return 1;
		}

		std::string		fileName;
		std::string     runType;
		if (ProcessArgv(argc, argv, Output::directory, fileName, runType) == 1) {
			Error::PrintStackTrace();
			exit(1);
		}

		char*	inputPath = 0;
		Tools::CreatePath(Output::directory, fileName, Output::directorySeparator, &inputPath);

		Simulation  simulation(runType);
		if (LoadInput(inputPath, simulation) == 1) {
			Error::PrintStackTrace();
			exit(1);
		}
	    simulation.binary = new BinaryFileAdapter(simulation.settings->output);
		simulation.binary->LogStartParameters(argc, argv);
		simulation.binary->Log("Simulation data was successfully loaded and deserialized.", true);

		if (simulation.Initialize() == 1) {
			Error::PrintStackTrace();
			exit(1);
		}
		simulation.binary->Log("Simulation was successfully initialized.", false);

		Simulator		simulator(&simulation);
		if (simulation.runType == "Continue" ) {
			if (simulator.Continue() == 1) {
				Error::PrintStackTrace();
				exit(1);
			}
		}

		if (simulator.Run() == 1) {
			Error::PrintStackTrace();
			exit(1);
		}

		return 0;
	}
}

int		main(int argc, char* argv[])
{
	bool failed = false;

#ifdef TEST_DUSTPARTICLE
	int result = TestDustParticle(argc, argv);
#endif

#ifdef TEST_CALCULATE
	failed = TestCalculate();
#endif

#ifdef TEST_CONSTANTS
	failed = TestConstants();
#endif

#ifdef TEST_TYPEI_MIGRATION_TIME
	failed = TestTypeIMigrationTime();
#endif

#ifdef TEST_MEAN_FREE_PATH
    failed = TestMeanFreePath();
#endif

#ifdef TEST_BODYGROUPLIST
    failed = TestBodyGroupList();
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