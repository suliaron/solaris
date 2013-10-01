#include <iostream>
#include <fstream>

#include "Error.h"
#include "FargoParameters.h"
#include "Tokenizer.h"
#include "Tools.h"

FargoParameters::FargoParameters()
{
	aspectRatio			= 0.05;            // Thickness over Radius in the disc
	sigma0				= 1.45315e-5;      // Surface Density at r=1
	alphaViscosity		= 0.001;           // Uniform kinematic viscosity
	sigmaSlope          = 0.5;             // Slope of surface density profile.
	flaringIndex        = 0.0;
	excludeHill			= "YES";

// Planet parameters

	planetConfig        = "./planets.cfg";
	thicknessSmoothing  = 0.6;             // Smoothing parameters in disk thickness

// Numerical method parameters

	transport           = "FARGO";
	innerBoundary       = "STOCKHOLM";     // choose : OPEN or RIGID or NONREFLECTING
	disk                = "YES"; 
	omegaFrame          = 0;
	frame               = "COROTATING";
	indirectTerm        = "YES";

// Mesh parameters

	nRad                = 256;             // Radial number of zones
	nSec                = 512;             // Azimuthal number of zones (sectors)
	rMin                = 0.2;             // Inner boundary radius
	rMax                = 15.0;            // Outer boundary radius
	radialSpacing       = "Logarithmic";   // Zone interfaces evenly spaced

// Output control parameters

	nTot                = 800000;          // Total number of time steps
	nInterm             = 2000;            // Time steps between outputs
	dT                  = 0.314159;        // Time step length. 2PI = 1 orbit
	outputDir           = "./";

// Viscosity damping due to a dead zone

	viscModR1           = 0.0;             // Inner radius of dead zone
	viscModDeltaR1      = 0.0;             // Width of viscosity transition at inner radius
	viscModR2           = 0.0;             // Outer radius of dead zone
	viscModDeltaR2      = 0.0;             // Width of viscosity transition at outer radius
	viscMod             = 0.0;             // Viscosity damp
}

int FargoParameters::ReadConfigFile(std::string& path)
{
	std::ifstream file(path);
	if (file) {
		std::string str;
		while (std::getline(file, str))
		{
			if (str.length() == 0)
				continue;
			if (str[0] == '#')
				continue;
			config += str;
			config.push_back('\n');
		} 	
	}
	else {
		std::cerr << "The file '" << path << "' could not opened!\r\nExiting to system!" << std::endl;
		exit(1);
	}
	file.close();

	return 0;
}

int	FargoParameters::ParseConfig(const bool verbose)
{
	// instantiate Tokenizer classes
	Tokenizer fileTokenizer;
	Tokenizer lineTokenizer;
	std::string line;

	fileTokenizer.set(config, "\n");
	while ((line = fileTokenizer.next()) != "") {
		lineTokenizer.set(line, " \t");
		std::string token;
		int tokenCounter = 1;

		std::string key; 
		std::string value;
		while ((token = lineTokenizer.next()) != "" && tokenCounter <= 2) {

			if (tokenCounter == 1)
				key = token;
			else if (tokenCounter == 2)
				value = token;

			tokenCounter++;
		}
		if (tokenCounter > 2) {
			if (SetParameter(key, value, verbose) == 1)
				return 1;
		}
		else {
			Error::_errMsg = "Invalid key/value pair: " + line;
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			return 1;
		}
	}

	return 0;
}

int FargoParameters::SetParameter(std::string& key, std::string& value, const bool verbose)
{
	if (     key == "AspectRatio") {
		if (!Tools::IsNumber(value)) {
			Error::_errMsg = "Invalid number: '" + value + "'!";
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			return 1;
		}
		aspectRatio = atof(value.c_str());
    } 
    else if (key == "Sigma0") {
		if (!Tools::IsNumber(value)) {
			Error::_errMsg = "Invalid number: '" + value + "'!";
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			return 1;
		}
		sigma0 = atof(value.c_str());
    }
    else if (key == "AlphaViscosity") {
		if (!Tools::IsNumber(value)) {
			Error::_errMsg = "Invalid number: '" + value + "'!";
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			return 1;
		}
		alphaViscosity = atof(value.c_str());
    }
    else if (key == "SigmaSlope") {
		if (!Tools::IsNumber(value)) {
			Error::_errMsg = "Invalid number: '" + value + "'!";
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			return 1;
		}
		sigmaSlope = atof(value.c_str());
    }
    else if (key == "FlaringIndex") {
		if (!Tools::IsNumber(value)) {
			Error::_errMsg = "Invalid number: '" + value + "'!";
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			return 1;
		}
		flaringIndex = atof(value.c_str());
    }
    else if (key == "ExcludeHill") {
		excludeHill = value;
    }
    else if (key == "PlanetConfig") {
		planetConfig = value;
    }
    else if (key == "ThicknessSmoothing") {
		if (!Tools::IsNumber(value)) {
			Error::_errMsg = "Invalid number: '" + value + "'!";
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			return 1;
		}
		thicknessSmoothing = atof(value.c_str());
    }
    else if (key == "Transport") {
		transport = value;
    }
    else if (key == "InnerBoundary") {
		innerBoundary = value;
    }
    else if (key == "Disk") {
		disk = value;
    }
    else if (key == "OmegaFrame") {
		if (!Tools::IsNumber(value)) {
			Error::_errMsg = "Invalid number: '" + value + "'!";
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			return 1;
		}
		omegaFrame = atof(value.c_str());
    }
    else if (key == "Frame") {
		frame = value;
    }
    else if (key == "IndirectTerm") {
		indirectTerm = value;
    }
    else if (key == "Nrad") {
		if (!Tools::IsNumber(value)) {
			Error::_errMsg = "Invalid number: '" + value + "'!";
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			return 1;
		}
		nRad = atoi(value.c_str());
    }
    else if (key == "Nsec") {
		if (!Tools::IsNumber(value)) {
			Error::_errMsg = "Invalid number: '" + value + "'!";
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			return 1;
		}
		nSec = atoi(value.c_str());
    }
    else if (key == "Rmin") {
		if (!Tools::IsNumber(value)) {
			Error::_errMsg = "Invalid number: '" + value + "'!";
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			return 1;
		}
		rMin = atof(value.c_str());
    }
    else if (key == "Rmax") {
		if (!Tools::IsNumber(value)) {
			Error::_errMsg = "Invalid number: '" + value + "'!";
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			return 1;
		}
		rMax = atof(value.c_str());
    }
    else if (key == "RadialSpacing") {
		radialSpacing = value;
    }
    else if (key == "Ntot") {
		if (!Tools::IsNumber(value)) {
			Error::_errMsg = "Invalid number: '" + value + "'!";
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			return 1;
		}
		nTot = atoi(value.c_str());
    }
    else if (key == "Ninterm") {
		if (!Tools::IsNumber(value)) {
			Error::_errMsg = "Invalid number: '" + value + "'!";
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			return 1;
		}
		nInterm = atoi(value.c_str());
    }
    else if (key == "DT") {
		if (!Tools::IsNumber(value)) {
			Error::_errMsg = "Invalid number: '" + value + "'!";
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			return 1;
		}
		dT = atof(value.c_str());
    }
    else if (key == "OutputDir") {
		outputDir = value;
    }
    else if (key == "ViscModR1") {
		if (!Tools::IsNumber(value)) {
			Error::_errMsg = "Invalid number: '" + value + "'!";
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			return 1;
		}
		viscModR1 = atof(value.c_str());
    }
    else if (key == "ViscModDeltaR1") {
		if (!Tools::IsNumber(value)) {
			Error::_errMsg = "Invalid number: '" + value + "'!";
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			return 1;
		}
		viscModDeltaR1 = atof(value.c_str());
    }
    else if (key == "ViscModR2") {
		if (!Tools::IsNumber(value)) {
			Error::_errMsg = "Invalid number: '" + value + "'!";
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			return 1;
		}
		viscModR2 = atof(value.c_str());
    }
    else if (key == "ViscModDeltaR2") {
		if (!Tools::IsNumber(value)) {
			Error::_errMsg = "Invalid number: '" + value + "'!";
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			return 1;
		}
		viscModDeltaR2 = atof(value.c_str());
    }
    else if (key == "ViscMod") {
		if (!Tools::IsNumber(value)) {
			Error::_errMsg = "Invalid number: '" + value + "'!";
			Error::PushLocation(__FILE__, __FUNCTION__, __LINE__);
			return 1;
		}
		viscMod = atof(value.c_str());
    }
    else {
        std::cerr << "Unrecoginzed key: '" << key << "'!" << std::endl;
		return 1;
    }

	if (verbose) {
		std::cout << key << " was assigned to " << value << std::endl;
	}

	return 0;
}
