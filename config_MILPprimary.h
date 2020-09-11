#pragma once
#include "Config.h"

class config_MILPprimary :	public Config {
	
	public:
		config_MILPprimary();
		config_MILPprimary(string fileName, double objFunctionStep, double fracFunctionStep);
		~config_MILPprimary();

		double objFunctionStep;
		double fracFunctionStep;

		// Make this property of solver instead?
		uint16_t vMin;
		uint16_t vMax;
};

