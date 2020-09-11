#include "config_MILPprimary.h"

config_MILPprimary::config_MILPprimary() {

}

config_MILPprimary::config_MILPprimary(string fileName, double objFunctionStep, double fracFunctionStep) {
	loadInstance(fileName);
	this->objFunctionStep	= objFunctionStep;
	this->fracFunctionStep	= fracFunctionStep;
}


config_MILPprimary::~config_MILPprimary() {

}
