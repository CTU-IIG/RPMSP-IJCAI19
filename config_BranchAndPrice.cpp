#include "config_BranchAndPrice.h"

config_BranchAndPrice::config_BranchAndPrice(){

}

config_BranchAndPrice::config_BranchAndPrice(string fileName, config_MILPprimary pricingProblemConfig) {
	loadInstance(fileName);

	this->pricingProblemConfig			= pricingProblemConfig;
	this->pricingProblem_solverVersion	= 1;
	this->pricingProblem_criterionCoef	= 0.000000005;

	this->use_heuristic			= false;
	this->use_regressionModel	= false;

	this->reg_cMinus			= 0.0;
	this->reg_cPlus				= 0.0;
	this->reg_epsilonTolerance	= 0.0;
	this->reg_lambda_lb			= 0.5;
	this->reg_lambda_ub			= 1;
	this->reg_iterationsNeeded	= 0;
	this->reg_type = 1;

	this->_allocate_threads = 0;
}


config_BranchAndPrice::~config_BranchAndPrice()
{
}
