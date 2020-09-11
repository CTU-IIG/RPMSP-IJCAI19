#pragma once
#include "generic.h"
#include "Config.h"
#include "config_MILPprimary.h"

class config_BranchAndPrice : public Config {
	public:
		config_BranchAndPrice();
		config_BranchAndPrice(string fileName, config_MILPprimary pricingProblemConfig);
		~config_BranchAndPrice();

		config_MILPprimary pricingProblemConfig;
		uint16_t pricingProblem_solverVersion;

		double pricingProblem_criterionCoef; // Rename

		bool use_heuristic;
		bool use_regressionModel;

		double reg_cMinus;
		double reg_cPlus;
		double reg_epsilonTolerance;
		double reg_lambda_lb;
		double reg_lambda_ub;
		uint16_t reg_iterationsNeeded;
		uint16_t reg_type;

		int _allocate_threads;
};

