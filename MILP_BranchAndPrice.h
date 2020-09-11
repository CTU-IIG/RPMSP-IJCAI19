#pragma once
#include "generic.h"
#include "MILPsolver.h"
#include "MILPsolver_primary.h"
#include "config_BranchAndPrice.h"
#include "config_MILPprimary.h"

#include <omp.h>
#include <atomic>

#define CNSTR_FORCE_PAIR 1
#define CNSTR_FORBID_PAIR 0

typedef vector<double> costs;
//typedef vector<double> pattern;
typedef vector<pattern> patterns;
typedef pair<uint16_t, uint16_t> member;
typedef pair<pattern, double> cPattern;
typedef pair<vector<pattern>, vector<double>> cPatterns;
typedef vector<vector<double>> matrix;

struct constraint {
	uint16_t first_member_ID;
	uint16_t second_member_ID;
	uint16_t relation;
};

struct iterationWrapper {
	vector<constraint>	applicableConstraints;
	vector<member>		memberDegree; // Member ID and count
	vector<uint16_t>	affiliationMap;
	uint16_t			nodeID;

	patterns	nodePatterns;
	costs		nodePatternCosts;
};

class MILP_BranchAndPrice : public MILPsolver {

	public:
		MILP_BranchAndPrice(config_BranchAndPrice c);
		~MILP_BranchAndPrice();

		void loadSolverConfig(config_BranchAndPrice c);
		Schedule solve();

	private:
		bool timed_out;
		stack<iterationWrapper> pendingIterations;
		config_BranchAndPrice conf;
		uint64_t iterationsPerformed;

		GRBModel* pricingModel;
		GRBVar*	  x; // Assignment variable in pricing model

		cPatterns runNextIteration(patterns* availablePatterns, costs* serviceLevels, double* objVal);
		cPatterns solveMasterProblem(vector<double> c, vector<pattern> patterns, vector<double>* jobPrices, double* gammaValue, double* objVal); // Master problem solver is regenerated with each new pattern
		cPattern  solvePricingProblem(); // Pricing problem model is generated only when branching occurs, otherwise prices are updated on the go
		cPattern  solvePricingProblem_enumeration(costs prices);
		cPattern  solvePricingProblem_heuristic(costs prices, double gamma);

		// Regression
		matrix			reg_features;
		vector<double>	reg_targets;
		vector<double>	reg_getFeatures(costs c);
		vector<double>	reg_train(vector<vector<double>> features, vector<double> targets);
		double			reg_getPrediction(vector<double> weights, vector<double> features);
		
		Schedule getScheduleFromMaster(cPatterns masterSolution);
		patterns getPatternsFromSchedule(Schedule s);
		costs	 getPricesFromSchedule(Schedule s);		

		int generateAdditionalPatterns(patterns* p, costs* c);
		void generatePricingModel();
		void updatePricesOnPricingModel(costs prices);		
		void checkPatternsConstraintViolation(patterns* p, costs* c);
		
		void interpretConstraint(constraint c, GRBVar* x, GRBModel* model);
		bool branchOnConstraints(patterns parentPatterns, costs patternCosts);
		
		cPatterns regeneratePatterns();

		vector<double>	 approxSteps_objFunction();
		vector<uint32_t> getFeasibleVariances();	

		// For version one of solver
		MILPsolver_primary* milp;

		// Diagnostics
		double totalRuntime_master;
		double totalRuntime_pricing;

		uint64_t totalPatternsGenerated;
		uint64_t heuristicPatternsGenerated;

		bool usePricingHeuristic;
		bool useRegression;
		ofstream timeOutput;
		ofstream graphFile;
		uint16_t nodeCounter;

		vector<uint16_t> jobsScheduledInPatterns;
};

