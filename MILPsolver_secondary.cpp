#include "MILPsolver_secondary.h"


MILPsolver_secondary::MILPsolver_secondary(){

}


MILPsolver_secondary::~MILPsolver_secondary(){

}

void MILPsolver_secondary::setupSolver(instance i) {
	this->loadedInstance = i;
}

Schedule MILPsolver_secondary::solve() {
	uint16_t numOfMachines	= loadedInstance.numOfMachines;
	JobList  jobs			= loadedInstance.jobs;
	
	GRBEnv	 env;
	GRBModel model(env);
	Schedule result(numOfMachines);

	uint32_t N = jobs.size();
	double one = 1.0;

	GRBVar*** X = new GRBVar**[numOfMachines];
	for (uint32_t i = 0; i < numOfMachines; i++) {
		X[i] = new GRBVar*[N];;
		for (uint32_t j = 0; j < N; j++) {
			X[i][j] = new GRBVar[N];
			for (uint32_t k = 0; k < N; k++) {
				X[i][j][k] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
			}
		}
	}

	// Constraints
	for (uint32_t j = 0; j < N; j++) {
		GRBLinExpr* cnstr_onlyOnce = new GRBLinExpr();
		for (uint32_t i = 0; i < numOfMachines; i++) {
			for (uint16_t k = 0; k < N; k++) {
				cnstr_onlyOnce->addTerms(&one, &X[i][j][k], 1);
			}
		}

		model.addConstr(*cnstr_onlyOnce, GRB_EQUAL, 1.0);
	}

	for (uint32_t i = 0; i < numOfMachines; i++) {
		for (uint32_t k = 0; k < N; k++) {
			GRBLinExpr* cnstr_eachJob = new GRBLinExpr();
			for (uint32_t j = 0; j < N; j++) {
				cnstr_eachJob->addTerms(&one, &X[i][j][k], 1);
			}

			model.addConstr(*cnstr_eachJob, GRB_LESS_EQUAL, 1.0);
		}
	}

	// Criterion
	uint32_t minVar = UINT32_MAX;
	for (Job j : jobs) {
		if (j.second < minVar)
			minVar = j.second;
	}

	double coef_h = ((N / numOfMachines)*(N / numOfMachines + 1)*numOfMachines) / 2;
	double coef_l = 1 / (sqrt(coef_h*minVar));
	double coef_mean = 1.0;
	double coef_var = 1.0;
	GRBLinExpr* constraint = new GRBLinExpr();

	for (uint32_t i = 0; i < numOfMachines; i++) {
		for (uint32_t j = 0; j < N; j++) {
			for (uint32_t k = 0; k < N; k++) {
				coef_mean = k * jobs[j].first;
				constraint->addTerms(&coef_mean, &X[i][j][k], 1);

				coef_var = coef_l * (k ^ 2)*jobs[j].second;
				constraint->addTerms(&coef_var, &X[i][j][k], 1);
			}
		}
	}

	model.setObjective(*constraint, GRB_MAXIMIZE);
	model.optimize();

	// Results
	this->expandedNodes = (uint64_t)model.get(GRB_DoubleAttr_NodeCount);
	this->runtime = model.get(GRB_DoubleAttr_Runtime);
	double val;
	for (uint32_t i = 0; i < numOfMachines; i++) {
		for (uint32_t j = 0; j < N; j++) {
			for (uint32_t k = 0; k < N; k++) {
				val = X[i][j][k].get(GRB_DoubleAttr_X);

				if (val >= 0.5)
					result[i].push_back(jobs[j]);
			}
		}
	}

	return result;
}
