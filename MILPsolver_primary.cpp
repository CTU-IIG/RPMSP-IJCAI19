#include "MILPsolver_primary.h"

// Make instance and config one
MILPsolver_primary::MILPsolver_primary(config_MILPprimary c, vector<double> prices) {
	loadSolverConfig(c);

	if (prices.empty()) { 
		prices = vector<double>(loadedInstance.numOfJobs, 0.0); 
	} else {
		// Branch and price initiated, change number of machines to 1
		this->loadedInstance.numOfMachines = 1;
		this->conf.loadedInstance.numOfMachines = 1;
	}

	buildModel(prices);
	buildConstraints();
	buildObjective();
}

MILPsolver_primary::~MILPsolver_primary() {
	flush();
}

void MILPsolver_primary::flush() {
	delete[] v, u, t;

	for (uint16_t i = 0; i < loadedInstance.numOfMachines; i++) {
		delete[] z[i], w[i], s[i], y[i];
	}
	delete[] z, w, s, y;

	delete model, env;
}

void MILPsolver_primary::buildModel(vector<double> prices) {
	// Move this to protected function build model?
	uint16_t M = loadedInstance.numOfMachines;
	uint16_t N = loadedInstance.numOfJobs;
	uint16_t L = fracApproxIntervalBounds.size() - 1;
	JobList jobs = loadedInstance.jobs;

	// Make model
	env = new GRBEnv();
	model = new GRBModel(*env);

	// Model variables
	v = new GRBVar[M];
	u = new GRBVar[M];
	t = new GRBVar[M];
	z = new GRBVar*[M];
	w = new GRBVar*[M];
	s = new GRBVar*[M];
	y = new GRBVar*[M];

	// Variable init
	for (uint16_t i = 0; i < M; i++) {
		v[i] = model->addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS, "v" + to_string(i));
		u[i] = model->addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS, "u" + to_string(i));
		t[i] = model->addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS, "t" + to_string(i));
		z[i] = new GRBVar[N];
		w[i] = new GRBVar[N];
		s[i] = new GRBVar[L];
		y[i] = new GRBVar[L];

		for (uint16_t j = 0; j < N; j++) {
			z[i][j] = model->addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS, "z" + to_string(i) + to_string(j));
			w[i][j] = model->addVar(0.0, 1.0, prices[j], GRB_BINARY, "w" + to_string(i) + to_string(j)); // Objective value set to price
		}

		for (uint16_t k = 0; k < L; k++) {
			s[i][k] = model->addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS, "s" + to_string(i) + to_string(k));
			y[i][k] = model->addVar(0.0, 1.0, 0.0, GRB_BINARY, "y" + to_string(i) + to_string(k));
		}
	}
}


void MILPsolver_primary::buildConstraints() {
	// Constants initialization	
	uint16_t M = loadedInstance.numOfMachines;
	uint16_t N = loadedInstance.numOfJobs;
	JobList jobs = loadedInstance.jobs;

	double deadline = (double)loadedInstance.deadline;
	double pOne = 1.0;
	double nOne = -1.0;

	bigM_1 = 10000;
	bigM_2 = bigM_1;

		vector<double> h_bounds = fracApproxIntervalBounds;
	vector<double> h; // Calculate middle points of each interval
	for (uint32_t i = 1; i < h_bounds.size(); i++) { h.push_back((h_bounds[i] + h_bounds[i - 1]) / 2); }
	uint16_t L = h.size();

	// Transformed mean and var assignemnt
	double currentMean, currentVar;

	for (uint16_t i = 0; i < M; i++) {
		GRBLinExpr* cnstr_meanAssign = new GRBLinExpr();
		GRBLinExpr* cnstr_varAssign = new GRBLinExpr();

		for (uint16_t j = 0; j < N; j++) {
			currentMean = (double)jobs[j].first;
			currentVar = (double)jobs[j].second;

			cnstr_meanAssign->addTerms(&currentMean, &z[i][j], 1);
			cnstr_varAssign->addTerms(&currentVar, &z[i][j], 1);
		}

		model->addConstr(*cnstr_meanAssign, GRB_EQUAL, v[i]);
		model->addConstr(*cnstr_varAssign, GRB_EQUAL, u[i]);
	}

	// Each machine has at least one job
	for (uint16_t i = 0; i < M; i++) {
		GRBLinExpr* cnstr_eachMachineHasJob = new GRBLinExpr();

		for (uint16_t j = 0; j < N; j++) {
			model->addConstr(t[i] - (1 - w[i][j])*bigM_1, GRB_LESS_EQUAL, z[i][j]);
			model->addConstr(t[i] + (1 - w[i][j])*bigM_1, GRB_GREATER_EQUAL, z[i][j]);
			model->addConstr(z[i][j], GRB_LESS_EQUAL, w[i][j] * bigM_1);

			cnstr_eachMachineHasJob->addTerms(&pOne, &z[i][j], 1);
		}
		model->addConstr(*cnstr_eachMachineHasJob, GRB_GREATER_EQUAL, t[i]);
	}

	// Each job is scheduled once
	// Constraint won't be added if the branch and price solver creates this model with single machine
	if (loadedInstance.numOfMachines > 1) {
		for (uint16_t j = 0; j < N; j++) {
			GRBLinExpr* cnstr_eachJobOnce = new GRBLinExpr();
			for (uint16_t i = 0; i < M; i++) {
				cnstr_eachJobOnce->addTerms(&pOne, &w[i][j], 1);
			}
			model->addConstr(*cnstr_eachJobOnce, GRB_EQUAL, 1.0);
		}
	}

	// Variance constraint transformation
	double helpCoef = 0.0;
	double deltaH_half;

	for (uint16_t i = 0; i < M; i++) {
		GRBLinExpr* cnstr_varValue = new GRBLinExpr();
		GRBLinExpr* cnstr_onlyOneInterval = new GRBLinExpr();

		for (uint16_t k = 0; k < L; k++) {
			cnstr_varValue->addTerms(&pOne, &s[i][k], 1);
			cnstr_onlyOneInterval->addTerms(&pOne, &y[i][k], 1);

			model->addConstr(s[i][k], GRB_GREATER_EQUAL, (1 / h[k]) - (1 / (h[k] * h[k]))*(t[i] - h[k]) - (1 - y[i][k])*bigM_2);
			model->addConstr(s[i][k], GRB_LESS_EQUAL, (1 / h[k]) - (1 / (h[k] * h[k]))*(t[i] - h[k]) + (1 - y[i][k])*bigM_2);
			model->addConstr(s[i][k], GRB_LESS_EQUAL, y[i][k] * bigM_2);

			// Enforce corresponding interval
			deltaH_half = ((h_bounds[k + 1] - h_bounds[k]) / 2);
			model->addConstr(t[i] - h[k], GRB_GREATER_EQUAL, -(1 - y[i][k]) * bigM_2 - deltaH_half);
			model->addConstr(t[i] - h[k], GRB_LESS_EQUAL, (1 - y[i][k]) * bigM_2 + deltaH_half);
		}
		model->addConstr(*cnstr_onlyOneInterval, GRB_EQUAL, 1.0);
		model->addConstr(*cnstr_varValue, GRB_EQUAL, u[i]);
	}

	// Symmetry Break 
	for (uint16_t i = 0; i < M - 1; i++) {
		model->addConstr(v[i], GRB_GREATER_EQUAL, v[i + 1]);
	}
}

void MILPsolver_primary::buildObjective() {
	// Criterion, make this custom function to work on model
	uint64_t intervalSize	= objFunctionApproxPoints.size();
	double*cdfAppr			= new double[intervalSize];
	double*x				= new double[intervalSize];

	for (uint32_t i = 0; i < intervalSize; i++) {
		x[i]		= objFunctionApproxPoints[i];
		cdfAppr[i]	= log10(cumulativeNormalDistribution(x[i]));
	}

	GRBVar* criterionVariable = new GRBVar[loadedInstance.numOfMachines];
	for (uint16_t i = 0; i < loadedInstance.numOfMachines; i++) {
		criterionVariable[i] = model->addVar(x[0], x[intervalSize - 1], 0.0, GRB_CONTINUOUS, "c" + to_string(i));
		model->addConstr(criterionVariable[i], GRB_EQUAL, loadedInstance.deadline*t[i] - v[i]);
		model->setPWLObj(criterionVariable[i], intervalSize, x, cdfAppr);
	}
	model->set(GRB_IntAttr_ModelSense, GRB_MAXIMIZE);

	// Set timeout
	model->set(GRB_DoubleParam_TimeLimit, this->timeout);

	delete x, cdfAppr;
}

void MILPsolver_primary::loadSolverConfig(config_MILPprimary c) {
	this->loadedInstance	= c.loadedInstance;
	this->conf				= c;
	this->conf.vMin			= calculate_vMin();
	this->conf.vMax			= calculate_vMax();
	this->timeout			= c.timeout;

	fracApproxIntervalBounds	= approxSteps_fraction();
	objFunctionApproxPoints		= approxSteps_objFunction();	
}

Schedule MILPsolver_primary::solve() {
	model->set(GRB_IntParam_OutputFlag, 0);
	model->update();
	model->optimize();

	// Check if solution found
	//if (model->get(GRB_IntAttr_SolCount) == 0) continue;

	this->expandedNodes = (uint64_t)model->get(GRB_DoubleAttr_NodeCount);
	this->nodesToBest	= 0;
	this->runtime		= model->get(GRB_DoubleAttr_Runtime);
	this->objGap		= model->get(GRB_DoubleAttr_MIPGap);

	// Return schedule
	Schedule result(loadedInstance.numOfMachines);
	for (uint16_t i = 0; i < loadedInstance.numOfMachines; i++) {
		for (uint16_t j = 0; j < loadedInstance.numOfJobs; j++) {
			if (w[i][j].get(GRB_DoubleAttr_X) > 0.6) {
				result[i].push_back(loadedInstance.jobs[j]);
			}
		}
	}

	return result;
}

vector<double> MILPsolver_primary::approxSteps_fraction() {
	// Here only boundary interval points are calcualted, the interpolation is done in model construction
	double startPoint, endPoint, stdSum;
	vector<uint16_t> variances;

	for (Job j : loadedInstance.jobs) { variances.push_back(j.second); }
	sort(variances.begin(), variances.end());

	stdSum = 0;
	for (uint16_t i = 0; i < conf.vMin; i++) { stdSum += sqrt(variances[i]); }
	endPoint = 1 / stdSum;

	stdSum = 0;
	for (uint16_t i = variances.size() - 1; i >= variances.size() - conf.vMax; i--) { stdSum += sqrt(variances[i]); }
	startPoint = 1 / stdSum;


	// Approximate interval 
	vector<double> approximationPoints;
	approximationPoints.push_back(startPoint);
	while (startPoint < endPoint) {		
		startPoint = startPoint / (1 - conf.fracFunctionStep * startPoint);
		approximationPoints.push_back(startPoint); // So the first one reaching beyond endPoint is added as well just in case
	}

	return approximationPoints;
}

vector<double> MILPsolver_primary::approxSteps_objFunction() {
	double startPoint, endPoint, stdSum;
	vector<uint16_t> sortedMean, sortedVars;
	vector<double> approximationPoints;
	int16_t meanSum, criterion;

	for (Job j : loadedInstance.jobs) {
		sortedMean.push_back(j.first);
		sortedVars.push_back(j.second);
	}
	sort(sortedMean.begin(), sortedMean.end());
	sort(sortedVars.begin(), sortedVars.end());

	// First min point
	stdSum = 0, meanSum = 0;
	for (uint16_t i = sortedMean.size() - 1; i >= sortedMean.size() - conf.vMax; i--) { meanSum += sortedMean[i]; }
	criterion = loadedInstance.deadline - meanSum;
	// Take either vMax lowest or highest variances/stds
	if (criterion > 0) { 
		for (uint16_t i = sortedMean.size() - 1; i >= sortedMean.size() - conf.vMax; i--) { stdSum += sqrt(sortedVars[i]); } 
	} else { 
		for (uint16_t i = 0; i < conf.vMax; i++) { stdSum += sqrt(sortedVars[i]); }
	}
	startPoint = (criterion / stdSum);

	// We can however calculate lower interval bound from the LB of the problem, by taking probit in the point of customer service level
	double startPoint_fromLB = probit(customerServiceLevel(calculateLB(loadedInstance), loadedInstance.deadline));
	startPoint = max(startPoint, startPoint_fromLB);

	// Now max point
	stdSum = 0, meanSum = 0;
	for (uint16_t i = 0; i < conf.vMin; i++) { meanSum += sortedMean[i]; }
	criterion = loadedInstance.deadline - meanSum;
	// Take either vMin lowest or highest variances/stds
	if (criterion > 0) {
		for (uint16_t i = 0; i < conf.vMin; i++) { stdSum += sqrt(sortedVars[i]); }		
	} else {
		for (uint16_t i = sortedMean.size(); i > sortedMean.size() - conf.vMin; i--) { stdSum += sqrt(sortedVars[i]); }
	}
	endPoint = (criterion / stdSum);

	// Create interval
	approximationPoints.push_back(startPoint);
	while (startPoint < endPoint) {
		startPoint += conf.objFunctionStep;
		approximationPoints.push_back(startPoint); // So that last one excees the end point just in case
	}

	return approximationPoints;
}

// In these functions remove previous settings as well

void MILPsolver_primary::preloadSolution(Schedule s) {
	if (s.empty()) {
		return;
	}

	for (uint16_t i = 0; i < loadedInstance.numOfMachines; i++) {
		for (Job jb : s[i]) {
				w[i][jb.id].set(GRB_DoubleAttr_Start, 1.0);
		}
	}
}

void MILPsolver_primary::enforceSolution(Schedule s) {
	if (s.empty()) {
		return;
	}

	for (uint16_t i = 0; i < loadedInstance.numOfMachines; i++) {
		for (Job j : s[i]) {
			w[i][j.id].set(GRB_DoubleAttr_LB, 1.0);
		}
	}
}
