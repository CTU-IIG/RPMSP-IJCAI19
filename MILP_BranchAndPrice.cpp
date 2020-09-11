#include "MILP_BranchAndPrice.h"

// If constraint cannot be found, 
// Compare two lower bounds

/*
	Nodes expanded, runtimes, timeouted instances
*/

// Check introduction
// 

MILP_BranchAndPrice::MILP_BranchAndPrice(config_BranchAndPrice c) {
	this->usePricingHeuristic	= c.use_heuristic;
	this->useRegression			= c.use_regressionModel;

	this->loadedInstance	= c.loadedInstance;	
	this->timeout			= c.timeout;
	this->timed_out			= false;
	this->timer				= Timer();

	this->iterationsPerformed = 0;
	
	c.pricingProblemConfig.vMax = calculate_vMax();
	c.pricingProblemConfig.vMin = calculate_vMin();
	
	this->conf					= c;
	this->totalRuntime_master	= 0.0;
	this->totalRuntime_pricing	= 0.0;
	this->nodesToBest			= 0;
	this->expandedNodes			= 0;
	this->objGap				= 0.0;

	this->totalPatternsGenerated		= 0;
	this->heuristicPatternsGenerated	= 0;
	
	jobsScheduledInPatterns = vector<uint16_t>(loadedInstance.numOfJobs, 0);

	// Build first wrapper
	vector<member> initialMemberStatus(loadedInstance.numOfJobs, member{ 0, 0 });
	vector<uint16_t> componentAffiliations(loadedInstance.numOfJobs);
	for (uint16_t i = 0; i < loadedInstance.numOfJobs; i++) {
		initialMemberStatus[i].first = i;
		componentAffiliations[i]	 = i;
	}
	pendingIterations.push(iterationWrapper{vector<constraint>(), initialMemberStatus, componentAffiliations, 0});

	// Enforce constraints
	// If i has to be with j, join them
	// If cant be together, take only first found
	timeOutput.open("pricingTimeInfo.csv");

	string graphFileName = "graphFile_" + string(this->conf.instancePath.end() - 15, this->conf.instancePath.end());
	nodeCounter	= 0;
	graphFile.open(graphFileName);
	graphFile << "graph constraintGraph {" << endl;
}

MILP_BranchAndPrice::~MILP_BranchAndPrice() {

}

void MILP_BranchAndPrice::loadSolverConfig(config_BranchAndPrice c) {
	this->conf = c;
}

Schedule MILP_BranchAndPrice::solve() {
	if (this->conf.pricingProblem_solverVersion == 2 && this->conf.use_regressionModel == true) {
		cerr << "Cannot run regression and enumeratin pricing at the same time!" << endl;
		return Schedule();
	}
	
	timer.start();

	cPatterns iterationResult;
	Schedule masterSchedule, currentBestSchedule;
	double currentBestObjVal, iterationObjVal;

	Schedule LBschedule				= calculateLB(loadedInstance);		
	patterns iterationPatterns		= getPatternsFromSchedule(LBschedule);
	costs iterationServiceLevels	= getPricesFromSchedule(LBschedule);

	pendingIterations.top().nodePatterns		= iterationPatterns;
	pendingIterations.top().nodePatternCosts	= iterationServiceLevels;

	// Initial best solution
	currentBestSchedule = LBschedule;
	currentBestObjVal	= customerServiceLevel(currentBestSchedule, this->loadedInstance.deadline);

	while (!pendingIterations.empty()) {
		iterationPatterns		= pendingIterations.top().nodePatterns;
		iterationServiceLevels	= pendingIterations.top().nodePatternCosts;

		// Check patterns if they dont violate constraints and if they do, regenerate them
		checkPatternsConstraintViolation(&iterationPatterns, &iterationServiceLevels);

		iterationResult = runNextIteration(&iterationPatterns, &iterationServiceLevels, &iterationObjVal);
		if (timed_out) { break; }
		if (iterationObjVal < currentBestObjVal) { 
			pendingIterations.pop(); 
			continue; 
		}	// Cut-off, dont branch
			
		masterSchedule	= getScheduleFromMaster(iterationResult);
				
		if (masterSchedule.size() == 0) {
			cout << "Branching..." << endl;
			if(!branchOnConstraints(iterationPatterns, iterationServiceLevels)){
				continue;
			} 
			// This adds 2 new iterations to the stack, pops from inside, they inherit patterns and costs from parent
			cout << "Done!" << endl;
		}
		else { // Feasible solution found
			if (customerServiceLevel(masterSchedule, this->loadedInstance.deadline) > currentBestObjVal) {
				currentBestSchedule = masterSchedule;
				currentBestObjVal	= customerServiceLevel(currentBestSchedule, this->loadedInstance.deadline);
			}			
			pendingIterations.pop();
		}		
	}

	timer.stop();
	this->runtime = timer.getElapsed();

	this->expandedNodes = totalPatternsGenerated;
	this->nodesToBest	= heuristicPatternsGenerated;
	timeOutput << this->totalRuntime_master << ";" << this->totalRuntime_pricing << endl;
	graphFile << "}" << endl;
	timeOutput.close();
	graphFile.close();

	// Patterns diagnostic
	/*
	uint32_t total = 0;
	for (uint16_t i = 0; i < loadedInstance.numOfJobs; i++) {
		cout << i + 1 << ": " << jobsScheduledInPatterns[i] << endl;
		total += jobsScheduledInPatterns[i];
	} cout << endl;
	cout << totalRuntime_pricing << "   " << total << endl;
	cout << totalRuntime_pricing / (double)total << endl;
	*/

	return currentBestSchedule;
}

cPatterns MILP_BranchAndPrice::runNextIteration(patterns* availablePatterns, costs* serviceLevels, double* objVal) {
		
	cPatterns result;
	cPattern  newPattern, heuristicPattern;
	costs	prices(this->loadedInstance.numOfJobs, 0.0);
	double	gamma, sum, UB;

	// Regression
	vector<vector<double>>	reg_features;
	vector<double>			reg_targets;

	if (this->conf.pricingProblem_solverVersion == 1) {
		generatePricingModel(); // Builds pricing model specified by config, applies current constraints
		pricingModel->set(GRB_IntParam_OutputFlag, 0);
	}
	

	while (true) {
		result	= solveMasterProblem(*serviceLevels, *availablePatterns, &prices, &gamma, objVal);

		// Check feasibility of master
		if (result.first.empty()) {
			// Need to generate additional patterns due to model infeasibility
			if (generateAdditionalPatterns(availablePatterns, serviceLevels) == GRB_INFEASIBLE) {
				*objVal = 0.0; // This forces the cutoff;
				return cPatterns();
			}
			else { continue; }			
		}

		// Heuristic
		// Change back to new pattern later
		if (usePricingHeuristic) {
			heuristicPattern = solvePricingProblem_heuristic(prices, gamma);

			// Check negative reduced cost
			sum = -gamma;
			for (uint16_t i = 0; i < prices.size(); i++) {
				sum += prices[i] * heuristicPattern.first[i];
			}

			if (sum + heuristicPattern.second > this->conf.pricingProblem_criterionCoef) {
				// Heuristic solution stands
				heuristicPattern.second = log10(customerServiceLevel(newPattern.first, this->conf.loadedInstance.deadline, this->loadedInstance.jobs));
				availablePatterns->push_back(heuristicPattern.first);
				serviceLevels->push_back(heuristicPattern.second);				
				continue;			
			}
		} 

		if (useRegression) {
			if (iterationsPerformed >= this->conf.reg_iterationsNeeded) {
				vector<double> weights = reg_train(reg_features, reg_targets);
				UB = reg_getPrediction(weights, reg_features[reg_features.size() - 1]);

				pricingModel->set(GRB_DoubleParam_Cutoff, UB);
				cout << "UB: " << UB << endl;
			}
			reg_features.push_back(reg_getFeatures(prices));
		}

		// Update prices
		if (this->conf.pricingProblem_solverVersion == 1) {
			updatePricesOnPricingModel(prices);
			newPattern = solvePricingProblem();
		}
		else if (this->conf.pricingProblem_solverVersion == 2) {
			newPattern = solvePricingProblem_enumeration(prices);
		}
		
		if (useRegression) { 
			reg_targets.push_back(pricingModel->get(GRB_DoubleAttr_ObjVal));

			if (iterationsPerformed >= this->conf.reg_iterationsNeeded) {				
				cout << "Real: " << pricingModel->get(GRB_DoubleAttr_ObjVal) << endl;
			}
		}

		// Check constraint
		sum = -gamma;
		for (uint16_t i = 0; i < prices.size(); i++) { sum += prices[i] * newPattern.first[i]; }
		if (sum + newPattern.second < this->conf.pricingProblem_criterionCoef) {
			break;
		}

		availablePatterns->push_back(newPattern.first);
		serviceLevels->push_back(newPattern.second);

		// Check for timeout
		if (timer.getElapsed() > this->timeout) {
			timed_out = true;
			break;
		}
		iterationsPerformed++;
	}

	return result;
}

cPatterns MILP_BranchAndPrice::solveMasterProblem(vector<double> c, vector<pattern> patterns, vector<double>* jobPrices, double* gammaValue, double* objVal) {
	// Model is being recreated, could be reused instead if implemented it dual version

	if (patterns.size() < this->conf.loadedInstance.numOfMachines) {
		return cPatterns(vector<pattern>(), costs());
	}

	GRBEnv	 env;
	GRBModel model(env);
	uint16_t numOfJobs = patterns[0].size();

	GRBVar* y		= new GRBVar[patterns.size()];
	GRBConstr* phi	= new GRBConstr[numOfJobs];
	GRBConstr gamma;
	
	double one = 1.0;
	GRBLinExpr* gammaConstr = new GRBLinExpr();
	for (uint16_t k = 0; k < patterns.size(); k++) {
		y[k] = model.addVar(0.0, GRB_INFINITY, c[k], GRB_CONTINUOUS);
		gammaConstr->addTerms(&one, &y[k], 1);
	}
	gamma = model.addConstr(*gammaConstr, GRB_LESS_EQUAL, loadedInstance.numOfMachines);

	for (uint16_t j = 0; j < numOfJobs; j++) {
		GRBLinExpr* firstConstr = new GRBLinExpr();

		for (uint16_t k = 0; k < patterns.size(); k++) {
			firstConstr->addTerms(&patterns[k][j], &y[k], 1);
		}
		phi[j] = model.addConstr(*firstConstr, GRB_GREATER_EQUAL, 1.0);
	}
	
	model.set(GRB_IntParam_OutputFlag, 0);
	model.set(GRB_IntAttr_ModelSense, GRB_MAXIMIZE);
	model.update();
	model.optimize();

	// If model infeasible, return empty
	if (model.get(GRB_IntAttr_Status) == GRB_INFEASIBLE || model.get(GRB_IntAttr_Status) == GRB_INF_OR_UNBD) {
		return cPatterns(vector<pattern>(), costs());
	}

	for (uint16_t i = 0; i < numOfJobs; i++) {
		(*jobPrices)[i] = (-phi[i].get(GRB_DoubleAttr_Pi)); 
	}
	*gammaValue = gamma.get(GRB_DoubleAttr_Pi);
	*objVal		= pow(10,model.get(GRB_DoubleAttr_ObjVal));

	vector<pattern> allPatterns;
	vector<double> activations;
	for (uint16_t i = 0; i < patterns.size(); i++) {
		allPatterns.push_back(patterns[i]);
		activations.push_back(y[i].get(GRB_DoubleAttr_X));
	}

	// Update diagnostics and return
	this->totalRuntime_master += model.get(GRB_DoubleAttr_Runtime);
	//this->expandedNodes		  += round(model.get(GRB_DoubleAttr_NodeCount));
	return pair<vector<pattern>, vector<double>>(allPatterns, activations);
}

cPattern MILP_BranchAndPrice::solvePricingProblem() {
	pricingModel->optimize();

	if (pricingModel->get(GRB_IntAttr_Status) == GRB_INFEASIBLE) {
		cout << "Pricing infeasible!" << endl;
		return cPattern();
	}

	// Extract pattern
	cPattern newPattern;
	pattern	 p(this->loadedInstance.numOfJobs, 0.0);
	double   patternPrice;
	JobList  j;

	uint16_t jobsScheduled = 0;
	for (uint16_t i = 0; i < this->loadedInstance.numOfJobs; i++) {
		if (x[i].get(GRB_DoubleAttr_X) >= 0.9) {
			p[i] = 1.0;
			jobsScheduled++;
			j.push_back(this->loadedInstance.jobs[i]);
		}
	}
	patternPrice = log10(customerServiceLevel(j, this->loadedInstance.deadline));
	newPattern	 = cPattern(p, patternPrice);

	jobsScheduledInPatterns[jobsScheduled]++;

	this->totalRuntime_pricing += pricingModel->get(GRB_DoubleAttr_Runtime);
	//this->expandedNodes		   += round(pricingModel->get(GRB_DoubleAttr_NodeCount));
	return newPattern;
}

cPattern MILP_BranchAndPrice::solvePricingProblem_heuristic(costs prices, double gamma) {
	double	patternCost, currentCost, phiAccumulated;
	pattern newPattern(prices.size(), 0);
	JobList tempJobList;
	vector<constraint> currentConstraints;

	struct virtualJob {
		uint8_t id;
		double weight;
	};
	vector<virtualJob> relativeCosts;

	for (uint8_t i = 0; i < prices.size(); i++) {
		currentCost = prices[i] / this->conf.loadedInstance.jobs[i].first; // Cost from master divided by mean
		relativeCosts.push_back(virtualJob{ i, currentCost });
	}
	sort(relativeCosts.begin(), relativeCosts.end(), [](virtualJob v1, virtualJob v2) {return  v1.weight > v2.weight; });

	double bestCost = -DBL_MAX;
	pattern bestPattern;

	currentConstraints	= pendingIterations.top().applicableConstraints;
	patternCost			= 0;
	currentCost			= 0;
	phiAccumulated		= 0;
	vector<bool> forbidden(this->conf.loadedInstance.numOfJobs, false);
	for (virtualJob vj : relativeCosts) {		
		if (forbidden[vj.id]) continue;
				
		newPattern[vj.id] = 1.0;

		// Check constraints
		for (constraint c : currentConstraints) {
			if (vj.id == c.first_member_ID) {
				if (c.relation = CNSTR_FORCE_PAIR) {
					newPattern[c.second_member_ID] = 1.0;
					phiAccumulated += prices[c.second_member_ID];
				}
				forbidden[c.second_member_ID] = true; // Applies to forbid pair constraint, also to the force pair so it is not scheduled twice
			}

			if (vj.id == c.second_member_ID) {
				if (c.relation = CNSTR_FORCE_PAIR) {
					newPattern[c.first_member_ID] = 1.0;
					phiAccumulated += prices[c.first_member_ID];
				}
				forbidden[c.first_member_ID] = true;
			}
		}

		phiAccumulated   += prices[vj.id];
 		tempJobList.push_back(this->conf.loadedInstance.jobs[vj.id]);

		currentCost = log10(customerServiceLevel(tempJobList, this->conf.loadedInstance.deadline));
		currentCost += phiAccumulated;

		// Take max
		if (currentCost > bestCost) {
			bestCost	= currentCost;
			bestPattern = newPattern;
		}
	}

	return cPattern(bestPattern, log10(customerServiceLevel(bestPattern, this->conf.loadedInstance.deadline, this->conf.loadedInstance.jobs)));
}


vector<pattern> MILP_BranchAndPrice::getPatternsFromSchedule(Schedule s) {
	vector<pattern> patterns(s.size(), pattern(loadedInstance.numOfJobs, 0));
	
	for (uint16_t i = 0; i < s.size(); i++) {
		for (Job j : s[i]) {
			patterns[i][j.id] = 1.0;
		}
	}

	return patterns;
}

vector<double> MILP_BranchAndPrice::getPricesFromSchedule(Schedule s) {
	vector<double> prices(s.size());
	for (uint16_t i = 0; i < s.size(); i++) {
		prices[i] = log10(customerServiceLevel(s[i], loadedInstance.deadline));
	}

	return prices;
}

Schedule MILP_BranchAndPrice::getScheduleFromMaster(cPatterns masterSolution) {
	Schedule result(this->loadedInstance.numOfMachines);
	uint16_t machineCounter = 0;
	
	for (uint16_t i = 0; i < masterSolution.second.size(); i++) {
		if (masterSolution.second[i] == 1.0) {
			// Activated pattern, extract
			for (uint16_t j = 0; j < masterSolution.first[i].size(); j++) {
				if (masterSolution.first[i][j] == 1.0) {
					result[machineCounter].push_back(this->loadedInstance.jobs[j]);
				}
			}
			machineCounter++;
		}
	}
	
	if (machineCounter != this->loadedInstance.numOfMachines) {
		// Incomplete solution, need to apply constraints, handled upstream
		return Schedule();
	}

	return result;
}

void MILP_BranchAndPrice::interpretConstraint(constraint c, GRBVar* x, GRBModel* model) {
	if (c.relation == CNSTR_FORBID_PAIR) {
		GRBLinExpr expression;
		double one = 1.0;

		expression.addTerms(&one, &x[c.first_member_ID], 1);
		expression.addTerms(&one, &x[c.second_member_ID], 1);
		model->addConstr(expression, GRB_LESS_EQUAL, 1.0);
	}
	else if (c.relation == CNSTR_FORCE_PAIR) {
		model->addConstr(x[c.first_member_ID], GRB_EQUAL, x[c.second_member_ID]);
	}
}

void MILP_BranchAndPrice::generatePricingModel() {
	// Make this into switch
	if (this->conf.pricingProblem_solverVersion == 1) {
		milp		 = new MILPsolver_primary(this->conf.pricingProblemConfig, costs(this->loadedInstance.numOfJobs, 0));
		x			 = milp->getAssignmentVariable();
		pricingModel = milp->getModel();
	}
	else if (this->conf.pricingProblem_solverVersion == 2) {
		// Generate the second pricing model

	}

	vector<constraint> currentConstraints = pendingIterations.top().applicableConstraints;
	for (constraint c : currentConstraints) {
		interpretConstraint(c, x, pricingModel);
	}
}

void MILP_BranchAndPrice::updatePricesOnPricingModel(costs prices) {
	for (uint16_t i = 0; i < this->loadedInstance.numOfJobs; i++) {
		x[i].set(GRB_DoubleAttr_Obj, prices[i]);
	}
}

void MILP_BranchAndPrice::checkPatternsConstraintViolation(patterns * p, costs * c) {
	patterns preservedPatterns;
	costs	 preservedCosts;	
	bool	 discard;

	pattern  pat;
	for (uint64_t i = 0; i < p->size(); i++) {
		pat		= (*p)[i];
		discard = false;
		// All constraints must hold for single pattern
		for (constraint cnstr : pendingIterations.top().applicableConstraints) {
			if (cnstr.relation == CNSTR_FORBID_PAIR) {
				if (pat[cnstr.first_member_ID] + pat[cnstr.second_member_ID] > 1.0 + 2*DBL_MIN) { // Tolerance
					discard = true;
					break;
				}
			}
			else if (cnstr.relation == CNSTR_FORCE_PAIR) {
				if (pat[cnstr.first_member_ID] != pat[cnstr.second_member_ID]) {
					discard = true;
					break;
				}
			}
		}

		if (!discard) {
			preservedPatterns.push_back(pat);
			preservedCosts.push_back((*c)[i]);
		}
	}

	// Check if enough patterns remained
	if (preservedPatterns.size() < this->loadedInstance.numOfMachines) {		
		cPatterns generatedPatterns = regeneratePatterns();

		preservedPatterns.insert(preservedPatterns.end(), generatedPatterns.first.begin(), generatedPatterns.first.end());
		preservedCosts.insert(preservedCosts.end(), generatedPatterns.second.begin(), generatedPatterns.second.end());
	}

	// Update
	(*p) = preservedPatterns;
	(*c) = preservedCosts;
}

int MILP_BranchAndPrice::generateAdditionalPatterns(patterns * p, costs * c){
	GRBEnv	 env;
	GRBModel m(env);
	double	 one = 1.0;

	GRBVar** X = new GRBVar*[this->loadedInstance.numOfMachines];
	for (uint16_t i = 0; i < this->loadedInstance.numOfMachines; i++) {
		GRBLinExpr le;

		X[i] = new GRBVar[this->loadedInstance.numOfJobs];
		for (uint16_t j = 0; j < this->loadedInstance.numOfJobs; j++) {
			X[i][j] = m.addVar(0.0, 1.0, 1.0, GRB_BINARY);
			le.addTerms(&one, &X[i][j], 1);
		}

		m.addConstr(le, GRB_LESS_EQUAL, this->conf.pricingProblemConfig.vMax);

		// Apply branching constraints
		for (constraint c : this->pendingIterations.top().applicableConstraints) {
			if (c.relation == CNSTR_FORBID_PAIR) {
				m.addConstr(X[i][c.first_member_ID] + X[i][c.second_member_ID], GRB_LESS_EQUAL, 1);
				//cout << "Adding " << c.first_member_ID << " cannot be with " << c.second_member_ID << endl;
			}
			else if (c.relation == CNSTR_FORCE_PAIR) {
				m.addConstr(X[i][c.first_member_ID], GRB_EQUAL, X[i][c.second_member_ID]);
				//cout << "Adding " << c.first_member_ID << " must be with " << c.second_member_ID << endl;
			}

			if (c.first_member_ID == 5 && c.second_member_ID == 0) {
				//cout << "Here!" << endl;
			}
		}
	}

	for (uint16_t j = 0; j < this->loadedInstance.numOfJobs; j++) {
		GRBLinExpr le;
		for (uint16_t i = 0; i < this->loadedInstance.numOfMachines; i++) {
			le.addTerms(&one, &X[i][j], 1);
		}

		m.addConstr(le, GRB_GREATER_EQUAL, 1);
	}

	m.set(GRB_IntParam_OutputFlag, 0);
	m.update();
	m.optimize();

	if (m.get(GRB_IntAttr_Status) == GRB_INFEASIBLE) {
		return GRB_INFEASIBLE;
	}

	JobList jl;
	pattern pat(this->loadedInstance.numOfJobs);
	
	for (uint16_t i = 0; i < this->loadedInstance.numOfMachines; i++) {
		jl = JobList();

		for (uint16_t j = 0; j < this->loadedInstance.numOfJobs; j++) {
			pat[j] = X[i][j].get(GRB_DoubleAttr_X);
			if (pat[j] > 0.9) {
				jl.push_back(this->loadedInstance.jobs[j]);
			}
		}

		// Add
		p->push_back(pat);
		c->push_back(log10(customerServiceLevel(jl, this->conf.loadedInstance.deadline)));
	}
	return 0;
}

bool MILP_BranchAndPrice::branchOnConstraints(patterns parentPatterns, costs patternCosts) {
	vector<member> degrees	    = pendingIterations.top().memberDegree;	
	vector<uint16_t> components = pendingIterations.top().affiliationMap;
	vector<constraint> constraints;

	// Sort the degrees with respect to member degree (second), lamba ftw
	sort(degrees.begin(), degrees.end(), [](member one, member two) -> bool {return one.second < two.second; });

	// Select appropriate variables
	uint16_t firstPosition	= 0;
	uint16_t secondPosition = 1;
	uint16_t firstMember, secondMember;
	bool constraintExists;

	// Check if constraint exists or if there already is some constraint "component" between those two
	do {
		cout << firstPosition << "   " << secondPosition << endl;
		if (firstPosition == 13 && secondPosition == 14) {
			cout << "Here!" << endl;
		}

		constraintExists	= false;
		firstMember			= degrees[firstPosition].first;
		secondMember		= degrees[secondPosition].first;			

		if (components[firstMember] - components[secondMember] == 0) {
			constraintExists = true;
		} else {
			for (constraint c : pendingIterations.top().applicableConstraints) {
				if (c.first_member_ID == firstMember && c.second_member_ID == secondMember) {
					constraintExists = true;
					break;
				}
				else if (c.first_member_ID == secondMember && c.second_member_ID == firstMember) {
					constraintExists = true;
					break;
				}
			}
		}
		
		if (constraintExists) {
			secondPosition++;

			if (secondPosition >= degrees.size()) {
				firstPosition++;
				secondPosition = firstPosition + 1;
			}

			if (firstPosition >= degrees.size()) {
				// Cannot generate additional constraints, add arbitrary constraint which will cause infeasibility
				return false;
			}
		}

		// So what happens if we cannot constraint on it more?
	} while (constraintExists);

	degrees[firstPosition].second++;
	degrees[secondPosition].second++;
	
	// Componentize
	uint16_t targetComponent = max(components[firstMember], components[secondMember]);
	uint16_t newComponentsNo = min(components[firstMember], components[secondMember]);
	for (uint16_t i = 0; i < components.size(); i++) {
		if (components[i] == targetComponent) { components[i] = newComponentsNo; }
	}

	// Add forbid pair constraint
	constraints				= pendingIterations.top().applicableConstraints;
	constraint c1			= {firstMember, secondMember, CNSTR_FORBID_PAIR};
	constraints.push_back(c1);
	iterationWrapper iw1	= { constraints, degrees, components, nodeCounter + 1, parentPatterns, patternCosts };

	// Add force pair constraint
	constraints				= pendingIterations.top().applicableConstraints;
	constraint c2			= { firstMember, secondMember, CNSTR_FORCE_PAIR };
	constraints.push_back(c2);
	iterationWrapper iw2	= { constraints, degrees, components, nodeCounter + 2, parentPatterns, patternCosts };
	
	// Record to graphViz	
	graphFile << to_string(nodeCounter + 1) << " [label = \"X" << to_string(c1.first_member_ID) << " + X" << to_string(c1.second_member_ID) << "<= 1\"];" << endl;
	graphFile << to_string(pendingIterations.top().nodeID) << " -- " << to_string(nodeCounter + 1) << ";" << endl;

	graphFile << to_string(nodeCounter + 2) << " [label = \"X" << to_string(c2.first_member_ID) << " = X" << to_string(c2.second_member_ID) << "\"];" << endl;
	graphFile << to_string(pendingIterations.top().nodeID) << " -- " << to_string(nodeCounter + 2) << ";" << endl;

	nodeCounter += 2;

	pendingIterations.pop();
	pendingIterations.push(iw1);
	pendingIterations.push(iw2);
	return true;
}

cPatterns MILP_BranchAndPrice::regeneratePatterns() {

	return cPatterns();
}


vector<double> MILP_BranchAndPrice::approxSteps_objFunction() {
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
	for (uint16_t i = sortedMean.size() - 1; i >= sortedMean.size() - conf.pricingProblemConfig.vMax; i--) { meanSum += sortedMean[i]; }
	criterion = loadedInstance.deadline - meanSum;
	// Take either vMax lowest or highest variances/stds
	if (criterion > 0) {
		for (uint16_t i = sortedMean.size() - 1; i >= sortedMean.size() - conf.pricingProblemConfig.vMax; i--) { stdSum += sqrt(sortedVars[i]); }
	}
	else {
		for (uint16_t i = 0; i < conf.pricingProblemConfig.vMax; i++) { stdSum += sqrt(sortedVars[i]); }
	}
	startPoint = (criterion / stdSum);

	// We can however calculate lower interval bound from the LB of the problem, by taking probit in the point of customer service level
	double startPoint_fromLB = probit(customerServiceLevel(calculateLB(loadedInstance), loadedInstance.deadline));
	startPoint = max(startPoint, startPoint_fromLB);

	// Now max point
	stdSum = 0, meanSum = 0;
	for (uint16_t i = 0; i < conf.pricingProblemConfig.vMin; i++) { meanSum += sortedMean[i]; }
	criterion = loadedInstance.deadline - meanSum;
	// Take either vMin lowest or highest variances/stds
	if (criterion > 0) {
		for (uint16_t i = 0; i < conf.pricingProblemConfig.vMin; i++) { stdSum += sqrt(sortedVars[i]); }
	}
	else {
		for (uint16_t i = sortedMean.size(); i > sortedMean.size() - conf.pricingProblemConfig.vMin; i--) { stdSum += sqrt(sortedVars[i]); }
	}
	endPoint = (criterion / stdSum);

	// Create interval
	approximationPoints.push_back(startPoint);
	while (startPoint < endPoint) {
		startPoint += conf.pricingProblemConfig.objFunctionStep;
		approximationPoints.push_back(startPoint); // So that last one excees the end point just in case
	}

	return approximationPoints;
}

vector<uint32_t> MILP_BranchAndPrice::getFeasibleVariances() {
	vector<uint32_t> availableVariances;
	vector<uint32_t> feasibleVariances;
	uint32_t minVar, maxVar;

	for (Job j : this->loadedInstance.jobs) {
		availableVariances.push_back(j.second);
	}
	sort(availableVariances.begin(), availableVariances.end());

	minVar = 0;
	for (uint16_t i = 0; i < conf.pricingProblemConfig.vMin; i++) {
		minVar += availableVariances[i];
	}

	maxVar = 0;
	for (uint16_t i = 0; i < conf.pricingProblemConfig.vMax; i++) {
		maxVar += availableVariances[availableVariances.size() - 1 - i];
	}

	for (uint16_t i = minVar; i <= maxVar; i++) {
		feasibleVariances.push_back(i);
	}

	return feasibleVariances;
}

vector<double> MILP_BranchAndPrice::reg_getFeatures(costs c) {
	uint16_t numOfJobs	= c.size();
	uint16_t h			= round(c.size() / this->loadedInstance.numOfMachines);
	vector<double> features(this->loadedInstance.numOfMachines, 0);

	double feature;
	for (uint16_t i = 0; i < this->loadedInstance.numOfMachines; i++) {
		feature = 0.0;
		for (uint16_t j = 0; j < h; j++) {
			if (i*h + j >= c.size()) { break; }
			
			feature += c[i*h + j] / this->loadedInstance.jobs[i*h + j].first;
		}
		features[i] = feature;
	}

	return features;
}

vector<double> MILP_BranchAndPrice::reg_train(vector<vector<double>> features, vector<double> targets) {
	uint16_t num_features = features[0].size();

	vector<double> weights;
	vector<double> lambda;

	// Build lamba vector
	for (uint64_t i = 0; i < features.size(); i++) {
		lambda.push_back(pow(10, ((double)i)/(features.size() - 1))); // Tweak this later
	}

	GRBEnv env;
	GRBModel regressionModel(env);

	GRBVar* r_minus = new GRBVar[targets.size()];
	GRBVar* r_plus	= new GRBVar[targets.size()];
	GRBVar* b_coef	= new GRBVar[targets.size()];
	GRBVar* w		= new GRBVar[features[0].size()];

	for (uint16_t i = 0; i < num_features; i++) {
		w[i] = regressionModel.addVar(-GRB_INFINITY, GRB_INFINITY, 0.0, GRB_CONTINUOUS);
	}

	for (uint16_t i = 0; i < targets.size(); i++) {
		r_minus[i]	= regressionModel.addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS);
		r_plus[i]	= regressionModel.addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS);
		b_coef[i]	= regressionModel.addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS);

		GRBLinExpr* constraint = new GRBLinExpr();
		for (uint16_t j = 0; j < num_features; j++) {
			constraint->addTerms(&features[i][j], &w[j], 1);
		}
		regressionModel.addConstr(*constraint, GRB_EQUAL, targets[i] + r_plus[i] - r_minus[i]);
		regressionModel.addGenConstrMax(b_coef[i], &r_plus[i], 1, 0.0);
	}

	GRBLinExpr* objective = new GRBLinExpr();
	GRBGenConstr* maxConstr;
	double alpha_1, alpha_2;
	for (uint16_t i = 0; i < features.size(); i++) {
		alpha_1 = lambda[i]*this->conf.reg_cMinus;
		alpha_2 = lambda[i]*this->conf.reg_cPlus;

		objective->addTerms(&alpha_1, &r_minus[i], 1);
		objective->addTerms(&alpha_2, &b_coef[i], 1);
	}
	regressionModel.setObjective(*objective, GRB_MINIMIZE);
	regressionModel.set(GRB_IntParam_OutputFlag, 0);
	regressionModel.update();
	regressionModel.optimize();

	for (uint16_t i = 0; i < num_features; i++) {
		weights.push_back(w[i].get(GRB_DoubleAttr_X));
	}
	return weights;
}

double MILP_BranchAndPrice::reg_getPrediction(vector<double> weights, vector<double> features) {
	double prediction = 0.0;
	for (uint16_t i = 0; i < weights.size(); i++) {
		prediction += weights[i] * features[i];
	}
	return prediction;
}

cPattern MILP_BranchAndPrice::solvePricingProblem_enumeration(costs prices) {
	uint16_t v_bar = 0;
	uint16_t v_min = UINT16_MAX;
	
	vector<double> sortedVar;
	for (Job j : this->loadedInstance.jobs) {
		v_bar += j.second; 
		sortedVar.push_back(j.second); 
		if (j.second < v_min) { v_min = j.second; }
	}
	this->expandedNodes = v_bar;

	v_bar = 0;
	sort(sortedVar.begin(), sortedVar.end(), greater<uint16_t>());
	for (uint16_t j = 0; j < this->conf.pricingProblemConfig.vMax; j++){
		v_bar += sortedVar[j];
	}
	this->nodesToBest = v_bar;

	patterns solutionPool_patterns(v_bar - v_min + 1, pattern());
	costs solutionPool_reducedCosts(v_bar - v_min + 1, DBL_MIN);

	omp_set_num_threads(this->conf._allocate_threads);
	#pragma omp parallel for
	for (uint16_t v = v_min; v <= v_bar; v++) {
		double sigma, mean, coef;

		coef = sqrt(v);
		GRBEnv env;
		GRBModel enumModel(env);

		GRBVar* x = new GRBVar[this->loadedInstance.jobs.size()];
		GRBVar mu;
		GRBLinExpr cnstr_v, cnstr_mu;

		mu = enumModel.addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS);
		for (uint16_t i = 0; i < this->loadedInstance.jobs.size(); i++) {
			x[i] = enumModel.addVar(0.0, 1.0, prices[i], GRB_BINARY);

			sigma	= (double)this->loadedInstance.jobs[i].second;
			mean	= (double)this->loadedInstance.jobs[i].first;
			cnstr_v.addTerms(&sigma, &x[i], 1);
			cnstr_mu.addTerms(&mean, &x[i], 1);
		}

		enumModel.addConstr(cnstr_v, GRB_EQUAL, v);
		enumModel.addConstr(cnstr_mu, GRB_EQUAL, mu);

		// Add piecewise objective
		vector<double> objPoints	= approxSteps_objFunction();
		uint64_t intervalSize		= objPoints.size();
		double*cdfAppr				= new double[intervalSize];
		double*t					= new double[intervalSize];

		for (uint32_t i = 0; i < intervalSize; i++) {
			t[i]		= objPoints[i];
			cdfAppr[i]	= log10(cumulativeNormalDistribution(t[i]));
		}

		GRBVar criterionVariable;
		criterionVariable = enumModel.addVar(t[0], t[intervalSize - 1], 0.0, GRB_CONTINUOUS);
		enumModel.addConstr(criterionVariable, GRB_EQUAL, (loadedInstance.deadline - mu) / coef);
		enumModel.setPWLObj(criterionVariable, intervalSize, t, cdfAppr);
		
		enumModel.set(GRB_IntAttr_ModelSense, GRB_MAXIMIZE);
		// Set timeout
		enumModel.set(GRB_DoubleParam_TimeLimit, this->timeout);
		enumModel.set(GRB_IntParam_OutputFlag, 0);

		// Add constraints
		vector<constraint> currentConstraints = pendingIterations.top().applicableConstraints;
		for (constraint c : currentConstraints){
		    interpretConstraint(c, x, &enumModel);
		}

		enumModel.optimize();

		if (enumModel.get(GRB_IntAttr_Status) != GRB_INFEASIBLE) {
			pattern output(prices.size(), 0.0);
			for (uint16_t i = 0; i < prices.size(); i++) {
				output[i] = x[i].get(GRB_DoubleAttr_X);
			}

			solutionPool_patterns[v - v_min] = output;
			solutionPool_reducedCosts[v - v_min] = enumModel.get(GRB_DoubleAttr_ObjVal);
		}		

		delete cdfAppr, t, x;
	}

	cPattern result;
	result.second = DBL_MIN;
	for (uint16_t i = 0; i < solutionPool_reducedCosts.size(); i++) {
		if (result.second < solutionPool_reducedCosts[i]) {
			result.first = solutionPool_patterns[i];
			result.second = solutionPool_reducedCosts[i];
		}
	}

	result.second = customerServiceLevel(result.first, this->conf.loadedInstance.deadline, this->conf.loadedInstance.jobs);
	return result;
}
