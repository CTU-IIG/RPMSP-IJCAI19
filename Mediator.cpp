#include "Mediator.h"

Mediator::Mediator() {
	this->configLoaded = false;
}

Mediator::~Mediator() {
	delete(slvr);
	outputFile.close();
}

void Mediator::createSolver(string solverType, config_MILPprimary c) {
	// Config passed as arugment to constructors

	if (solverType == "MILP") {
		slvr = new MILPsolver_primary(c);
		this->loadedConfig = c;
		this->configLoaded = true;
	} else {
		cerr << "Wrong combination of solver type and config type!" << endl;
	}

	outputFile.open("output.csv", ofstream::out | ofstream::app);
	if (!outputFile) { cerr << "Unable to open output file!"; }
}

void Mediator::createSolver(string solverType, config_BranchAndPrice c) {
	if (solverType == "BP") {
		slvr = new MILP_BranchAndPrice(c);
		this->loadedConfig = c;
		this->configLoaded = true;
	}
	else {
		cerr << "Wrong combination of solver type and config type!" << endl;
	}

	outputFile.open("output.csv", ofstream::out | ofstream::app);
	if (!outputFile) { cerr << "Unable to open output file!"; }
}

Schedule Mediator::solve() {
	if (configLoaded) {
		result = slvr->solve();

		optVal			= customerServiceLevel(result, loadedConfig.loadedInstance.deadline);
		expandedNodes	= slvr->getExpandedNodes();
		nodesToBest		= slvr->getNodesToBest();
		runtime			= slvr->getRuntime();
		objGap			= slvr->getObjGap();
	} else {
		cerr << "You must load some solver configuration first!" << endl;
		return Schedule();
	}

	return result;
}

void Mediator::printSchedule(Schedule s) {
	for (uint16_t i = 0; i < s.size(); i++) {
		cout << "M" << i << ": ";
		for (Job j : s[i]) {
			cout << "J" << j.id << " ";
		} cout << endl;
	}
}

void Mediator::printResult() {
	cout << endl;
	printSchedule(result);
	cout << endl;
	cout << "Optimal value: " << optVal << endl;
	cout << "Total nodes expanded: " << expandedNodes << endl;
	cout << "Nodes to best: " << nodesToBest << endl;
	cout << "Total runtime: " << runtime << endl;
}

void Mediator::storeResult() {
	outputFile << loadedConfig.instancePath << ";";
	outputFile << optVal << ";";
	outputFile << expandedNodes << ";";
	outputFile << nodesToBest << ";";
	outputFile << runtime << ";";
	outputFile << objGap;
	outputFile << endl;
}
