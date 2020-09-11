#include "Mediator.h"

Mediator::Mediator(string solverType) {
	if (solverType == "MILP") {
		slvr = new MILPsolver_primary();
	}
	else if (solverType == "MILP2") {
		slvr = new MILPsolver_secondary();
	}

	outputFile.open("output.csv", ofstream::out | ofstream::app);
	if (!outputFile) { cerr << "Unable to open output file!"; }
}

Mediator::Mediator(string instanceFileName, string solverType) {	
	if (solverType == "MILP") {
		slvr = new MILPsolver_primary();
	}
	else if (solverType == "MILP2") {
		slvr = new MILPsolver_secondary();
	}
	loadInstance(instanceFileName);

	outputFile.open("output.csv", ofstream::out | ofstream::app);
	if (!outputFile) { cerr << "Unable to open output file!"; }
}

Mediator::~Mediator() {
	delete(slvr);
	outputFile.close();
}

void Mediator::loadInstance(string fileName) {
	instance inst;
	inst.fileName = fileName;
	cout << "Reading instance: Reading file..." << endl;	

	ifstream inputFile(fileName);
	if (inputFile) {
		uint16_t mean, var;
		inputFile >> inst.numOfJobs >> inst.numOfMachines >> inst.deadline;
		
		for (uint16_t i = 0; i < inst.numOfJobs; i++) {
			inputFile >> mean >> var;
			if (mean == 0) mean = 1;
			if (var == 0) var = 1;

			inst.jobs.push_back(Job{ mean,var,i });
		}
		instanceLoaded = true;
		cout << "Reading instance: Done!" << endl;
	} else {
		cerr << "Unable to open input file!" << endl;
		cin.get();
		inst.jobs = JobList();
	}

	this->currentInstace = inst;
	slvr->setupSolver(currentInstace);
}

Schedule Mediator::solve() {
	if (instanceLoaded) {
		result = slvr->solve();

		optVal			= customerServiceLevel(result, currentInstace.deadline);
		expandedNodes	= slvr->getExpandedNodes();
		nodesToBest		= slvr->getNodesToBest();
		runtime			= slvr->getRuntime();
		fracSteps		= slvr->getFracSteps();
		objSteps		= slvr->getObjSteps();
	}
	else {
		cerr << "You must load some instance first!" << endl;
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
	outputFile << currentInstace.fileName << ";";
	outputFile << optVal << ";";
	outputFile << expandedNodes << ";";
	outputFile << nodesToBest << ";";
	outputFile << runtime << ";";
	outputFile << fracSteps << ";";
	outputFile << objSteps << ";";
	outputFile << slvr->getGap() << ";";
	outputFile << endl;
}
