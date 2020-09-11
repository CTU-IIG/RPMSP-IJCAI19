#include "Config.h"

Config::Config() {

}

Config::Config(string fileName) {
	this->instancePath = fileName;
	loadInstance(fileName);

	timeout = DBL_MAX;
}

Config::~Config() {

}

void Config::loadInstance(string fileName) {
	instance inst;
	inst.fileName = fileName;
	cout << "Reading instance: Reading file..." << endl;

	ifstream inputFile(fileName);
	if (inputFile) {
		uint16_t mean, var;
		inputFile >> inst.numOfJobs >> inst.numOfMachines >> inst.deadline;

		for (uint16_t i = 0; i < inst.numOfJobs; i++) {
			inputFile >> mean >> var;

			// Safeguard for jobs with mean or variance equal to zero
			if (mean == 0) { mean = 1; }
			if (var == 0) { var = 1; }

			inst.jobs.push_back(Job{ mean,var,i });
		}
		cout << "Reading instance: Done!" << endl;
	}
	else {
		cerr << "Unable to open input file!" << endl;
		cin.get();
		inst.jobs = JobList();
	}

	instancePath	= fileName;
	loadedInstance	= inst;
}
