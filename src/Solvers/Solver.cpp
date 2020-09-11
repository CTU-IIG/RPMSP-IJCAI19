#include "Solver.h"
#include "statiscticsLib.h"

Solver::Solver() {

}

Solver::~Solver() {

}

void Solver::preloadSolution(Schedule s) {

}

void Solver::enforceSolution(Schedule s) {

}

void Solver::loadConfig(MILPconfig conf) {

}

Schedule Solver::calculateLB(instance i) {
	uint16_t numOfMachines = i.numOfMachines;
	JobList	 jobs = i.jobs;
	//sort(jobs.begin(), jobs.end(), [](Job one, Job two) -> bool { return one.second*one.first > two.second*two.first; });
	// Push sorting to setupInstance and activate depedning on config file

	vector<uint16_t> cumulativeMeans(numOfMachines, 0);
	vector<uint16_t> cumulativeVars(numOfMachines, 0);
	Schedule tempSchedule = Schedule(numOfMachines, JobList());

	uint16_t sum_mi = 0;
	uint16_t rhoIndex = 0;
	double averageLoad = 0.0;
	double rhoMin = 0.0;

	vector<double> rhoVal = vector<double>(numOfMachines, 0.0);

	for (uint16_t k = 0; k < numOfMachines; k++) {
		if (jobs.empty())
			break;

		tempSchedule[k].push_back(jobs.front());
		cumulativeMeans[k] += jobs.front().first;
		cumulativeVars[k] += jobs.front().second;

		sum_mi += jobs.front().first;
		jobs.erase(jobs.begin());
	}

	while (!jobs.empty()) {
		averageLoad = (double)sum_mi / (double)numOfMachines;
		rhoMin = DBL_MAX;

		for (uint16_t i = 0; i < numOfMachines; i++) {
			rhoVal[i] = (cumulativeMeans[i] - averageLoad) / sqrt(cumulativeVars[i]);

			if (rhoVal[i] < rhoMin) {
				rhoMin = rhoVal[i];
				rhoIndex = i;
			}
		}

		sum_mi += jobs.front().first;
		cumulativeMeans[rhoIndex] += jobs.front().first;
		cumulativeVars[rhoIndex] += jobs.front().second;
		tempSchedule[rhoIndex].push_back(jobs.front());

		jobs.erase(jobs.begin());
	}
	
	return tempSchedule;
}

