#include "Solver.h"

Solver::Solver() {

}

Solver::~Solver() {

}

void Solver::loadSolverConfig(Config c) {

}

void Solver::preloadSolution(Schedule s) {

}

void Solver::enforceSolution(Schedule s) {

}

Schedule Solver::calculateLB(instance i) {
	uint16_t numOfMachines = i.numOfMachines;
	JobList	 jobs = i.jobs;
	sort(jobs.begin(), jobs.end(), [](Job one, Job two) -> bool { return one.second*one.first > two.second*two.first; });

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

uint16_t Solver::calculate_vMin() {
	// Implement later
	return 1;
}

uint16_t Solver::calculate_vMax() {
	Schedule LBschedule = Solver::calculateLB(loadedInstance);
	double LBval = customerServiceLevel(LBschedule, loadedInstance.deadline);

	vector<uint16_t> means, vars;
	for (Job j : loadedInstance.jobs) {
		means.push_back(j.first);
		vars.push_back(j.second);
	}
	sort(means.begin(), means.end());
	sort(vars.begin(), vars.end());

	Schedule tempSchedule = Schedule(loadedInstance.numOfMachines);
	uint16_t vMax = 0;
	double	 currentOptVal = 1.0;
	Job		 virtualJob = { 0,0,0, };

	while (currentOptVal >= LBval) {
		virtualJob.first = means[vMax];
		virtualJob.second = vars[vMax];
		tempSchedule[0].push_back(virtualJob);

		currentOptVal = customerServiceLevel(tempSchedule, loadedInstance.deadline);
		vMax++;
	}

	// Subtract one since last added already broke the inequality
	cout << "vMax set to:" << vMax - 1 << endl;
	return vMax - 1;
}

