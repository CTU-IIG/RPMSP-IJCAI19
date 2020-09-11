#include "generic.h"
#include <ctime>
#include <fstream>

#include "statiscticsLib.h"
#include "Mediator.h"


int main(int argc, char* argv[]) {
	string solverType;
	vector<string> solverTypes = {"BB","MILP","MILP2"}, files;

	if (argc == 1) {
		string fileName = "instances/instance_typeM3N12C100_9.txt";
		files.push_back(fileName);
		solverType = solverTypes[1];
	} else {
		solverType	= argv[1];
		for (uint16_t i = 2; i < argc; i++) {
			files.push_back(argv[i]);
		}
	}
	
	//Mediator* medExact = new Mediator("BB");
	
	vector<double> objStep = { 0.005 };
	vector<double> fracStep = { 0.5 };
	MILPconfig conf = { 0,0,0,0 };

	for (double objs : objStep) {
		for (double fracs : fracStep) {
			conf.objFunctionStep	= objs;
			conf.fracFuncStep		= fracs;

			Mediator* med = new Mediator(solverType);
			for (string f : files) {
				med->loadInstance(f);
				//medExact->loadInstance(f);
				med->preloadSolution(Solver::calculateLB(med->getLoadedInstance()));
				//med->preloadSolution(medExact->solve());
				//med->enforceSolution(medExact->solve());

				med->loadConfig(conf);
				med->solve();
				med->printResult();
				med->storeResult();
			}
			delete(med);
		}
	}


	cin.get();
	return 0;
}