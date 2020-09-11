#include "generic.h"
#include "statiscticsLib.h"
#include "Mediator.h"


int main(int argc, char* argv[]) {
	vector<string> solverTypes = {"MILP","MILP2", "BP"};
	string fileName, solverType;
	
	if (argc == 1) {
		fileName = "taskset2/instance_typeM4N14C075_3.txt";

		solverType = solverTypes[3];
	} else {
		solverType	= argv[1];
		fileName	= argv[2];
	}	
		
	config_MILPprimary cmp(fileName, 0.01, 0.5);
	cmp.timeout = 3600.0;
	config_BranchAndPrice c(fileName, cmp);
	c.timeout		= 3600.0;
	c.use_heuristic = false;
	
	c.use_regressionModel = false;
	c.reg_iterationsNeeded = 6;
	c.reg_epsilonTolerance = 0.02;
	c.reg_cPlus		= 10;
	c.reg_cMinus	= 2;
	c.reg_lambda_lb = 1;
	c.reg_lambda_ub = 2;
	c._allocate_threads = 4;
	c.reg_type = 4;
	
	c.pricingProblem_solverVersion = 2;
	
	//Config c(fileName);
	//c.timeout = 3600.0;
	
	Mediator med;
	med.createSolver(solverType, c);
	med.solve();
	med.printResult();
	med.storeResult();
		
	cin.get();
	return 0;
}
