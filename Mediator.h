#pragma once
#include "generic.h"
#include "Config.h"
#include "MILPsolver_primary.h"
#include "MILPsolver_secondary.h"
#include "MILP_BranchAndPrice.h"

class Mediator {
	public:
		Mediator();
		~Mediator();

		void createSolver(string solverType, Config c);
		void createSolver(string solverType, config_MILPprimary c);
		void createSolver(string solverType, config_BranchAndPrice c);

		void printSchedule(Schedule s); // Make this static
		void printResult();
		void storeResult();

		void preloadSolution(Schedule s)	{ slvr->preloadSolution(s); }
		void enforceSolution(Schedule s)	{ slvr->enforceSolution(s); }
		void loadSolverConfig(Config c)		{ slvr->loadSolverConfig(c); }

		instance getLoadedInstance() { return loadedConfig.loadedInstance; };
		Schedule solve();	

	private:
		bool	 configLoaded;
		Solver*  slvr;
		Config	 loadedConfig;
		
		ofstream	outputFile;
		Schedule	result;
		double		optVal;
		double		runtime;
		uint64_t	expandedNodes;
		uint64_t	nodesToBest;
		double		objGap;
};

