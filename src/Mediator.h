#pragma once
#include "generic.h"
#include "MILPsolver_primary.h"
#include "MILPsolver_secondary.h"

class Mediator {
	public:
		Mediator(string solverType);
		Mediator(string instanceFileName, string solverType);
		~Mediator();

		void loadInstance(string fileName);
		void printSchedule(Schedule s);
		void printResult();
		void storeResult();

		void preloadSolution(Schedule s) { slvr->preloadSolution(s); }
		void enforceSolution(Schedule s) { slvr->enforceSolution(s); }
		void loadConfig(MILPconfig c)	 { slvr->loadConfig(c); }

		instance getLoadedInstance() { return currentInstace; }
		Schedule solve();	

	private:
		bool	 instanceLoaded = false;
		Solver*  slvr;
		instance currentInstace;
		
		ofstream	outputFile;
		Schedule	result;
		double		optVal;
		double		runtime;
		uint64_t	expandedNodes;
		uint64_t	nodesToBest;
		uint64_t	fracSteps;
		uint64_t	objSteps;
};

