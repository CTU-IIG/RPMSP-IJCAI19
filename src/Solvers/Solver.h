#pragma once
#include "generic.h"

struct MILPconfig {
	double objFunctionStep;
	double fracFuncStep;
	uint16_t vMin;
	uint16_t vMax;
};

class Solver {
	public:
		Solver();
		~Solver();

		virtual void setupSolver(instance i) = 0;
		virtual Schedule solve() = 0;

		virtual void preloadSolution(Schedule s);
		virtual void enforceSolution(Schedule s);
		virtual void loadConfig(MILPconfig conf);

		uint64_t getExpandedNodes()	{ return expandedNodes; }
		uint64_t getNodesToBest()	{ return nodesToBest; }
		uint64_t getFracSteps()		{ return fracSteps; }
		uint64_t getObjSteps()		{ return objSteps; }
		double	 getRuntime()		{ return runtime; }
		double	 getGap()			{ return gap; }
		
		static Schedule calculateLB(instance i);

	protected:
		instance loadedInstance;

		uint64_t expandedNodes;
		uint64_t nodesToBest;
		uint64_t fracSteps;
		uint64_t objSteps;
		double	 runtime;	
		double	 gap;
};

