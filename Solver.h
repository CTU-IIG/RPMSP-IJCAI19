#pragma once
#include "generic.h"
#include "Config.h"
#include "Timer.h"
#include "statiscticsLib.h"

class Solver {
	public:
		Solver();
		~Solver();

		// Remove setup solver function and move everything to constructor
		virtual Schedule solve() = 0;		

		// This functions do nothing unless reimplemented, so maybe make them pure virutal since I will only need them to do nothing in BB?
		virtual void loadSolverConfig(Config c);
		virtual void preloadSolution(Schedule s);
		virtual void enforceSolution(Schedule s);

		uint64_t getExpandedNodes()	{ return expandedNodes; }
		uint64_t getNodesToBest()	{ return nodesToBest; }
		double	 getRuntime()		{ return runtime; }
		double	 getObjGap()		{ return objGap; }
		
		static Schedule calculateLB(instance i);
		
		uint16_t calculate_vMin();
		uint16_t calculate_vMax();

	protected:
		instance loadedInstance;

		Timer timer;
		double timeout; // Set Time Limit on solver execution

		uint64_t expandedNodes;
		uint64_t nodesToBest;
		double	 runtime;		
		double	 objGap;
};

