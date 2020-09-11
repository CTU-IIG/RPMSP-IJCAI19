#pragma once
#include "generic.h"
#include "MILPsolver.h"
#include "config_MILPprimary.h"

class MILPsolver_primary : public MILPsolver {

	public:
		MILPsolver_primary(config_MILPprimary c, vector<double> prices = vector<double>() ); // Second parameter used for branch and price
		~MILPsolver_primary();

		void loadSolverConfig(config_MILPprimary c);
		Schedule solve();		

		void preloadSolution(Schedule s);
		void enforceSolution(Schedule s);

		GRBModel* getModel()				{ return model; }
		GRBVar*	  getAssignmentVariable()	{ return w[0]; }
		void buildObjective();

	private:
		vector<double> objFunctionApproxPoints; // Used in objective function approximation
		vector<double> fracApproxIntervalBounds; // Used in fraction approximation

		vector<double> approxSteps_fraction();
		vector<double> approxSteps_objFunction();
		
		void buildModel(vector<double> prices);
		void buildConstraints();

		void flush();

		config_MILPprimary conf;
		
		// Model Variables
		GRBVar* v, *u, *t;
		GRBVar** z, **s, **y, **w;

		double bigM_1, bigM_2; // Make function to init those two, they are passed as reference so they have to be here
};

