#pragma once

#include <math.h>
#include "generic.h"
#include "Solver.h"
#include "statiscticsLib.h"
#include "gurobi_c++.h"

class MILPsolver : public Solver {

	public:
		MILPsolver();
		~MILPsolver();	
		
	protected:
		GRBEnv*		env;
		GRBModel*	model;
};

