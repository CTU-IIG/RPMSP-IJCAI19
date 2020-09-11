#pragma once
#include "generic.h"
#include "MILPsolver.h"

class MILPsolver_secondary : public MILPsolver {

public:
	MILPsolver_secondary();
	~MILPsolver_secondary();

	void setupSolver(instance i);
	Schedule solve();
	
};