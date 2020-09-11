#include "statiscticsLib.h"

double cumulativeNormalDistribution(double x) {
	// Approximation works only for positive x, so adjustment is needed
	bool xNegative = 0;
	if (x < 0) {
		xNegative = true;
		x = abs(x);
	}

	double b[]	= { 0.2316419, 0.31938153, -0.356563782, 1.781477937, -1.821255978, 1.330274429 };
	double t	= 1 / (1 + b[0]*x);
	double sum	= 0.0;

	for (uint8_t i = 1; i < 6; i++) 
		sum += b[i] * pow(t,i);

	if (!xNegative)
		return 1 - sum * standardNormalDistribution(x);
	else
		return sum * standardNormalDistribution(x);
}

double standardNormalDistribution(double x) {
	x = pow(x, 2);
	return (1/sqrt(2*M_PI)) * exp(-0.5*x);
}

double probit(double x) {
	return M_SQRT2 * errorInverseFunction(2 * x - 1);
}

double errorInverseFunction(double z) {
	double halfSqrtPi	= 0.886226925452758;
	double newZ			= z * halfSqrtPi;
	double result;

	vector<double> coeficients(12, 1);
	coeficients[1]	= 0.333333333333333;
	coeficients[2]	= 0.233333333333333;
	coeficients[3]	= 0.201587301587302;
	coeficients[4]	= 0.192636684303351;
	coeficients[5]	= 0.195325476992144;
	coeficients[6]	= 0.205935864546976;
	coeficients[7]	= 0.223209757418752;
	coeficients[8]	= 0.246970233142755;
	coeficients[9]	= 0.277653825603224;
	coeficients[10] = 0.316142623553117;
	coeficients[11] = 0.363717587039692;

	if (z >= 0.95) {
		result = (coeficients[0] * pow(newZ, 1.00)) +
				(coeficients[1] * pow(newZ, 2.00)) +
				(coeficients[2] * pow(newZ, 5.00)) +
				(coeficients[3] * pow(newZ, 7.00)) +
				(coeficients[4] * pow(newZ, 9.00)) +
				(coeficients[5] * pow(newZ, 11.00)) +
				(coeficients[6] * pow(newZ, 13.00)) +
				(coeficients[7] * pow(newZ, 15.00)) +
				(coeficients[8] * pow(newZ, 17.00)) +
				(coeficients[9] * pow(newZ, 19.00)) +
				(coeficients[10] * pow(newZ, 21.00)) +
				(coeficients[11] * pow(newZ, 23.00));
	} else {
		result = (coeficients[0] * pow(newZ, 1.00)) +
				(coeficients[1] * pow(newZ, 3.00)) +
				(coeficients[2] * pow(newZ, 5.00)) +
				(coeficients[3] * pow(newZ, 7.00)) +
				(coeficients[4] * pow(newZ, 9.00)) +
				(coeficients[5] * pow(newZ, 11.00)) +
				(coeficients[6] * pow(newZ, 13.00)) +
				(coeficients[7] * pow(newZ, 15.00)) +
				(coeficients[8] * pow(newZ, 17.00)) +
				(coeficients[9] * pow(newZ, 19.00)) +
				(coeficients[10] * pow(newZ, 21.00)) +
				(coeficients[11] * pow(newZ, 23.00));
	}

	return result;
}

double customerServiceLevel(Schedule schedule, uint16_t deadline) {
	if (schedule.empty()) {
		return 1.0;
	}

	uint16_t numOfMachines = schedule.size();
	double	 x, lowerBound = 1.0;
	vector<uint16_t> cumulativeMeans(numOfMachines, 0);
	vector<uint16_t> cumulativeVars(numOfMachines, 0);

	for (uint16_t m = 0; m < schedule.size(); m++) {
		for (Job j : schedule[m]) {
			cumulativeMeans[m] += j.first;
			cumulativeVars[m] += j.second;
		}

		x = ((int32_t(deadline) - cumulativeMeans[m]) / sqrt(cumulativeVars[m]));
		lowerBound *= cumulativeNormalDistribution(x);
	}
	return lowerBound;
}

// Checked
double customerServiceLevel(JobList jobList, uint16_t deadline) {
	uint16_t cumulativeMean = 0;
	uint16_t cumulativeVar = 0;

	for (Job j : jobList) {
		cumulativeMean += j.first;
		cumulativeVar += j.second;
	}
	return cumulativeNormalDistribution((int32_t(deadline) - cumulativeMean) / sqrt(cumulativeVar));
}

// Checked
double customerServiceLevel(vector<uint16_t> cumulativeMeans, vector<uint16_t> cumulativeVars, uint16_t deadline) {
	uint16_t numOfMachines = cumulativeMeans.size();
	double x, lowerBound = 1.0;

	for (uint16_t k = 0; k < numOfMachines; k++) {
		x = (int32_t(deadline) - cumulativeMeans[k]) / sqrt(cumulativeVars[k]);
		lowerBound *= cumulativeNormalDistribution(x);
	}

	return lowerBound;
}

double customerServiceLevel(pattern pat, uint16_t deadline, JobList jobList) {
	JobList jl;

	for (uint16_t i = 0; i < pat.size(); i++) {
		if (pat[i] == 1.0) {
			jl.push_back(jobList[i]);
		}
	}

	return customerServiceLevel(jl, deadline);
}
