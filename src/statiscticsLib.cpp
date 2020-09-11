#include "statiscticsLib.h"

// Both functions are yet imprecise!
double probit(double x) {
	return sqrt(2) * errorInverseFunction(2 * x - 1);
}

// Should be precise up to 4th decimal point
double errorInverseFunction(double x) {
	// Taken from https://alohaprog.wordpress.com/category/cc/
	double coef[] = { 1.0, 0.333333333333333, 0.233333333333333, 0.201587301587302, 0.192636684303351, 0.195325476992144, 0.205935864546976, 0.223209757418752, 0.246970233142755, 0.277653825603224, 0.316142623553117,  0.363717587039692 };
	double halfSqrtPi = 0.886226925452758;
	x = x * halfSqrtPi;

	double approximation = (coef[0] * pow(x, 1)) + (coef[2] * pow(x, 5)) +
		(coef[3] * pow(x, 7)) + (coef[4] * pow(x, 9)) + (coef[5] * pow(x, 11)) +
		(coef[6] * pow(x, 13)) + (coef[7] * pow(x, 15)) + (coef[8] * pow(x, 17)) +
		(coef[9] * pow(x, 19)) + (coef[10] * pow(x, 21)) + (coef[11] * pow(x, 23));

	if (x >= 0.95) {
		approximation += (coef[1] * pow(x, 2));
	} else {
		approximation += (coef[1] * pow(x, 3));		
	}

	return approximation;

}

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

uint16_t meanSum(JobList j) {
	uint16_t sum = 0;
	for (Job jb : j) {
		sum += jb.first;
	}
	return sum;
}

uint16_t varSum(JobList j) {
	uint16_t sum = 0;
	for (Job jb : j) {
		sum += jb.second;
	}
	return sum;
}
