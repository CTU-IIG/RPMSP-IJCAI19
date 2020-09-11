#pragma once
#define _USE_MATH_DEFINES

#include "generic.h"
#include <math.h>

typedef vector<double> pattern;

// Move customerServiceLevel function here!

double cumulativeNormalDistribution(double x);
double standardNormalDistribution(double x);

// Following two functions are taken from https://alohaprog.wordpress.com/category/cc/ (with slight changes)
double probit(double x);
double errorInverseFunction(double z);

double customerServiceLevel(Schedule schedule, uint16_t deadline);
double customerServiceLevel(JobList jobList, uint16_t deadline);
double customerServiceLevel(vector<uint16_t> cumulativeMeans, vector<uint16_t> cumulativeVars, uint16_t deadline);
double customerServiceLevel(pattern pat, uint16_t deadline, JobList jobList);