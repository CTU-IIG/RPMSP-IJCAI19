#pragma once
#define _USE_MATH_DEFINES

#include "generic.h"
#include <math.h>

// Move customerServiceLevel function here!

double cumulativeNormalDistribution(double x);
double standardNormalDistribution(double x);

double probit(double x); // Quanitle function of normal distribution
double errorInverseFunction(double x);

double customerServiceLevel(Schedule schedule, uint16_t deadline);
double customerServiceLevel(JobList jobList, uint16_t deadline);
double customerServiceLevel(vector<uint16_t> cumulativeMeans, vector<uint16_t> cumulativeVars, uint16_t deadline);

uint16_t meanSum(JobList j);
uint16_t varSum(JobList j);