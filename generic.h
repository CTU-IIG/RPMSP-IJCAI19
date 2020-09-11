#pragma once

#include <ctime>
#include <time.h>
#include <iostream>
#include <string>
#include <map>
#include <vector>
#include <queue>
#include <stack>
#include <fstream>
#include <limits>
#include <cfloat>
#include <algorithm>

using namespace std;

// Reformat later to mean and variance
struct Job {
	uint16_t first;
	uint16_t second;
	uint16_t id;
};

typedef vector<Job>		JobList;
typedef vector<JobList>	Schedule;

struct instance {
	string	 fileName;
	uint16_t numOfJobs;
	uint16_t numOfMachines;
	uint16_t deadline;
	JobList	 jobs;
};
