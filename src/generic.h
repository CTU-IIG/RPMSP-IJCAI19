#pragma once

#include <iostream>
#include <vector>
#include <queue>
#include <stack>

using namespace std;

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
