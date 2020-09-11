#pragma once
#include "generic.h"

// This can be struct?
// Is current machine and currentJobIndex the same?

struct Node {
	uint32_t level;
	JobList pending;
	Schedule schedule;
	uint16_t currentMachine;
	uint16_t currentJobIndex;
	double localUB;
};

/*
class Node {
	public:
		Node(uint32_t level, JobList jobList, Schedule currentSchedule, uint16_t currentMachine, uint16_t currentJobIndex, double localUpperBound);
		~Node();

		bool isFinal()						{ return pending.empty(); };
		int16_t	 getCurrentMachineNo()		{ return currentMachine; };		
		Schedule getSchedule()				{ return schedule; };
		JobList	 getRemainingJobs()			{ return pending;  };
		uint16_t getRemainingJobsCount()	{ return pending.size(); }
		uint16_t getcurrentJobIndex()		{ return currentJobIndex; }
		uint32_t getLevel()					{ return level; };
		double getLocalUpperBound()			{ return localUpperBound; };

		void setCurrentJobIndex(uint16_t newJobIndex) { currentJobIndex = newJobIndex; };

	private:
		uint32_t	level;		
		uint16_t	currentJobIndex;
		int16_t		currentMachine;
		JobList		pending;
		Schedule	schedule;

		double partialSchedule = 1.0;
		double localUpperBound;

};
*/
