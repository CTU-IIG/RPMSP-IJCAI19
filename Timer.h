#pragma once
#include "generic.h"
#include <chrono>

using namespace std::chrono;

class Timer {

	public:
		Timer();
		~Timer();

		void start();
		void stop();
		
		double getElapsed();

	private:
		bool isRunning;

		steady_clock::time_point startTime;
		steady_clock::time_point middleTime;
		steady_clock::time_point endTime;
		duration<double> timeSpan;
};

