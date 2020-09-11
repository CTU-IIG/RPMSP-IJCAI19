#include "Timer.h"

Timer::Timer() {
	isRunning = false;
}

Timer::~Timer() {

}

void Timer::start() {
	startTime = steady_clock::now();
	isRunning = true;
}

void Timer::stop() {
	endTime	 = steady_clock::now();
	timeSpan = duration_cast<duration<double>>(endTime - startTime);

	isRunning = false;
}

double Timer::getElapsed() {
	if (isRunning) {
		middleTime	= steady_clock::now();
		timeSpan	= duration_cast<duration<double>>(middleTime - startTime);
	}

	return timeSpan.count();
}
