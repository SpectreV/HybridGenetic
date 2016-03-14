/*
 * Scheduler.h
 *
 *  Created on: 5 Jan 2016
 *      Author: liang
 */

#ifndef SCHEDULER_H_
#define SCHEDULER_H_

#include "Task.h"
#include "Kernel.h"
#include <vector>

enum SchedulePolicy {
	EDF, PATS, CUSTOM
};

enum PlacementPolicy {
	IMMEDIATE, SLACK
};

namespace std {

class Scheduler {
private:
	vector<Task> arr_list, exe_list, fin_list, wait_list;
	vector<vector<int> > conflict, place, depth;
	vector<int> offchip, kernel_in_use, exe_count;
	vector<double> delay, average_exe_time;
	vector<int> wcet, bcet;
	int numKernels, numTasks, maxW, maxH, maxD;
	Scheduler();
	bool Place(int index, bool replace);
	void PrintDepth();
	bool Overlap(int x1, int y1, int w1, int h1, int x2, int y2, int w2, int h2);
public:
	Scheduler(const int numKernels, const int numTasks);
	virtual ~Scheduler();
	void TaskGenerator(const int MaxDepth, const float arrivalRate);
	void Schedule(SchedulePolicy sp, PlacementPolicy pp, bool replace, int debug);

};

}

#endif /* SCHEDULER_H_ */
