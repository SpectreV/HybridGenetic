/*
 * Task.h
 *
 *  Created on: 30 Dec 2015
 *      Author: liang
 */

#ifndef TASK_H_
#define TASK_H_

class Task {
private:
	int k_index;
	int arr_time;
	int deadline;

	int start_time;
	int fin_time;
	double priority;
	int id;
public:
	Task();
	Task(int id, int index, int arrT, int ddl);
	Task(const Task& t);
	virtual ~Task();
	int getArrTime() const;
	void setArrTime(int arrTime);
	int getDeadline() const;
	void setDeadline(int deadline);
	int getIndex() const;
	void setIndex(int index);
	double getPriority() const;
	void setPriority(double priority);
	bool operator <(const Task& t);
	int getStartTime() const;
	void setStartTime(int exeTime);
	int getFinTime() const;
	void setFinTime(int finTime);
	int getId() const;
	void setId(int id);
};

#endif /* TASK_H_ */
