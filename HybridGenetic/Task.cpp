/*
 * Task.cpp
 *
 *  Created on: 30 Dec 2015
 *      Author: liang
 */

#include "Task.h"

Task::Task() {
	// TODO Auto-generated constructor stub
	k_index = -1;
	arr_time = -1;
	deadline = -1;
	start_time = -1;
	fin_time = -1;
	priority = -1;
	id = -1;
}

Task::Task(int id, int index, int arrT, int ddl) {
	this->id = id;
	k_index = index;
	arr_time = arrT;
	deadline = ddl;

	start_time = -1;
	fin_time = -1;
	priority = -1;
}

Task::Task(const Task& t) {
	this->id = t.id;
	this->k_index = t.k_index;
	this->arr_time = t.arr_time;
	this->deadline = t.deadline;
	this->start_time = t.start_time;
	this->fin_time = t.fin_time;
	this->priority = t.priority;
}

int Task::getArrTime() const {
	return arr_time;
}

void Task::setArrTime(int arrTime) {
	arr_time = arrTime;
}

int Task::getDeadline() const {
	return deadline;
}

void Task::setDeadline(int deadline) {
	this->deadline = deadline;
}

int Task::getIndex() const {
	return k_index;
}

void Task::setIndex(int index) {
	k_index = index;
}

Task::~Task() {
	// TODO Auto-generated destructor stub
}

double Task::getPriority() const {
	return priority;
}

void Task::setPriority(double priority) {
	this->priority = priority;
}

int Task::getStartTime() const {
	return start_time;
}

void Task::setStartTime(int startTime) {
	start_time = startTime;
}

int Task::getFinTime() const {
	return fin_time;
}

void Task::setFinTime(int finTime) {
	fin_time = finTime;
}

bool Task::operator <(const Task& t) {
	return this->priority < t.priority;
}

int Task::getId() const {
	return id;
}

void Task::setId(int id) {
	this->id = id;
}
