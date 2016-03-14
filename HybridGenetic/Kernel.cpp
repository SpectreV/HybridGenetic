/*
 * kernel.cpp
 *
 *  Created on: 15 Dec 2015
 *      Author: liang
 */

#include "Kernel.h"

#include <vector>
#include <cmath>

using namespace std;

int RandU(int nMin, int nMax) {
	return nMin + (rand() % (nMax - nMin));
}

Kernel::Kernel() {
	widths = vector<int>();
	depths = vector<int>();
	heights = vector<int>();
	delays = vector<double>();
	priority = -1;
	numImp = 0;
	index = -1;
	curImp = -1;
}

Kernel::Kernel(int priority) {
	widths = vector<int>();
	depths = vector<int>();
	heights = vector<int>();
	delays = vector<double>();
	this->priority = priority;
	numImp = 0;
	index = -1;
	curImp = -1;
}

Kernel::Kernel(const Kernel &k) {
	this->priority = k.priority;
	this->widths = k.widths;
	this->heights = k.heights;
	this->depths = k.depths;
	this->delays = k.delays;
	numImp = k.numImp;
	index = k.index;
	curImp = k.curImp;
}

Kernel::Kernel(const Kernel &k, float scale, float variation) {
	if (scale <= 0)
		Kernel(k);
	this->priority = k.priority;
	this->widths = k.widths;
	this->heights = k.heights;
	this->depths = k.depths;
	this->delays = k.delays;
	numImp = k.numImp;
	for (unsigned i = 0; i < numImp; ++i) {
		this->widths[i] = k.widths[i] * scale
				* RandU(100 * (1 - variation), 100 * (1 + variation)) / 100.0;
		this->heights[i] = k.heights[i] * scale
				* RandU(100 * (1 - variation), 100 * (1 + variation)) / 100.0;
		this->depths[i] = k.depths[i] * scale
				* RandU(100 * (1 - variation), 100 * (1 + variation)) / 100.0;
		this->delays[i] = k.delays[i] * scale
				* RandU(100 * (1 - variation), 100 * (1 + variation)) / 100.0;
		if (this->widths[i] < 1)
			this->widths[i] = 1;
		if (this->heights[i] < 1)
			this->heights[i] = 1;
		if (this->depths[i] < 1)
			this->depths[i] = 1;
		if (this->delays[i] < 1)
			this->delays[i] = 1;
	}

	index = -1;
	curImp = -1;
}

bool Kernel::Equals(Kernel kernel) {
	if (priority != kernel.priority || index != kernel.index
			|| numImp != kernel.numImp)
		return false;

	for (int i = 0; i < numImp; ++i) {
		if (widths[i] != kernel.widths[i] || heights[i] != kernel.heights[i]
				|| depths[i] != kernel.depths[i]
				|| delays[i] != kernel.delays[i])
			return false;
	}

	return true;
}

void Kernel::AddImplementation(int w, int h, int d, double delay, int RU_size) {
	widths.push_back((int) ceil(w * 1.0 / RU_size));
	heights.push_back((int) ceil(h * 1.0 / RU_size));
	depths.push_back(d);
	delays.push_back(delay);
	numImp++;
}
