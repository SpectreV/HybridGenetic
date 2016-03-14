/*
 * kernel.h
 *
 *  Created on: 15 Dec 2015
 *      Author: liang
 */

#ifndef KERNEL_H_
#define KERNEL_H_

#include <vector>
#include <cstdlib>

using namespace std;

int RandU(int nMin, int nMax);

class Kernel {
private:
	int priority;
	int index;
	int numImp;
	int curImp;
	vector<int> widths, heights, depths;
	vector<double> delays;
public:
	Kernel();

	Kernel(int priority);

	Kernel(const Kernel&);

	Kernel(const Kernel &k, float scale, float variation);

	bool Equals(Kernel kernel);

	void AddImplementation(int w, int h, int d, double delay, int RU_size);

	inline int GetPriority() {
		return priority;
	}

	inline void SetCurrentImplementation(int current) {
		curImp = current;
	}

	inline void SetIndex(int index) {
		this->index = index;
	}

	inline int GetIndex() {
		return index;
	}

	inline int GetWidth() {
		if (curImp < 0 || curImp >= numImp)
			return 0;
		else
			return widths[curImp];
	}

	inline int GetHeight() {
		if (curImp < 0 || curImp >= numImp)
			return 0;
		else
			return heights[curImp];
	}

	inline int GetDepth() {
		if (curImp < 0 || curImp >= numImp)
			return 0;
		else
			return depths[curImp];
	}

	inline double GetDelay() {
		if (curImp < 0 || curImp >= numImp)
			return 0;
		else
			return delays[curImp];
	}

	inline int GetWidth(int imp) {
		if (imp < 0 || imp >= numImp)
			return 0;
		else
			return widths[imp];
	}

	inline int GetHeight(int imp) {
		if (imp < 0 || imp >= numImp)
			return 0;
		else
			return heights[imp];
	}

	inline int GetDepth(int imp) {
		if (imp < 0 || imp >= numImp)
			return 0;
		else
			return depths[imp];
	}

	inline double GetDelay(int imp) {
		if (imp < 0 || imp >= numImp)
			return 0;
		else
			return delays[imp];
	}

	inline int GetNumImplements() {
		return numImp;
	}

	inline int GetCurrentImplementation() {
		return curImp;
	}
};

#endif /* KERNEL_H_ */
