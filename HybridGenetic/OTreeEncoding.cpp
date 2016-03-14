/*
 * OTreeEncoding.cpp
 *
 *  Created on: 16 Dec 2015
 *      Author: liang
 */

#include <vector>
#include <limits>
#include <cstdlib>
#include <cstdio>
#include <algorithm>
#include <string>
#include <cmath>

#include "OTreeEncoding.h"

using namespace std;

vector<Kernel> OTreeEncoding::kernels = vector<Kernel>();
int OTreeEncoding::MAX_DEPTH = 32;
int OTreeEncoding::MAX_PHY_DEPTH = 32;
int OTreeEncoding::MAX_HEIGHT = 32;
int OTreeEncoding::MAX_WIDTH = 32;

OTreeEncoding::OTreeEncoding() {
	traversal = vector<bool>();
	sequence = vector<int>();
	implements = vector<int>();
	depth = vector<vector<int> >();
	fitness = std::numeric_limits<double>::max();
	width = 0;
	length = 0;
}

OTreeEncoding::OTreeEncoding(const OTreeEncoding &code) {
	this->traversal = code.traversal;
	this->sequence = code.sequence;
	this->implements = code.implements;
	this->depth = code.depth;
	this->fitness = code.fitness;
	this->width = code.width;
	this->length = code.length;
}

OTreeEncoding::OTreeEncoding(vector<bool> tra, vector<int> seq, vector<int> imp) {
	traversal = tra;
	sequence = seq;
	implements = imp;
	fitness = std::numeric_limits<double>::max();
	width = 0;
	length = 0;
}

void OTreeEncoding::SetFitness(double fit) {
	fitness = fit;
}

double OTreeEncoding::GetFitness() {
	return fitness;
}

void OTreeEncoding::SetKernels(vector<Kernel> list) {
	kernels = list;
}

void OTreeEncoding::InsertSequence(int location, int s) {
	sequence.insert(sequence.begin() + location, s);
}

void OTreeEncoding::RemoveSequenceAt(int location) {
	sequence.erase(sequence.begin() + location);
}

void OTreeEncoding::InsertTraversal(int location, bool b) {
	traversal.insert(traversal.begin() + location, b);
}

void OTreeEncoding::RemoveTraversalAt(int location) {
	traversal.erase(traversal.begin() + location);
}

bool OTreeEncoding::GetTraversalAt(int location) {
	return traversal[location];
}

int OTreeEncoding::GetTraversalCount() {
	return traversal.size();
}

void OTreeEncoding::AddSequence(int s) {
	sequence.push_back(s);
}

void OTreeEncoding::AddImplement(int imp) {
	implements.push_back(imp);
}

int OTreeEncoding::GetImplementAt(int index) {
	return implements[index];
}

void OTreeEncoding::SetImplementAt(int index, int imp) {
	implements[index] = imp;
}

OTreeEncoding OTreeEncoding::GetMutation(int mutImpProb, int mutSeqProb) {
	OTreeEncoding mutation = *this;

	int i = RandU(0, kernels.size());
	int start, length, newStart;

	if (RandU(0, 100) < mutImpProb)
		mutation.implements[i] = RandU(0, kernels[i].GetNumImplements());

	else {
		start = RandU(0, sequence.size());
		length = RandU(0, sequence.size() - start);
		newStart = RandU(0, sequence.size() - length);

		mutation.sequence.erase(mutation.sequence.begin() + start, mutation.sequence.begin() + start + length);
		mutation.sequence.insert(mutation.sequence.begin() + newStart, sequence.begin() + start,
				sequence.begin() + start + length);
	}

	return mutation;
}

void OTreeEncoding::GetPartialTree(vector<int> indecies, vector<int> *sequence, vector<bool> *traversal) {
	*sequence = this->sequence;
	*traversal = this->traversal;

	unsigned i, k, j;
	unsigned zeros, ones;

	for (i = 0; i < sequence->size(); ++i) {
		if (find(indecies.begin(), indecies.end(), sequence->at(i)) != indecies.end()) {
			continue;
		}

		zeros = -1;
		ones = -1;

		for (k = 0; k < traversal->size(); ++k) {
			if (!traversal->at(k))
				zeros++;

			if (zeros == i) {
				traversal->erase(traversal->begin() + k);
				//zeros++; // Prevent further removing zeros
				ones = 0;
				for (j = k; j < traversal->size(); ++j) {
					if (!traversal->at(j))
						ones--;
					else
						ones++;
					if (ones == 1) {
						traversal->erase(traversal->begin() + j);
						break;
					}
				}
				break;
			}

		}

		sequence->erase(sequence->begin() + i);
		i--;
	}
}

double OTreeEncoding::Encoding2Mapping(vector<int> &x, vector<int> &y, ObjectiveType type) {
	unsigned size = kernels.size();
	int i, j, parentIdx, admissableMinX, admissableMinY, admissableMaxX, admissableMaxY;
	bool code;
	int maxW = 0;
	int maxH = 0;
	int w, h, d;
	int m, n;

	int currentWidth, currentLength;
	currentWidth = currentLength = 0;

	vector<vector<int> > depth(currentLength, vector<int>(currentWidth, 0));

	for (i = 0; i < size; ++i) {
		x[i] = -1;
	}

	for (i = 0; i < size; ++i) {
		y[i] = -1;
	}

	parentIdx = -1;

	int sequenceIdx = 0;
	int actualIdx = -1;

	int zeroCount, temp;
	bool parentFound;

	for (i = 0; i < implements.size(); ++i) {
		kernels[i].SetCurrentImplementation(implements[i]);
		if (kernels[i].GetDepth() > MAX_PHY_DEPTH)
			return numeric_limits<double>::max();
	}

	for (i = 0; i < traversal.size(); ++i) {
		code = traversal[i];

		if (!code) {
			w = kernels[sequence[sequenceIdx]].GetWidth();
			h = kernels[sequence[sequenceIdx]].GetHeight();
			d = kernels[sequence[sequenceIdx]].GetDepth();
			actualIdx = sequence[sequenceIdx];
			if (parentIdx == -1)
				x[actualIdx] = 0;
			else
				x[actualIdx] = x[parentIdx] + kernels[parentIdx].GetWidth();

			admissableMinY = admissableMaxY = 0;
			admissableMinX = x[actualIdx];
			admissableMaxX = std::min(x[actualIdx] + w, currentWidth);
			while (admissableMaxY - admissableMinY < h && admissableMaxY < currentLength) {
				for (j = admissableMinX; j < admissableMaxX; ++j) {
					if (d + depth[admissableMaxY][j] > MAX_PHY_DEPTH) {
						admissableMinY = admissableMaxY + 1;
						break;
					}
				}
				admissableMaxY++;
			}
			y[actualIdx] = admissableMinY;

			maxW = std::max(currentWidth, x[actualIdx] + w);
			maxH = std::max(currentLength, y[actualIdx] + h);

			depth.resize(maxH, vector<int>(currentWidth, 0));
			for (n = 0; n < maxH; ++n) {
				depth[n].resize(maxW);
			}

			for (m = 0; m < w; ++m) {
				for (n = 0; n < h; ++n) {
					depth[admissableMinY + n][admissableMinX + m] += d;
				}
			}

			parentIdx = actualIdx;
			sequenceIdx++;
			currentLength = maxH;
			currentWidth = maxW;
		} else {
			zeroCount = 0;
			parentFound = false;

			temp = sequenceIdx - 1;

			for (j = i; j >= 0; j--) {
				if (traversal[j])
					zeroCount--;
				else if (zeroCount < 0) {
					zeroCount++;
					temp--;
				} else {
					parentFound = true;
					break;
				}
			}
			if (parentFound)
				parentIdx = sequence[temp];
			else
				parentIdx = -1;
		}
	}

	width = maxW;
	length = maxH;

	int conflict;
	double utility = 0;

	switch (type) {

	case MinConflict:
		conflict = 0;
		for (i = 0; i < size; ++i) {
			for (j = i + 1; j < size; ++j) {
				if (((x[i] + kernels[i].GetWidth() > x[j] && x[j] >= x[i])
						|| (x[j] + kernels[j].GetWidth() > x[i] && x[i] >= x[j]))
						&& ((y[i] + kernels[i].GetHeight() > y[j] && y[j] >= y[i])
								|| (y[j] + kernels[j].GetHeight() > y[i] && y[i] >= y[j])))
						//Overlapping Condition
						{
					conflict += 1;
				}
			}
		}
		return conflict + sqrt(width * length) + (width + length) / 2.0;
	case MaxUtility:
		if (maxW > MAX_WIDTH || maxH > MAX_HEIGHT)
			return sqrt(width * length) + (width + length) / 2.0;
		else {
			for (i = 0; i < maxH; ++i) {
				for (j = 0; j < maxW; ++j) {
					utility += depth[i][j] <= MAX_PHY_DEPTH ? depth[i][j] : MAX_PHY_DEPTH;
				}
			}
			return 1 - utility / MAX_WIDTH / MAX_HEIGHT / MAX_PHY_DEPTH;
		}
	default:
		return sqrt(width * length) + (width + length) / 2.0;
		//return width * length;
	}
}

int OTreeEncoding::GetWidth() {
	return width;
}

int OTreeEncoding::GetLength() {
	return length;
}

string OTreeEncoding::ToString() {

	string s = "T:\n[";
	char buf[21];
	unsigned i = 0;
	for (i = 0; i < traversal.size(); ++i) {

		if (traversal[i])
			s += " 1 ";
		else
			s += " 0 ";
	}
	s += "] \nS:\n[";
	for (i = 0; i < sequence.size(); ++i) {
		sprintf(buf, " %d", sequence[i]);
		s += buf;
	}

	s += "]\nI:\n[";
	for (i = 0; i < implements.size(); ++i) {
		sprintf(buf, " %d", implements[i]);
		s += buf;
	}
	s += "]";
	return s;

}

bool OTreeEncoding::operator ==(const OTreeEncoding &encode) {
	if (encode.sequence.size() != sequence.size())
		return false;
	if (encode.traversal.size() != traversal.size())
		return false;
	if (encode.implements.size() != implements.size())
		return false;

	if (width == encode.width && length == encode.length)
		return true;

	unsigned i;
	for (i = 0; i < encode.sequence.size(); ++i) {
		if (encode.sequence[i] != sequence[i])
			return false;
	}
	for (i = 0; i < encode.traversal.size(); ++i) {
		if (encode.traversal[i] != traversal[i])
			return false;
	}
	for (i = 0; i < encode.implements.size(); ++i) {
		if (encode.implements[i] != implements[i])
			return false;
	}
	return true;
}
