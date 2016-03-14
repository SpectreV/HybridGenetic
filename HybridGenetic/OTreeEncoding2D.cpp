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

#include "OTreeEncoding2D.h"

using namespace std;

OTreeEncoding2D::OTreeEncoding2D() {
	traversalH = vector<bool>();
	traversalV = vector<bool>();
	sequenceH = vector<int>();
	sequenceV = vector<int>();
	implements = vector<int>();
	depth = vector<vector<int> >();
	fitness = std::numeric_limits<double>::max();
	width = 0;
	length = 0;
}

OTreeEncoding2D::OTreeEncoding2D(const OTreeEncoding2D &code) {
	this->traversalH = code.traversalH;
	this->traversalV = code.traversalV;
	this->sequenceH = code.sequenceH;
	this->sequenceV = code.sequenceV;
	this->implements = code.implements;
	this->depth = code.depth;
	this->fitness = code.fitness;
	this->width = code.width;
	this->length = code.length;
}

OTreeEncoding2D::OTreeEncoding2D(vector<bool> traH, vector<bool> traV, vector<int> seqH, vector<int> seqV,
		vector<int> imp) {
	traversalH = traH;
	sequenceH = seqH;
	traversalV = traV;
	sequenceV = seqV;
	implements = imp;
	fitness = std::numeric_limits<double>::max();
	width = 0;
	length = 0;
}

void OTreeEncoding2D::SetFitness(double fit) {
	fitness = fit;
}

double OTreeEncoding2D::GetFitness() {
	return fitness;
}

void OTreeEncoding2D::SetKernels(vector<Kernel> list) {
	kernels = list;
}

void OTreeEncoding2D::InsertSequence(int location, int s, Direction d) {
	if (d == H)
		sequenceH.insert(sequenceH.begin() + location, s);
	else
		sequenceV.insert(sequenceV.begin() + location, s);
}

void OTreeEncoding2D::RemoveSequenceAt(int location, Direction d) {
	if (d == H)
		sequenceH.erase(sequenceH.begin() + location);
	else
		sequenceV.erase(sequenceV.begin() + location);
}

void OTreeEncoding2D::InsertTraversal(int location, bool b, Direction d) {
	if (d == H)
		traversalH.insert(traversalH.begin() + location, b);
	else
		traversalV.insert(traversalV.begin() + location, b);
}

void OTreeEncoding2D::RemoveTraversalAt(int location, Direction d) {
	if (d == H)
		traversalH.erase(traversalH.begin() + location);
	else
		traversalV.erase(traversalV.begin() + location);
}

bool OTreeEncoding2D::GetTraversalAt(int location, Direction d) {
	if (d == H)
		return traversalH[location];
	else
		return traversalV[location];
}

int OTreeEncoding2D::GetTraversalCount(Direction d) {
	if (d == H)
		return traversalH.size();
	else
		return traversalV.size();
}

void OTreeEncoding2D::AddSequence(int s, Direction d) {
	if (d == H)
		sequenceH.push_back(s);
	else
		sequenceV.push_back(s);
}

void OTreeEncoding2D::AddImplement(int imp) {
	implements.push_back(imp);
}

int OTreeEncoding2D::GetImplementAt(int index) {
	return implements[index];
}

void OTreeEncoding2D::SetImplementAt(int index, int imp) {
	implements[index] = imp;
}

OTreeEncoding2D OTreeEncoding2D::GetMutation(int mutImpProb, int mutSeqProb) {
	OTreeEncoding2D mutation = *this;

	int i = RandU(0, kernels.size());

	if (RandU(0, 100) < mutImpProb)
		mutation.implements[i] = RandU(0, kernels[i].GetNumImplements());

	int start, length, newStart;
	if (RandU(0, 100) < mutSeqProb) {
		start = RandU(0, sequenceH.size());
		length = RandU(0, sequenceH.size() - start);
		newStart = RandU(0, sequenceH.size() - length);

		mutation.sequenceH.erase(mutation.sequenceH.begin() + start, mutation.sequenceH.begin() + start + length);
		mutation.sequenceH.insert(mutation.sequenceH.begin() + newStart, sequenceH.begin() + start,
				sequenceH.begin() + start + length);

		start = RandU(0, sequenceV.size());
		length = RandU(0, sequenceV.size() - start);
		newStart = RandU(0, sequenceV.size() - length);

		mutation.sequenceV.erase(mutation.sequenceV.begin() + start, mutation.sequenceV.begin() + start + length);
		mutation.sequenceV.insert(mutation.sequenceV.begin() + newStart, sequenceV.begin() + start,
				sequenceV.begin() + start + length);
	}

	return mutation;
}

void OTreeEncoding2D::GetPartialTree(vector<int> indecies, vector<int> *sequence, vector<bool> *traversal,
		Direction d) {
	if (d == H) {
		*sequence = this->sequenceH;
		*traversal = this->traversalH;
	} else {
		*sequence = this->sequenceV;
		*traversal = this->traversalV;
	}
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

double OTreeEncoding2D::Encoding2Mapping(vector<int> &x, vector<int> &y, ObjectiveType type) {
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

	int zeroCount, temp;
	bool parentFound;

	for (i = 0; i < implements.size(); ++i) {
		kernels[i].SetCurrentImplementation(implements[i]);
		if (kernels[i].GetDepth() > MAX_PHY_DEPTH)
			return numeric_limits<double>::max();
	}

	int sequenceIdx = 0;
	int actualIdx = -1;

	for (i = 0; i < traversalH.size(); ++i) {
		code = traversalH[i];

		if (!code) {
			w = kernels[sequenceH[sequenceIdx]].GetWidth();
//			h = kernels[sequence[sequenceIdx]].GetHeight();
//			d = kernels[sequence[sequenceIdx]].GetDepth();
			actualIdx = sequenceH[sequenceIdx];
			if (parentIdx == -1)
				x[actualIdx] = 0;
			else
				x[actualIdx] = x[parentIdx] + kernels[parentIdx].GetWidth();

			/*admissableMinY = admissableMaxY = 0;
			 admissableMinX = x[actualIdx];
			 admissableMaxX = std::min(x[actualIdx] + w, currentWidth);
			 while (admissableMaxY - admissableMinY < h && admissableMaxY < currentLength) {
			 for (j = admissableMinX; j < admissableMaxX; ++j) {
			 if (d + depth[admissableMaxY][j] > MAX_DEPTH) {
			 admissableMinY = admissableMaxY + 1;
			 break;
			 }
			 }
			 admissableMaxY++;
			 }
			 y[actualIdx] = admissableMinY;*/

			maxW = std::max(currentWidth, x[actualIdx] + w);
			/*maxH = std::max(currentLength, y[actualIdx] + h);

			 depth.resize(maxH, vector<int>(currentWidth, 0));
			 for (n = 0; n < maxH; ++n) {
			 depth[n].resize(maxW);
			 }

			 for (m = 0; m < w; ++m) {
			 for (n = 0; n < h; ++n) {
			 depth[admissableMinY + n][admissableMinX + m] += d;
			 }
			 }*/

			parentIdx = actualIdx;
			sequenceIdx++;
//			currentLength = maxH;
			currentWidth = maxW;
		} else {
			zeroCount = 0;
			parentFound = false;

			temp = sequenceIdx - 1;

			for (j = i; j >= 0; j--) {
				if (traversalH[j])
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
				parentIdx = sequenceH[temp];
			else
				parentIdx = -1;
		}
	}

	sequenceIdx = 0;
	actualIdx = -1;

	for (i = 0; i < traversalV.size(); ++i) {
		code = traversalV[i];

		if (!code) {
//			w = kernels[sequence[sequenceIdx]].GetWidth();
			h = kernels[sequenceV[sequenceIdx]].GetHeight();
//			d = kernels[sequence[sequenceIdx]].GetDepth();
			actualIdx = sequenceV[sequenceIdx];
			if (parentIdx == -1)
				y[actualIdx] = 0;
			else
				y[actualIdx] = x[parentIdx] + kernels[parentIdx].GetHeight();

			/*admissableMinY = admissableMaxY = 0;
			 admissableMinX = x[actualIdx];
			 admissableMaxX = std::min(x[actualIdx] + w, currentWidth);
			 while (admissableMaxY - admissableMinY < h && admissableMaxY < currentLength) {
			 for (j = admissableMinX; j < admissableMaxX; ++j) {
			 if (d + depth[admissableMaxY][j] > MAX_DEPTH) {
			 admissableMinY = admissableMaxY + 1;
			 break;
			 }
			 }
			 admissableMaxY++;
			 }
			 y[actualIdx] = admissableMinY;*/

//			maxW = std::max(currentWidth, x[actualIdx] + w);
			maxH = std::max(currentLength, y[actualIdx] + h);

			/*depth.resize(maxH, vector<int>(currentWidth, 0));
			 for (n = 0; n < maxH; ++n) {
			 depth[n].resize(maxW);
			 }

			 for (m = 0; m < w; ++m) {
			 for (n = 0; n < h; ++n) {
			 depth[admissableMinY + n][admissableMinX + m] += d;
			 }
			 }*/

			parentIdx = actualIdx;
			sequenceIdx++;
			currentLength = maxH;
			//currentWidth = maxW;
		} else {
			zeroCount = 0;
			parentFound = false;

			temp = sequenceIdx - 1;

			for (j = i; j >= 0; j--) {
				if (traversalV[j])
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
				parentIdx = sequenceV[temp];
			else
				parentIdx = -1;
		}
	}

	width = maxW;
	length = maxH;

	depth.resize(maxH, vector<int>(currentWidth, 0));
	for (n = 0; n < maxH; ++n) {
		depth[n].resize(maxW);
	}

	for (i = 0; i < size; ++i) {
		w = kernels[i].GetWidth();
		h = kernels[i].GetHeight();
		d = kernels[i].GetDepth();
		for (m = 0; m < w; ++m) {
			for (n = 0; n < h; ++n) {
				if (depth[y[i] + n][x[i] + m] + d <= MAX_PHY_DEPTH)
					depth[y[i] + n][x[i] + m] += d;
				else
					return numeric_limits<double>::max();
			}
		}
	}

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
		//return sqrt(width * length) + (width + length) / 2.0;
		return width * length;
	}
}

int OTreeEncoding2D::GetWidth() {
	return width;
}

int OTreeEncoding2D::GetLength() {
	return length;
}

string OTreeEncoding2D::ToString() {

	string s = "T:\n[";
	char buf[21];
	unsigned i = 0;
	for (i = 0; i < traversalH.size(); ++i) {

		if (traversalH[i])
			s += "1";
		else
			s += "0";
	}
	s += ",";
	for (i = 0; i < traversalV.size(); ++i) {

		if (traversalV[i])
			s += "1";
		else
			s += "0";
	}
	s += "] \nS:\n[";
	for (i = 0; i < sequenceH.size(); ++i) {
		sprintf(buf, " %d", sequenceH[i]);
		s += buf;
	}
	s += ",";
	for (i = 0; i < sequenceV.size(); ++i) {
		sprintf(buf, " %d", sequenceV[i]);
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

bool OTreeEncoding2D::operator ==(const OTreeEncoding2D &encode) {
	if (encode.sequenceH.size() != sequenceH.size())
		return false;
	if (encode.traversalH.size() != traversalH.size())
		return false;
	if (encode.sequenceV.size() != sequenceV.size())
		return false;
	if (encode.traversalV.size() != traversalV.size())
		return false;
	if (encode.implements.size() != implements.size())
		return false;

	if(width == encode.width && length == encode.length)
			return true;

	unsigned i;
	for (i = 0; i < encode.sequenceH.size(); ++i) {
		if (encode.sequenceH[i] != sequenceH[i])
			return false;
	}
	for (i = 0; i < encode.traversalH.size(); ++i) {
		if (encode.traversalH[i] != traversalH[i])
			return false;
	}
	for (i = 0; i < encode.sequenceV.size(); ++i) {
		if (encode.sequenceV[i] != sequenceV[i])
			return false;
	}
	for (i = 0; i < encode.traversalV.size(); ++i) {
		if (encode.traversalV[i] != traversalV[i])
			return false;
	}
	for (i = 0; i < encode.implements.size(); ++i) {
		if (encode.implements[i] != implements[i])
			return false;
	}

	return true;
}
