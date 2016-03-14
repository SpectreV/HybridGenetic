/*
 * OTreeEncoding.h
 *
 *  Created on: 15 Dec 2015
 *      Author: liang
 */

#ifndef OTREEENCODING_H_
#define OTREEENCODING_H_

#include <vector>
#include <limits>
#include <cstdlib>
#include <algorithm>
#include <string>

#include "Kernel.h"

using namespace std;

enum ObjectiveType {
	MinArea, MaxUtility, MinConflict
};

class OTreeEncoding {

private:
	vector<bool> traversal;
	vector<int> sequence;

protected:
	vector<int> implements;
	int width;
	int length;
	vector<vector<int> > depth;
	int x[], y[];
	double fitness;

public:
	static int MAX_DEPTH, MAX_WIDTH, MAX_HEIGHT, MAX_PHY_DEPTH;
	static vector<Kernel> kernels;

	OTreeEncoding();

	OTreeEncoding(const OTreeEncoding&);

	OTreeEncoding(vector<bool> tra, vector<int> seq, vector<int> imp);

	void SetFitness(double fit);

	double GetFitness();

	static void SetKernels(vector<Kernel> list);

	void InsertSequence(int location, int s);

	void RemoveSequenceAt(int location);

	void InsertTraversal(int location, bool b);

	void RemoveTraversalAt(int location);

	bool GetTraversalAt(int location);

	int GetTraversalCount();

	void AddSequence(int s);

	void AddImplement(int imp);

	int GetImplementAt(int index);

	void SetImplementAt(int index, int imp);

	OTreeEncoding GetMutation(int mutImpProb, int mutSeqProb);

	void GetPartialTree(vector<int> indecies, vector<int> *sequence, vector<bool> *traversal);

	double Encoding2Mapping(vector<int> &x, vector<int> &y, ObjectiveType type = MinArea);

	bool operator <(const OTreeEncoding& encode) {
		return this->fitness < encode.fitness;
	}

	bool operator ==(const OTreeEncoding &encode);

	int GetWidth();

	int GetLength();

	string ToString();

};

#endif /* OTREEENCODING_H_ */
