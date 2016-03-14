/*
 * OTreeEncoding.h
 *
 *  Created on: 15 Dec 2015
 *      Author: liang
 */

#ifndef OTREEENCODING2D_H_
#define OTREEENCODING2D_H_

#include <vector>
#include <limits>
#include <cstdlib>
#include <algorithm>
#include <string>

#include "Kernel.h"
#include "OTreeEncoding.h"

using namespace std;

enum Direction {
	H, V
};

class OTreeEncoding2D: public OTreeEncoding {

private:
	vector<bool> traversalH;
	vector<bool> traversalV;
	vector<int> sequenceH;
	vector<int> sequenceV;

	/*
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
	 */
public:

	OTreeEncoding2D();

	OTreeEncoding2D(const OTreeEncoding2D&);

	OTreeEncoding2D(vector<bool> traH, vector<bool> traV, vector<int> seqH, vector<int> seqV, vector<int> imp);

	void SetFitness(double fit);

	double GetFitness();

	static void SetKernels(vector<Kernel> list);

	void InsertSequence(int location, int s, Direction d);

	void RemoveSequenceAt(int location, Direction d);

	void InsertTraversal(int location, bool b, Direction d);

	void RemoveTraversalAt(int location, Direction d);

	bool GetTraversalAt(int location, Direction d);

	int GetTraversalCount(Direction d);

	void AddSequence(int s, Direction d);

	void AddImplement(int imp);

	int GetImplementAt(int index);

	void SetImplementAt(int index, int imp);

	OTreeEncoding2D GetMutation(int mutImpProb, int mutSeqProb);

	void GetPartialTree(vector<int> indecies, vector<int> *sequence, vector<bool> *traversal, Direction d);

	double Encoding2Mapping(vector<int> &x, vector<int> &y, ObjectiveType type = MinArea);

	bool operator <(const OTreeEncoding2D& encode) {
		return this->fitness < encode.fitness;
	}

	bool operator ==(const OTreeEncoding2D &encode);

	int GetWidth();

	int GetLength();

	string ToString();

};

#endif /* OTREEENCODING_H_ */
