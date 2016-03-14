/*
 * InsertOrder.h
 *
 *  Created on: 16 Dec 2015
 *      Author: liang
 */

#ifndef INSERTORDER_H_
#define INSERTORDER_H_

#include <cstdlib>

#include "OTreeEncoding.h"
#include "OTreeEncoding2D.h"

class InsertOrder {

private:
	vector<int> order;
	vector<int> impls;
	OTreeEncoding encod;
	double fitness;
	int width, length;

public:
	InsertOrder();

	InsertOrder(vector<int> preoder, vector<int> impl);

	InsertOrder(const InsertOrder&);

	InsertOrder& operator=(const InsertOrder& rhs);

	bool operator <(const InsertOrder& o2);

	int GetBestWidth();

	int GetBestLength();

	double GetFitness();

	bool operator ==(const InsertOrder &o);

	vector<int> GetOrder();

	vector<int> GetImplementation();

	InsertOrder Mutation(int mutImpProb, int mutSeqProb);

	InsertOrder CrossOver(InsertOrder p2);

	OTreeEncoding GetGreedyTree(vector<int> &x, vector<int> &y, ObjectiveType type);

	OTreeEncoding2D GetGreedyTree2D(vector<int> &x, vector<int> &y, ObjectiveType type);
};

#endif /* INSERTORDER_H_ */
