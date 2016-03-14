/*
 * InsertOrder.cpp
 *
 *  Created on: 16 Dec 2015
 *      Author: liang
 */

#include <cstdlib>
#include <vector>

#include "InsertOrder.h"
#include "OTreeEncoding.h"
#include "OTreeEncoding2D.h"

using namespace std;

InsertOrder::InsertOrder() {
	order = vector<int>();
	impls = vector<int>();
	fitness = 0;
	width = 0;
	length = 0;
}

InsertOrder::InsertOrder(vector<int> preoder, vector<int> impl) {
	order = preoder;
	impls = impl;
	fitness = 0;
	width = 0;
	length = 0;

}

InsertOrder::InsertOrder(const InsertOrder &o) {
	order = o.order;
	impls = o.impls;
	fitness = o.fitness;
	width = o.width;
	length = o.length;
}

InsertOrder& InsertOrder::operator=(const InsertOrder& rhs) {
	this->order = rhs.order;
	this->impls = rhs.impls;
	this->fitness = rhs.fitness;
	this->width = rhs.width;
	this->length = rhs.length;
	return *this;
}

bool InsertOrder::operator <(const InsertOrder& o2) {
	return fitness < o2.fitness;
}

int InsertOrder::GetBestWidth() {
	return width;
}

int InsertOrder::GetBestLength() {
	return length;
}

double InsertOrder::GetFitness() {
	return fitness;
}

bool InsertOrder::operator ==(const InsertOrder &o) {
	if (o.order.size() != order.size())
		return false;
	if (o.impls.size() != impls.size())
		return false;
	for (int i = 0; i < order.size(); ++i) {
		if (order[i] != o.order[i])
			return false;
		if (impls[i] != o.impls[i])
			return false;
	}
	return true;
}

vector<int> InsertOrder::GetOrder() {
	return order;
}

vector<int> InsertOrder::GetImplementation() {
	return impls;
}

InsertOrder InsertOrder::Mutation(int mutImpProb, int mutSeqProb) {
	InsertOrder mutation = *this;
	for (unsigned i = 0; i < order.size(); ++i) {
		if (RandU(0, 100) < mutImpProb)
			do {
				mutation.impls[i] = (RandU(0, OTreeEncoding::kernels[i].GetNumImplements()));
				OTreeEncoding::kernels[i].SetCurrentImplementation(mutation.impls[i]);
			} while (OTreeEncoding::kernels[i].GetDepth() > OTreeEncoding::MAX_DEPTH);

	}

	int start, length, newStart;
	if (RandU(0, 100) < mutSeqProb) {
		start = RandU(0, order.size());
		length = RandU(0, order.size() - start);
		newStart = RandU(0, order.size() - length);

		mutation.order.erase(mutation.order.begin() + start, mutation.order.begin() + start + length);
		mutation.order.insert(mutation.order.begin() + newStart, order.begin() + start, order.begin() + start + length);
	}

	return mutation;
}

InsertOrder InsertOrder::CrossOver(InsertOrder p2) {
	vector<int> order_p1(order);
	vector<int> order_p2(p2.GetOrder());
	vector<int> impl;
	vector<int> indecies;
	int i;

	for (i = 0; i < order.size(); ++i) {
		if (RandU(0, 100) > 50)
			indecies.push_back(i);
		impl.push_back(0);
	}

	for (i = 0; i < indecies.size(); ++i) {
		order_p2.erase(find(order_p2.begin(), order_p2.end(), indecies[i]));
	}

	for (i = 0; i < order_p1.size(); ++i) {
		if (find(indecies.begin(), indecies.end(), order_p1[i]) == indecies.end()) {
			order_p1.erase(order_p1.begin() + i);
			i--;
		}
	}

	order_p1.insert(order_p1.end(), order_p2.begin(), order_p2.end());

	for (i = 0; i < order.size(); ++i) {
		if (find(indecies.begin(), indecies.end(), i) != indecies.end())
			impl[i] = this->impls[i];
		else
			impl[i] = p2.impls[i];
	}

	if (order_p1.size() == order.size() && impl.size() == impls.size()) {
		InsertOrder o(order_p1, impl);
		return o;
	} else
		exit(-1);
}

OTreeEncoding InsertOrder::GetGreedyTree(vector<int> &x, vector<int> &y, ObjectiveType type) {

	int i;
	int location;
	double cost, minCost;
	OTreeEncoding greed;

	//depth = new int[1, 1];
	/*int[] x, y;
	 x = new int[order.size()];
	 y = new int[order.size()];*/

	for (i = 0; i < order.size(); ++i) {
		greed.AddImplement(-1);
	}

	int bestLoc = -1;
	int seq, bestSeq, bestImp;

	for (i = 0; i < order.size(); ++i) {

		bestLoc = -1;
		bestSeq = -1;
		minCost = std::numeric_limits<double>::max();

		greed.SetImplementAt(order[i], impls[i]);

		seq = 0;

		for (location = 0; location <= greed.GetTraversalCount(); ++location) {
			greed.InsertSequence(seq, order[i]);
			greed.InsertTraversal(location, true);
			greed.InsertTraversal(location, false);

			cost = greed.Encoding2Mapping(x, y, type);
			if (cost < minCost) {
				minCost = cost;
				bestLoc = location;
				bestSeq = seq;
			}

			greed.RemoveTraversalAt(location);
			greed.RemoveTraversalAt(location);
			greed.RemoveSequenceAt(seq);

			if (location < greed.GetTraversalCount() && !greed.GetTraversalAt(location))
				seq++;

		}

		greed.InsertTraversal(bestLoc, true);
		greed.InsertTraversal(bestLoc, false);
		greed.InsertSequence(bestSeq, order[i]);
	}
	fitness = greed.Encoding2Mapping(x, y, type);
	width = greed.GetWidth();
	length = greed.GetLength();
	return greed;
}

OTreeEncoding2D InsertOrder::GetGreedyTree2D(vector<int> &x, vector<int> &y, ObjectiveType type) {

	int i;
	int locationH, locationV;
	double cost, minCost;
	OTreeEncoding2D greed;

	//depth = new int[1, 1];
	/*int[] x, y;
	 x = new int[order.size()];
	 y = new int[order.size()];*/

	for (i = 0; i < order.size(); ++i) {
		greed.AddImplement(-1);
	}

	int bestLocH, bestLocV;
	int seqH, seqV, bestSeqH, bestSeqV, bestImp;

	for (i = 0; i < order.size(); ++i) {

		bestLocH = bestLocV = -1;
		bestSeqH = bestSeqV = -1;
		minCost = std::numeric_limits<double>::max();

		greed.SetImplementAt(order[i], impls[i]);

		seqH = 0;

		for (locationH = 0; locationH <= greed.GetTraversalCount(H); ++locationH) {
			greed.InsertSequence(seqH, order[i], H);
			greed.InsertTraversal(locationH, true, H);
			greed.InsertTraversal(locationH, false, H);

			seqV = 0;

			for (locationV = 0; locationV <= greed.GetTraversalCount(V); ++locationV) {
				greed.InsertSequence(seqV, order[i], V);
				greed.InsertTraversal(locationV, true, V);
				greed.InsertTraversal(locationV, false, V);

				cost = greed.Encoding2Mapping(x, y, type);
				if (cost < minCost) {
					minCost = cost;
					bestLocH = locationH;
					bestLocV = locationV;
					bestSeqH = seqH;
					bestSeqV = seqV;
				}

				greed.RemoveTraversalAt(locationV, V);
				greed.RemoveTraversalAt(locationV, V);
				greed.RemoveSequenceAt(seqV, V);

				if (locationV < greed.GetTraversalCount(V) && !greed.GetTraversalAt(locationV, V))
					seqV++;

			}

			greed.RemoveTraversalAt(locationH, H);
			greed.RemoveTraversalAt(locationH, H);
			greed.RemoveSequenceAt(seqH, H);

			if (locationH < greed.GetTraversalCount(H) && !greed.GetTraversalAt(locationH, H))
				seqH++;

		}

		greed.InsertTraversal(bestLocH, true, H);
		greed.InsertTraversal(bestLocH, false, H);
		greed.InsertTraversal(bestLocV, true, V);
		greed.InsertTraversal(bestLocV, false, V);
		greed.InsertSequence(bestSeqH, order[i], H);
		greed.InsertSequence(bestSeqV, order[i], V);
	}
	fitness = greed.Encoding2Mapping(x, y, type);
	width = greed.GetWidth();
	length = greed.GetLength();
	return greed;
}

