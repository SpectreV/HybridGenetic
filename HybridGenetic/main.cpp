/*
 * main.cpp
 *
 *  Created on: 16 Dec 2015
 *      Author: liang
 */

#include <iostream>
#include <fstream>
#include <ctime>
#include <string>
#include <algorithm>
#include <omp.h>
#include <random>

#include "Kernel.h"
#include "OTreeEncoding.h"
#include "OTreeEncoding2D.h"
#include "InsertOrder.h"
#include "Task.h"
#include "Scheduler.h"

using namespace std;

class GeneticAlgorithm {
private:
	vector<Kernel> kernels;
	vector<OTreeEncoding> population;
	vector<OTreeEncoding2D> population2D;
	vector<double> fitness;
public:
	GeneticAlgorithm(int numKernels, int RU_size) {
		kernels = vector<Kernel>();

		Kernel ARF(RandU(0, 10));
		ARF.AddImplementation(15, 15, 25, 64, RU_size);
		ARF.AddImplementation(16, 16, 13, 32.11, RU_size);
		ARF.AddImplementation(23, 23, 7, 20.79, RU_size);
		kernels.push_back(ARF);

		Kernel DIFFEQ(RandU(1, 10));
		DIFFEQ.AddImplementation(6, 6, 32, 54.72, RU_size);
		DIFFEQ.AddImplementation(7, 7, 17, 32.81, RU_size);
		DIFFEQ.AddImplementation(8, 8, 9, 25.11, RU_size);
		kernels.push_back(DIFFEQ);

		Kernel DCT(RandU(1, 10));
		DCT.AddImplementation(14, 14, 26, 54.08, RU_size);
		DCT.AddImplementation(15, 15, 13, 32.5, RU_size);
		DCT.AddImplementation(20, 20, 7, 18.83, RU_size);
		kernels.push_back(DCT);

		Kernel EWF(RandU(1, 10));
		EWF.AddImplementation(11, 11, 25, 52.5, RU_size);
		EWF.AddImplementation(13, 13, 13, 33.41, RU_size);
		EWF.AddImplementation(17, 17, 7, 17.71, RU_size);
		kernels.push_back(EWF);

		Kernel FIR1(RandU(1, 10));
		FIR1.AddImplementation(12, 12, 25, 52.25, RU_size);
		FIR1.AddImplementation(14, 14, 13, 31.85, RU_size);
		FIR1.AddImplementation(19, 19, 7, 17.22, RU_size);
		kernels.push_back(FIR1);

		Kernel FIR2(RandU(1, 10));
		FIR2.AddImplementation(10, 10, 25, 46.25, RU_size);
		FIR2.AddImplementation(12, 12, 13, 29.9, RU_size);
		FIR2.AddImplementation(17, 17, 7, 22.4, RU_size);
		kernels.push_back(FIR2);

		Kernel HAL(RandU(1, 10));
		HAL.AddImplementation(9, 9, 25, 30, RU_size);
		HAL.AddImplementation(10, 10, 13, 20, RU_size);
		HAL.AddImplementation(14, 14, 7, 17.01, RU_size);
		kernels.push_back(HAL);

		Kernel PAULIN(RandU(1, 10));
		PAULIN.AddImplementation(6, 6, 30, 32.1, RU_size);
		PAULIN.AddImplementation(7, 7, 16, 24.16, RU_size);
		PAULIN.AddImplementation(8, 8, 8, 18.64, RU_size);
		kernels.push_back(PAULIN);

		Kernel Wavelet(RandU(1, 10));
		Wavelet.AddImplementation(13, 13, 26, 52.26, RU_size);
		Wavelet.AddImplementation(14, 14, 13, 32.63, RU_size);
		Wavelet.AddImplementation(19, 19, 7, 22.4, RU_size);
		kernels.push_back(Wavelet);

		Kernel SmoothTriangle(RandU(1, 10));
		SmoothTriangle.AddImplementation(11, 11, 15, 30.9, RU_size);
		SmoothTriangle.AddImplementation(14, 14, 8, 18.16, RU_size);
		SmoothTriangle.AddImplementation(18, 18, 5, 17, RU_size);
		kernels.push_back(SmoothTriangle);

		float scale = 1.5, variation = 0.2;

		kernels.push_back(Kernel(ARF, scale, variation));
		kernels.push_back(Kernel(DIFFEQ, scale, variation));
		kernels.push_back(Kernel(EWF, scale, variation));
		kernels.push_back(Kernel(FIR1, scale, variation));
		kernels.push_back(Kernel(FIR2, scale, variation));
		kernels.push_back(Kernel(DCT, scale, variation));
		kernels.push_back(Kernel(HAL, scale, variation));
		kernels.push_back(Kernel(PAULIN, scale, variation));
		kernels.push_back(Kernel(Wavelet, scale, variation));
		kernels.push_back(Kernel(SmoothTriangle, scale, variation));

		/*kernels.push_back(Kernel(ARF, scale, variation));
		 kernels.push_back(Kernel(DIFFEQ));
		 kernels.push_back(Kernel(EWF));
		 kernels.push_back(Kernel(FIR1));
		 kernels.push_back(Kernel(FIR2));
		 kernels.push_back(Kernel(DCT));
		 kernels.push_back(Kernel(HAL));
		 kernels.push_back(Kernel(PAULIN));
		 kernels.push_back(Kernel(Wavelet));
		 kernels.push_back(Kernel(SmoothTriangle));*/

		/*Kernel random1 = new Kernel(RandU(1, 10));
		 random1.AddImplementation(4, 4, 23, 11.06, RU_size);
		 random1.AddImplementation(5, 5, 11, 10.09, RU_size);
		 random1.AddImplementation(7, 7, 6, 9, RU_size);
		 kernels.Add(random1);
		 kernels.Add(new Kernel(random1));
		 Kernel random2 = new Kernel(RandU(1, 10));
		 random2.AddImplementation(14, 14, 10, 11.06, RU_size);
		 random2.AddImplementation(15, 15, 5, 10.09, RU_size);
		 random2.AddImplementation(17, 17, 3, 9, RU_size);
		 kernels.Add(random2);
		 kernels.Add(new Kernel(random2));*/

		for (unsigned int i = 0; i < kernels.size(); ++i) {
			kernels[i].SetIndex(i);
		}
		kernels.erase(kernels.begin() + numKernels, kernels.end());
		OTreeEncoding::SetKernels(kernels);
	}

	void ProceedRandomGreedy(ObjectiveType type) {
		cout << "Random Greedy Algorithm Begins..." << endl;
		clock_t start;
		start = clock();
		unsigned i, j;
		unsigned size = kernels.size();
		double t;
		vector<int> order = vector<int>();
		vector<int> impls = vector<int>();

		vector<int> x = vector<int>(size);
		vector<int> y = vector<int>(size);

		for (i = 0; i < size; ++i) {
			order.push_back(i);
			impls.push_back(0);
		}
		i = 0;

		InsertOrder child;
		double minCost = numeric_limits<double>::max();

		FILE *f = fopen("montecarlo.txt", "w");

		do {
			for (j = 0; j < kernels.size(); ++j) {
				impls[j] = (RandU(0, kernels[j].GetNumImplements()));
			}
			child = InsertOrder(order, impls);
			child.GetGreedyTree(x, y, type);
			fprintf(f, "%d,%d\n", child.GetBestWidth(), child.GetBestLength());
			if (child.GetFitness() < minCost) {
				minCost = child.GetFitness();
				t = (clock() - start) / (double) CLOCKS_PER_SEC;
				printf("Better Individual Found at %d: %6.2f\t%dx%d Time: %6.2f\n", i, minCost, child.GetBestWidth(),
						child.GetBestLength(), t);
			}
			if (i % 100 == 0) {
				t = (clock() - start) / (double) CLOCKS_PER_SEC;
				printf("Individual #%d: Time: %6.2f\n", i, t);
			}
			if (t > 3600) {
				return;
			}
			i++;
		} while (true);
		//while (next_permutation(order.begin(), order.end()));

		fclose(f);

	}

	void ProceedHybridGA(unsigned PopSize, int MaxGen, int mutImpProb, int mutSeqProb, int crossProb,
			ObjectiveType type) {
		FILE* log = fopen("hybridGA.log", "a");
		cout << "HybridGA Initial Population Generation Begins..." << endl;
		clock_t start = clock();
		int i, j;
		int size = kernels.size();
		vector<int> order = vector<int>();
		vector<int> impls = vector<int>();

		for (i = 0; i < size; ++i) {
			order.push_back(i);
			impls.push_back(0);
		}

		vector<int> x = vector<int>(size);
		vector<int> y = vector<int>(size);

		double t;
		int best, worst;
		i = 0;
		vector<InsertOrder> hybridPop = vector<InsertOrder>();
		InsertOrder child;
		while (hybridPop.size() < PopSize) {
			random_shuffle(order.begin(), order.end());
			for (j = 0; j < size; ++j) {
				do {
					/*if (type == MaxUtility)
					 impls[j] = RandU(0, kernels[j].GetNumImplements() + 1) - 1;
					 else*/
					impls[j] = RandU(0, kernels[j].GetNumImplements());
					kernels[j].SetCurrentImplementation(impls[j]);
				} while (kernels[j].GetDepth() > OTreeEncoding::MAX_PHY_DEPTH);
			}
			child = InsertOrder(order, impls);
			//child.GetGreedyTree(x, y, type);
			if (find(hybridPop.begin(), hybridPop.end(), child) != hybridPop.end())
				continue;
			else {
				child.GetGreedyTree(x, y, type);
				hybridPop.push_back(child);
			}
			if (hybridPop.size() % 10 == 0)
				cout << hybridPop.size() << endl;
		}
		t = (clock() - start) / (double) CLOCKS_PER_SEC;
		cout << "HybridGA Initial Population Generation Time: " << t << "s" << endl;

		int generation = 0, p1, p2;

		double cost = numeric_limits<double>::max();
		double bestCost = numeric_limits<double>::max();

		cout << "Hybrid Genetic Algorithm Starts..." << endl;
		start = clock();
		while (generation < MaxGen) {
			generation++;
			while (hybridPop.size() < 3 * PopSize) {
				p1 = RandU(0, PopSize);
				if (RandU(0, 100) > crossProb) {
					child = hybridPop[p1].Mutation(mutImpProb, mutSeqProb);
				} else {
					p2 = RandU(0, PopSize);
					child = hybridPop[p1].CrossOver(hybridPop[p2]);
				}
				if (find(hybridPop.begin(), hybridPop.end(), child) != hybridPop.end())
					continue;

				child.GetGreedyTree(x, y, type);
				cost = child.GetFitness();
				if (cost < numeric_limits<double>::max()) {
					hybridPop.push_back(InsertOrder(child));
				}
			}

			sort(hybridPop.begin(), hybridPop.end());
			hybridPop.erase(hybridPop.begin() + PopSize, hybridPop.end());
			t = (clock() - start) / (double) CLOCKS_PER_SEC;
			best = 0;
			worst = hybridPop.size() - 1;
			OTreeEncoding encode;
			if (hybridPop[best].GetFitness() < bestCost) {
				printf("Generation %d\t Best Individual:  %6.2f\t%dx%d\tTime: %8f s\n", generation,
						hybridPop[best].GetFitness(), hybridPop[best].GetBestWidth(), hybridPop[best].GetBestLength(),
						t);
				bestCost = hybridPop[best].GetFitness();
				encode = hybridPop[best].GetGreedyTree(x, y, type);
				printf("Encoding: %s\n", encode.ToString().c_str());
			} //else if (generation % 10 == 0) {

			fprintf(log, "%d, %6.2f, %6.2f, %d, %d, %8f\n", generation, hybridPop[best].GetFitness(),
					hybridPop[worst].GetFitness(), hybridPop[best].GetBestWidth(), hybridPop[best].GetBestLength(), t);

			/*fprintf(log, "Generation %d\t Best Individual:  %6.2f\t%dx%d\tTime: %8f s\n", generation,
			 hybridPop[best].GetFitness(), hybridPop[best].GetBestWidth(), hybridPop[best].GetBestLength(), t);
			 fprintf(log, "Generation %d\t Worst Individual: %6.2f\t%dx%d\tTime: %8f s\n", generation,
			 hybridPop[worst].GetFitness(), hybridPop[worst].GetBestWidth(), hybridPop[worst].GetBestLength(),
			 t);*/
			/*if (hybridPop[best].GetFitness() == hybridPop[worst].GetFitness()
			 && hybridPop[best].GetBestWidth() == hybridPop[worst].GetBestWidth()
			 && hybridPop[best].GetBestLength() == hybridPop[worst].GetBestLength()) {
			 InsertOrder o1 = hybridPop[best];
			 InsertOrder o2 = hybridPop[worst];
			 if (o1 == o2) {
			 printf("worst equals best, something is wrong!\n");
			 exit(-1);
			 } else {
			 printf("worst is not best, something is magical!\n");
			 }
			 }*/

		}

		child = hybridPop[best];
		OTreeEncoding encode = child.GetGreedyTree(x, y, type);
		encode.Encoding2Mapping(x, y, type);
		int implement, index;
		FILE *file = fopen("placement.hybrid.csv", "w");
		for (i = 0; i < size; ++i) {
			index = child.GetOrder()[i];
			implement = child.GetImplementation()[index];
			kernels[index].SetCurrentImplementation(implement);
			fprintf(file, "%d %d %d %d %d %d\n", index, x[index], y[index], kernels[index].GetWidth(),
					kernels[index].GetHeight(), kernels[index].GetDepth());
		}

		fclose(file);
		fclose(log);
	}

	void ProceedHybridGAOrig(unsigned PopSize, int MaxGen, int mutImpProb, int mutSeqProb, int crossProb,
			ObjectiveType type) {
		cout << "HybridGA Initial Population Generation Begins..." << endl;
		FILE* log = fopen("hybridGA-2D.log", "a");
		clock_t start = clock();
		int i, j;
		int size = kernels.size();
		vector<int> order = vector<int>();
		vector<int> impls = vector<int>();

		for (i = 0; i < size; ++i) {
			order.push_back(i);
			impls.push_back(0);
		}

		vector<int> x = vector<int>(size);
		vector<int> y = vector<int>(size);

		double t;
		int best, worst;
		i = 0;
		vector<InsertOrder> hybridPop = vector<InsertOrder>();
		InsertOrder child;
		while (hybridPop.size() < PopSize) {
			random_shuffle(order.begin(), order.end());
			for (j = 0; j < size; ++j) {
				do {
					/*if (type == MaxUtility)
					 impls[j] = RandU(0, kernels[j].GetNumImplements() + 1)
					 - 1;
					 else
					 */impls[j] = RandU(0, kernels[j].GetNumImplements());
					kernels[j].SetCurrentImplementation(impls[j]);
				} while (kernels[j].GetDepth() > OTreeEncoding::MAX_DEPTH);
			}
			child = InsertOrder(order, impls);
			if (find(hybridPop.begin(), hybridPop.end(), child) != hybridPop.end())
				continue;
			else {
				child.GetGreedyTree2D(x, y, type);
				hybridPop.push_back(child);
			}
			if (hybridPop.size() % 10 == 0)
				cout << hybridPop.size() << endl;
		}
		t = (clock() - start) / (double) CLOCKS_PER_SEC;
		cout << "HybridGA Initial Population Generation Time: " << t << "s" << endl;

		int generation = 0, p1, p2;

		double cost = numeric_limits<double>::max();
		double bestCost = numeric_limits<double>::max();

		cout << "Hybrid Genetic Algorithm Starts..." << endl;
		start = clock();
		while (generation < MaxGen) {
			generation++;
			while (hybridPop.size() < 2 * PopSize) {
				p1 = RandU(0, PopSize);
				if (RandU(0, 100) > crossProb) {
					child = hybridPop[p1].Mutation(mutImpProb, mutSeqProb);
				} else {
					p2 = RandU(0, PopSize);
					child = hybridPop[p1].CrossOver(hybridPop[p2]);
				}
				if (find(hybridPop.begin(), hybridPop.end(), child) != hybridPop.end())
					continue;

				child.GetGreedyTree2D(x, y, type);
				cost = child.GetFitness();
				if (cost < numeric_limits<double>::max()) {
					hybridPop.push_back(child);
				}
			}

			sort(hybridPop.begin(), hybridPop.end());
			hybridPop.erase(hybridPop.begin() + PopSize, hybridPop.end());
			t = (clock() - start) / (double) CLOCKS_PER_SEC;
			best = 0;
			worst = hybridPop.size() - 1;
			OTreeEncoding2D encode;
			if (hybridPop[best].GetFitness() < bestCost) {
				printf("Generation %d\t Best Individual:  %6.2f\t%dx%d\tTime: %8f s\n", generation,
						hybridPop[best].GetFitness(), hybridPop[best].GetBestWidth(), hybridPop[best].GetBestLength(),
						t);
				bestCost = hybridPop[0].GetFitness();
				encode = hybridPop[best].GetGreedyTree2D(x, y, type);
				printf("Encoding: %s\n", encode.ToString().c_str());
			}
			fprintf(log, "%d, %6.2f, %6.2f, %d, %d, %8f\n", generation, hybridPop[best].GetFitness(),
					hybridPop[worst].GetFitness(), hybridPop[best].GetBestWidth(), hybridPop[best].GetBestLength(), t);
			/*else if (generation % 10 == 0) {
			 printf("Generation %d\t Best Individual:  %6.2f\t%dx%d\tTime: %8f s\n", generation,
			 hybridPop[best].GetFitness(), hybridPop[best].GetBestWidth(), hybridPop[best].GetBestLength(),
			 t);
			 printf("Generation %d\t Worst Individual: %6.2f\t%dx%d\tTime: %8f s\n", generation,
			 hybridPop[worst].GetFitness(), hybridPop[worst].GetBestWidth(),
			 hybridPop[worst].GetBestLength(), t);
			 if (hybridPop[best].GetFitness() == hybridPop[worst].GetFitness()
			 && hybridPop[best].GetBestWidth() == hybridPop[worst].GetBestWidth()
			 && hybridPop[best].GetBestLength() == hybridPop[worst].GetBestLength()) {
			 InsertOrder o1 = hybridPop[best];
			 InsertOrder o2 = hybridPop[worst];
			 if (o1 == o2) {
			 exit(-1);
			 }
			 }
			 }*/
		}

		child = hybridPop[best];
		OTreeEncoding2D encode = child.GetGreedyTree2D(x, y, type);
		encode.Encoding2Mapping(x, y, type);
		int implement, index;
		FILE *file = fopen("placement.hybrid.csv", "w");
		for (i = 0; i < size; ++i) {
			index = child.GetOrder()[i];
			implement = child.GetImplementation()[index];
			kernels[index].SetCurrentImplementation(implement);
			fprintf(file, "%d %d %d %d %d %d\n", index, x[index], y[index], kernels[index].GetWidth(),
					kernels[index].GetHeight(), kernels[index].GetDepth());
		}

		fclose(file);
		fclose(log);
	}

	void ProceedGA(unsigned PopSize, int MaxGen, int mutImpProb, int mutSeqProb, int crossProb, ObjectiveType type) {
		printf("GA Initial Population Generation Begins...\n");
		FILE *log = fopen("GA.log", "a");
		clock_t start = clock();
		this->GenerateInitialPopulation(PopSize, type);
		double t = (clock() - start) / (double) CLOCKS_PER_SEC;
		printf("Initial Population Generation Time: %8fs", t);

		int generation = 0, p1, p2;

		double cost = numeric_limits<double>::max();
		double bestCost = numeric_limits<double>::max();
		int i;
		OTreeEncoding child;
		int size = kernels.size();
		InsertOrder o;

		vector<int> order = vector<int>();
		vector<int> impls = vector<int>();

		for (i = 0; i < size; ++i) {
			order.push_back(i);
			impls.push_back(0);
		}

		vector<int> x = vector<int>(size);
		vector<int> y = vector<int>(size);

		printf("General Genetic Algorithm Starts...\n");
		start = clock();

		while (generation < MaxGen) {
			generation++;

			while (population.size() < PopSize * 2) {
				p1 = RandU(0, PopSize);
				if (RandU(0, 100) < crossProb) {
					p2 = RandU(0, PopSize);
					child = this->CrossOver(population[p1], population[p2]);
				} else {
					child = population[p1].GetMutation(mutImpProb, mutSeqProb);
				}

				cost = child.Encoding2Mapping(x, y, type);
				if (find(population.begin(), population.end(), child) != population.end())
					continue;

				if (cost < numeric_limits<double>::max()) {
					child.SetFitness(cost);
					population.push_back(OTreeEncoding(child));
				}

			}

			sort(population.begin(), population.end());
			population.erase(population.begin() + PopSize, population.end());
			int worst = population.size() - 1;
			int best = 0;
			t = (clock() - start) / (double) CLOCKS_PER_SEC;
			if (population[best].GetFitness() < bestCost) {
				printf("Generation %d\t Best Individual:  %6.2f\t%dx%d\tTime: %6.2f s\n", generation,
						population[best].GetFitness(), population[best].GetWidth(), population[best].GetLength(), t);
				bestCost = population[0].GetFitness();
			}
			fprintf(log, "%d, %6.2f, %6.2f, %d, %d, %8f\n", generation, population[best].GetFitness(),
					population[worst].GetFitness(), population[best].GetWidth(), population[best].GetLength(), t);
			/*else if (generation % 10 == 0) {
			 printf("Generation %d\t Best Individual:  %6.2f\t%dx%d\tTime: %6.2f s\n", generation,
			 population[best].GetFitness(), population[best].GetWidth(), population[best].GetLength(), t);
			 printf("Generation %d\t Worst Individual: %6.2f\t%dx%d\tTime: %6.2f s\n", generation,
			 population[worst].GetFitness(), population[worst].GetWidth(), population[worst].GetLength(), t);
			 //				printf("%d\n", population[best] == population[worst]);
			 }*/
			random_shuffle(population.begin(), population.end());
		}
		fclose(log);
	}

	void ProceedGAOrig(unsigned PopSize, int MaxGen, int mutImpProb, int mutSeqProb, int crossProb,
			ObjectiveType type) {
		printf("GAOrig Initial Population Generation Begins...\n");
		FILE *log = fopen("GA-2D.log", "a");
		clock_t start = clock();
		this->GenerateInitialPopulation2D(PopSize, type);
		double t = (clock() - start) / (double) CLOCKS_PER_SEC;
		printf("Initial Population Generation Time: %8fs", t);

		int generation = 0, p1, p2;

		double cost = numeric_limits<double>::max();
		double bestCost = numeric_limits<double>::max();
		int i;
		OTreeEncoding2D child;
		int size = kernels.size();
		InsertOrder o;

		vector<int> order = vector<int>();
		vector<int> impls = vector<int>();

		for (i = 0; i < size; ++i) {
			order.push_back(i);
			impls.push_back(0);
		}

		vector<int> x = vector<int>(size);
		vector<int> y = vector<int>(size);


		printf("General Genetic Algorithm Starts...\n");
		start = clock();



		while (generation < MaxGen) {
			generation++;

			while (population2D.size() < PopSize * 2) {
				p1 = RandU(0, PopSize);
				if (RandU(0, 100) < crossProb) {
					p2 = RandU(0, PopSize);
					child = this->CrossOver2D(population2D[p1], population2D[p2]);
				} else {
					child = population2D[p1].GetMutation(mutImpProb, mutSeqProb);
				}

				cost = child.Encoding2Mapping(x, y, type);
				if (find(population2D.begin(), population2D.end(), child) != population2D.end())
					continue;

				if (cost < numeric_limits<double>::max()) {
					child.SetFitness(cost);
					population2D.push_back(OTreeEncoding2D(child));
				}

			}

			sort(population2D.begin(), population2D.end());
			population2D.erase(population2D.begin() + PopSize, population2D.end());
			int worst = population2D.size() - 1;
			int best = 0;
			t = (clock() - start) / (double) CLOCKS_PER_SEC;
			if (population2D[best].GetFitness() < bestCost) {
				printf("Generation %d\t Best Individual:  %6.2f\t%dx%d\tTime: %6.2f s\n", generation,
						population2D[best].GetFitness(), population2D[best].GetWidth(), population2D[best].GetLength(),
						t);
				bestCost = population2D[0].GetFitness();
			}
			fprintf(log, "%d, %6.2f, %6.2f, %d, %d, %8f\n", generation, population2D[best].GetFitness(),
					population2D[worst].GetFitness(), population2D[best].GetWidth(), population2D[best].GetLength(), t);
			/*else if (generation % 10 == 0) {
			 printf("Generation %d\t Best Individual:  %6.2f\t%dx%d\tTime: %6.2f s\n", generation,
			 population2D[best].GetFitness(), population2D[best].GetWidth(), population2D[best].GetLength(),
			 t);
			 printf("Generation %d\t Worst Individual: %6.2f\t%dx%d\tTime: %6.2f s\n", generation,
			 population2D[worst].GetFitness(), population2D[worst].GetWidth(),
			 population2D[worst].GetLength(), t);
			 //				printf("%d", population2D[best] == population2D[worst]);
			 }*/
			random_shuffle(population2D.begin(), population2D.end());
		}

		fclose(log);

	}

	OTreeEncoding Mutation(OTreeEncoding p, int mutImpProb, int mutSeqProb) {
		return p.GetMutation(mutImpProb, mutSeqProb);
	}

	OTreeEncoding CrossOver(OTreeEncoding p1, OTreeEncoding p2) {
		vector<int> part1 = vector<int>();
		vector<int> part2 = vector<int>();

		int size = kernels.size();

		int seperate = size / 2;

		vector<int> imp = vector<int>();
		for (int i = 0; i < size; ++i) {
			imp.push_back(-1);
		}

		for (int i = 0; i < size; ++i) {
			if (RandU(0, size) > seperate) {
				part1.push_back(i);
				imp[i] = p1.GetImplementAt(i);
			} else {
				part2.push_back(i);
				imp[i] = p2.GetImplementAt(i);
			}
		}

		vector<int> seq_p1, seq_p2;
		vector<bool> tra_p1, tra_p2;

		p1.GetPartialTree(part1, &seq_p1, &tra_p1);
		p2.GetPartialTree(part2, &seq_p2, &tra_p2);

		seq_p1.insert(seq_p1.end(), seq_p2.begin(), seq_p2.end());
		tra_p1.insert(tra_p1.end(), tra_p2.begin(), tra_p2.end());

		OTreeEncoding child = OTreeEncoding(tra_p1, seq_p1, imp);

		return child;
	}

	OTreeEncoding2D CrossOver2D(OTreeEncoding2D p1, OTreeEncoding2D p2) {
		vector<int> part1 = vector<int>();
		vector<int> part2 = vector<int>();

		int size = kernels.size();

		int seperate = size / 2;

		vector<int> imp = vector<int>();
		for (int i = 0; i < size; ++i) {
			imp.push_back(-1);
		}

		for (int i = 0; i < size; ++i) {
			if (RandU(0, size) > seperate) {
				part1.push_back(i);
				imp[i] = p1.GetImplementAt(i);
			} else {
				part2.push_back(i);
				imp[i] = p2.GetImplementAt(i);
			}
		}

		vector<int> seq_p1H, seq_p2H;
		vector<int> seq_p1V, seq_p2V;
		vector<bool> tra_p1H, tra_p2H;
		vector<bool> tra_p1V, tra_p2V;

		p1.GetPartialTree(part1, &seq_p1H, &tra_p1H, H);
		p2.GetPartialTree(part2, &seq_p2H, &tra_p2H, H);
		p1.GetPartialTree(part1, &seq_p1V, &tra_p1V, V);
		p2.GetPartialTree(part2, &seq_p2V, &tra_p2V, V);

		seq_p1H.insert(seq_p1H.end(), seq_p2H.begin(), seq_p2H.end());
		seq_p1V.insert(seq_p1V.end(), seq_p2V.begin(), seq_p2V.end());
		tra_p1H.insert(tra_p1H.end(), tra_p2H.begin(), tra_p2H.end());
		tra_p1V.insert(tra_p1V.end(), tra_p2V.begin(), tra_p2V.end());

		OTreeEncoding2D child = OTreeEncoding2D(tra_p1H, tra_p1V, seq_p1H, seq_p1V, imp);

		return child;
	}

	void GenerateInitialPopulation(unsigned popsize, ObjectiveType type) {
		int i, j;
		vector<Kernel> leftKernels = vector<Kernel>(kernels);
		double fit;

		int size = kernels.size();
		vector<int> seq = vector<int>();
		vector<int> impls = vector<int>();
		vector<bool> traversal = vector<bool>();

		for (i = 0; i < size; ++i) {
			seq.push_back(i);
			impls.push_back(0);
			traversal.push_back(true);
		}

		for (i = 0; i < size; ++i) {
			traversal.push_back(false);
		}

		vector<int> x = vector<int>(size);
		vector<int> y = vector<int>(size);

		OTreeEncoding child;

		int validation = 0;

		population.clear();

		while (population.size() < popsize) {
			random_shuffle(seq.begin(), seq.end());
			do {
				random_shuffle(traversal.begin() + validation, traversal.end());
			} while ((validation = validTraversal(traversal)) != traversal.size());
			for (j = 0; j < size; ++j) {
				impls[j] = (RandU(0, kernels[j].GetNumImplements()));
			}
			child = OTreeEncoding(traversal, seq, impls);
			if (find(population.begin(), population.end(), child) != population.end())
				continue;
			else {
				fit = child.Encoding2Mapping(x, y, type);
				child.SetFitness(fit);
				population.push_back(child);
			}
		}

	}

	void GenerateInitialPopulation2D(unsigned popsize, ObjectiveType type) {
		int i, j;
		vector<Kernel> leftKernels = vector<Kernel>(kernels);
		double fit;

		int size = kernels.size();
		vector<int> seqH = vector<int>();
		vector<int> seqV = vector<int>();
		vector<int> impls = vector<int>();
		vector<bool> traversalH = vector<bool>();
		vector<bool> traversalV = vector<bool>();

		for (i = 0; i < size; ++i) {
			seqH.push_back(i);
			seqV.push_back(i);
			impls.push_back(0);
			traversalH.push_back(true);
			traversalV.push_back(true);
		}

		for (i = 0; i < size; ++i) {
			traversalH.push_back(false);
			traversalV.push_back(false);
		}

		vector<int> x = vector<int>(size);
		vector<int> y = vector<int>(size);

		OTreeEncoding2D child;

		int validation = 0;

		population2D.clear();

		while (population2D.size() < popsize) {
			random_shuffle(seqH.begin(), seqH.end());
			random_shuffle(seqV.begin(), seqV.end());
			do {
				random_shuffle(traversalH.begin() + validation, traversalH.end());
			} while ((validation = validTraversal(traversalH)) != traversalH.size());
			validation = 0;
			do {
				random_shuffle(traversalV.begin() + validation, traversalV.end());
			} while ((validation = validTraversal(traversalV)) != traversalV.size());
			for (j = 0; j < size; ++j) {
				impls[j] = (RandU(0, kernels[j].GetNumImplements()));
			}
			child = OTreeEncoding2D(traversalH, traversalV, seqH, seqV, impls);
			if (find(population2D.begin(), population2D.end(), child) != population2D.end())
				continue;
			else {
				fit = child.Encoding2Mapping(x, y, type);
				child.SetFitness(fit);
				population2D.push_back(child);
			}
		}

	}

	int validTraversal(vector<bool> traversal) {
		int valid = 0, i;
		for (i = 0; i < traversal.size(); ++i) {
			if (traversal[i])
				valid--;
			else
				valid++;
			if (valid < 0)
				return i;
		}
		return traversal.size();
	}

};

int main(int argc, char** argv) {
	int numKernels = 10;
	int PopSize = 5 * numKernels;
	int MaxGen = 500;
	int mutImpProb = 90;
	int mutSeqProb = 90;
	int crossProb = 10;

	OTreeEncoding::MAX_DEPTH = 32;
	OTreeEncoding::MAX_PHY_DEPTH = 32;
	OTreeEncoding::MAX_WIDTH = 20;
	OTreeEncoding::MAX_HEIGHT = 20;
	int numTasks = 1000;
	float arrivalRate = 0.5;

	srand(time(NULL));

	std::cerr << std::fixed;
	std::cerr.precision(5);

	ObjectiveType type = MinArea;

	GeneticAlgorithm ga = GeneticAlgorithm(numKernels, 1);
	for (int i = 0; i < 6; ++i) {
		ga.ProceedHybridGA(PopSize, MaxGen, mutImpProb, mutSeqProb, crossProb, type);
		ga.ProceedHybridGAOrig(PopSize, MaxGen, mutImpProb, mutSeqProb, crossProb, type);

		ga.ProceedGAOrig(PopSize, MaxGen, mutImpProb, mutSeqProb, crossProb, type);
		ga.ProceedGA(PopSize, MaxGen, mutImpProb, mutSeqProb, crossProb, type);
	}
	//ga.ProceedRandomGreedy(type);

	/*Scheduler s = Scheduler(numKernels, numTasks);

	 s.TaskGenerator(OTreeEncoding::MAX_DEPTH, arrivalRate);
	 s.Schedule(EDF, IMMEDIATE, 0, 0);
	 s.Schedule(EDF, IMMEDIATE, 1, 0);
	 s.Schedule(EDF, SLACK, 0, 0);
	 s.Schedule(EDF, SLACK, 1, 0);
	 s.TaskGenerator(OTreeEncoding::MAX_DEPTH, arrivalRate);
	 s.Schedule(PATS, IMMEDIATE, 0, 0);
	 s.Schedule(PATS, IMMEDIATE, 1, 0);
	 s.Schedule(PATS, SLACK, 0, 0);
	 s.Schedule(PATS, SLACK, 1, 0);*/
	/*for (int i = 0; i < 11; ++i) {
	 s.TaskGenerator(OTreeEncoding::MAX_DEPTH, 0.75 + 0.05 * i, arrivalRate);
	 s.Schedule(PATS, 0);
	 }*/
	/*s.TaskGenerator(OTreeEncoding::MAX_DEPTH, 1.2, arrivalRate);
	 s.Schedule(EDF, 0);
	 s.TaskGenerator(OTreeEncoding::MAX_DEPTH, 0.9, arrivalRate);
	 s.Schedule(EDF, 0);*/

	/*arr_list[0] = Task(0, 0, 0, 18);
	 arr_list[1] = Task(1, 2, 0, 14);
	 arr_list[2] = Task(2, 6, 3, 29);
	 arr_list[3] = Task(3, 8, 4, 30);
	 arr_list[4] = Task(4, 7, 9, 25);*/

	return 0;
}

