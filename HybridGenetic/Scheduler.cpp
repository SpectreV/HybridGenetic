/*
 * Scheduler.cpp
 *
 *  Created on: 5 Jan 2016
 *      Author: liang
 */

#include "Scheduler.h"
#include "OTreeEncoding.h"
#include <algorithm>
#include <random>
#include <iostream>
#include <cassert>

namespace std {

Scheduler::Scheduler() {
	this->numKernels = 0;
	this->numTasks = 0;
	this->maxW = this->maxH = this->maxD = 0;
}

Scheduler::Scheduler(const int numKernels, const int numTasks) {
	arr_list = vector<Task>();
	wait_list = vector<Task>();
	exe_list = vector<Task>();
	fin_list = vector<Task>();
	this->numKernels = numKernels;
	this->numTasks = numTasks;
	this->maxW = this->maxH = this->maxD = 0;
}

Scheduler::~Scheduler() {
	// TODO Auto-generated destructor stub
}

void Scheduler::TaskGenerator(const int MaxDepth, const float arrivalRate) {
	int i, j, k, size = numKernels;
	FILE *file = fopen("placement.hybrid.csv", "r");
	place = vector<vector<int> >(size, vector<int>(6, 0));
	delay = vector<double>(numKernels);
	maxD = OTreeEncoding::MAX_PHY_DEPTH;
	maxH = OTreeEncoding::MAX_HEIGHT;
	maxW = OTreeEncoding::MAX_WIDTH;

	for (i = 0; i < size; ++i) {
		if (fscanf(file, "%d %d %d %d %d %d", &place[i][0], &place[i][1], &place[i][2], &place[i][3], &place[i][4],
				&place[i][5]) != 6) {
			cerr << "Error reading placement result!" << endl;
			exit(-1);
		}
		delay[place[i][0]] = place[i][5];
	}
	fclose(file);

	offchip = vector<int>(size, 0);

	vector<int> temp;
	for (i = 0; i < size; ++i) {
		for (j = i; j < size; ++j) {
			if (place[j][0] == i) {
				temp = place[i];
				place[i] = place[j];
				place[j] = temp;
			}
		}
	}

	depth = vector<vector<int> >(maxH, vector<int>(maxW));
	for (i = 0; i < size; ++i) {
		for (j = place[i][2]; j < place[i][2] + place[i][4]; ++j) {
			for (k = place[i][1]; k < place[i][1] + place[i][3]; ++k) {
				depth[j][k] += place[i][5];
				if (depth[j][k] > maxD)
					offchip[i] = 1;
			}
		}
		if (offchip[i]) {
			for (j = place[i][2]; j < place[i][2] + place[i][4]; ++j) {
				for (k = place[i][1]; k < place[i][1] + place[i][3]; ++k) {
					depth[j][k] -= place[i][5];
				}
			}
		}
	}

	/*cout << "Initial Depth: " << endl;
	 PrintDepth();*/

	conflict = vector<vector<int> >(size, vector<int>(size, 0));

	for (i = 0; i < size; ++i) {
		if (offchip[i])
			continue;
		for (j = i + 1; j < size; ++j) {
			if (Overlap(place[i][1], place[i][2], place[i][3], place[i][4], place[j][1], place[j][2], place[j][3],
					place[j][4]) && !offchip[j]) {
				conflict[i][j] = 1;
				conflict[j][i] = 1;
			}
		}
	}

	/*for (i = 0; i < size; ++i) {
	 for (j = 0; j < size; ++j) {
	 if (conflict[i][j]) {
	 cout << "Kernel " << place[i][0] << " overlaps kernel " << place[j][0] << endl;
	 }
	 }
	 }*/

	arr_list = vector<Task>();

	wcet = vector<int>(numKernels, 0);
	bcet = vector<int>(numKernels, 0);

	for (i = 0; i < numKernels; ++i) {
		for (j = 0; j < OTreeEncoding::kernels[i].GetNumImplements(); ++j) {
			if (wcet[i] < OTreeEncoding::kernels[i].GetDepth(j))
				wcet[i] = OTreeEncoding::kernels[i].GetDepth(j);
			if (bcet[i] > OTreeEncoding::kernels[i].GetDepth(j))
				bcet[i] = OTreeEncoding::kernels[i].GetDepth(j);
		}
	}

	int t = 0;

	/*i = 0;
	 while (i < numTasks) {
	 for (j = 0; j < numKernels; ++j) {
	 if (t % (3 * (int) wcet[j]) == 0) {
	 arr_list.push_back(Task(i, j, t, t + wcet[j] * 3));
	 i++;
	 }
	 }
	 t++;
	 }*/

	default_random_engine generator(time(NULL));
	exponential_distribution<double> interval(arrivalRate);
	int index, arrT, ddl;
	for (int i = 0; i < numTasks; ++i) {
		index = RandU(0, numKernels);
		arrT = t;
		ddl = t + delay[index] * 3;
		arr_list.push_back(Task(i, index, arrT, ddl));
		t = t + interval(generator);
	}
}

void Scheduler::Schedule(SchedulePolicy sp, PlacementPolicy pp, bool replace, int debug) {
	kernel_in_use = vector<int>(numKernels, 0);
	average_exe_time = delay;
	exe_count = vector<int>(numKernels, 0);

	wait_list.erase(wait_list.begin(), wait_list.end());
	exe_list.erase(exe_list.begin(), exe_list.end());
	fin_list.erase(fin_list.begin(), fin_list.end());

	unsigned i, j, k, tempx, tempy;
	int index, t, evictions = 0;

	int schedule = 0, temp;
	double priority = 0, alpha = 2, max_wait_time;
	bool blocked;

	int rejected = 0;
	bool placable;

	t = 0;

	while (fin_list.size() != numTasks) {
		for (i = 0; i < arr_list.size(); ++i) {
			if (arr_list[i].getArrTime() == t) {
				if (debug)
					cout << "Task " << arr_list[i].getId() << " arrived requiring kernel " << arr_list[i].getIndex()
							<< " at cycle " << t << endl;
				wait_list.push_back(Task(arr_list[i]));
				schedule = 1;
			}
			if (arr_list[i].getArrTime() > t)
				break;
		}
		for (i = 0; i < exe_list.size(); ++i) {
			if (exe_list[i].getFinTime() == t) {
				if (debug)
					cout << "Task " << exe_list[i].getId() << " released kernel " << exe_list[i].getIndex()
							<< " at cycle " << t << endl;
				fin_list.push_back(Task(exe_list[i]));
				kernel_in_use[exe_list[i].getIndex()] = 0;
				exe_list.erase(exe_list.begin() + i);
				i--;
				schedule = 1;

			}
		}
		if (schedule) {
			if (sp == EDF) {
				for (i = 0; i < wait_list.size(); ++i) {
					wait_list[i].setPriority(wait_list[i].getDeadline() - t);
				}
			} else if (sp == PATS) {
				for (i = 0; i < wait_list.size(); ++i) {
					index = wait_list[i].getIndex();
					if (!kernel_in_use[i])
						priority = alpha * 1;
					else {
						max_wait_time = 0;
						for (k = 0; k < exe_list.size(); ++k) {
							temp = exe_list[k].getIndex();
							if (conflict[index][temp] && temp != index
									&& max_wait_time < exe_list[k].getFinTime() - t) {
								max_wait_time = exe_list[k].getFinTime() - t;
							}
						}
						priority = alpha * delay[i] / (delay[i] + max_wait_time);
					}

					priority -= (wait_list[i].getDeadline() - t) / average_exe_time[index];
					wait_list[i].setPriority(priority);
				}
			}

			sort(wait_list.begin(), wait_list.end());
			k = 0;

			while (k < wait_list.size()) {
				index = wait_list[k].getIndex();

				if (wait_list[k].getDeadline() - t < bcet[index]) {
					if (debug)
						cout << "Task " << index << " rejected exceeding its slack time." << endl;
					rejected++;
					wait_list.erase(wait_list.begin() + k);
					fin_list.push_back(Task(wait_list[k]));
					k--;
					goto next_task;
				}

				// Exact kernel occupied by another task
				if (kernel_in_use[index] && pp == IMMEDIATE) {
					rejected++;
					if (debug)
						cout << "Task " << wait_list[k].getId() << " rejected as kernel " << index << " is in use."
								<< endl;
					fin_list.push_back(Task(wait_list[k]));
					wait_list.erase(wait_list.begin() + k);
					k--;
					goto next_task;
				} else if (kernel_in_use[index] && pp == SLACK) {
					goto next_task;
				} else {
					// Kernel overlap with other active kernels
					blocked = false;
					for (i = 0; i < numKernels; ++i) {
						if (conflict[index][i] && kernel_in_use[i] && i != index) {
							blocked = true;
							break;
						}
					}

					// Kernel available on chip
					if (!offchip[index] && !blocked) {
						wait_list[k].setStartTime(t);
						wait_list[k].setFinTime(t + 2 * delay[index]);
						kernel_in_use[index] = 1;
					}
					// Kernel available offchip
					else {
						placable = Place(index, replace);
						if (!placable && pp == SLACK) {
							goto next_task;
						} else if (!placable && pp == IMMEDIATE) {
							rejected++;
							if (debug)
								cout << "Task " << wait_list[k].getId() << " rejected as kernel " << index
										<< " is blocked by other kernels." << endl;
							fin_list.push_back(Task(wait_list[k]));
							wait_list.erase(wait_list.begin() + k);
							k--;
							goto next_task;
						}

						wait_list[k].setStartTime(t + 1);
						wait_list[k].setFinTime(t + 1 + delay[index]);
						kernel_in_use[index] = 1;
						offchip[index] = 0;

						make_room:

						for (i = place[index][1]; i < place[index][1] + place[index][3]; ++i) {
							for (j = place[index][2]; j < place[index][2] + place[index][4]; ++j) {
								if (depth[j][i] > maxD) {
									for (temp = 0; temp < numKernels; ++temp) {
										if (conflict[index][temp] && !offchip[temp] && temp != index
												&& i >= place[temp][1] && i < place[temp][1] + place[temp][3]
												&& j >= place[temp][2] && j < place[temp][2] + place[temp][4]) {
											for (tempx = place[temp][1]; tempx < place[temp][1] + place[temp][3];
													++tempx) {
												for (tempy = place[temp][2]; tempy < place[temp][2] + place[temp][4];
														++tempy) {
													depth[tempy][tempx] -= delay[temp];
												}
											}
											offchip[temp] = 1;
											evictions++;
											/*cout << "Kernel " << temp << " evicted from " << place[temp][1] << "\t"
											 << place[temp][2] << "\t";
											 cout << place[temp][3] << "\t" << place[temp][4] << "\t" << place[temp][5]
											 << endl;
											 PrintDepth();*/
											for (tempx = 0; tempx < numKernels; ++tempx) {
												if (temp != tempx
														&& Overlap(place[tempx][1], place[tempx][2], place[tempx][3],
																place[tempx][4], place[temp][1], place[temp][2],
																place[temp][3], place[temp][4]))
													conflict[temp][tempx] = conflict[tempx][temp] = 0;
											}
											if (debug)
												cout << "Kernel " << temp << " evicted by kernel " << index
														<< " at cycle " << t << endl;
											goto make_room;
										}
									}
								}
							}
						}
					}
					if (debug)
						cout << "Task " << wait_list[k].getId() << " occupied kernel " << wait_list[k].getIndex()
								<< " at cycle " << wait_list[k].getStartTime() << endl;
					average_exe_time[index] = average_exe_time[index] * exe_count[index] + wait_list[k].getFinTime()
							- wait_list[k].getArrTime();
					exe_count[index]++;
					average_exe_time[index] /= exe_count[index];
					exe_list.push_back(Task(wait_list[k]));
					wait_list.erase(wait_list.begin() + k);
					k--;
				}
				next_task: k++;
			}
			schedule = 0;
		}
		t++;
	}

	string policy;

	switch (sp) {
	case EDF:
		policy = "EDF";
		break;
	case PATS:
		policy = "PATS";
		break;
	case CUSTOM:
		policy = "CUSTOM";
		break;
	default:
		policy = "UNKNOWN";
		break;
	}

	FILE *placer = fopen("online_placer.csv", "a");

	cerr << "Total execution cycles of " << policy << ": " << t << endl;
	/*cerr << "Kernel execution count: ";
	 for (i = 0; i < numKernels; ++i) {
	 cerr << exe_count[i] << "\t";
	 }
	 cerr << "\nKernel average response time: ";
	 for (i = 0; i < numKernels; ++i) {
	 cerr << average_exe_time[i] << "\t";
	 }*/
	temp = 0;
	for (i = 0; i < fin_list.size(); ++i) {
		if (fin_list[i].getIndex() == 2) {
			if (debug)
				cout << "Task " << i << " arrived at " << fin_list[i].getArrTime() << " finished at "
						<< fin_list[i].getFinTime() << endl;
		}
		if (fin_list[i].getFinTime() > fin_list[i].getDeadline()) {
			temp += fin_list[i].getFinTime() - fin_list[i].getDeadline();
		}
	}
	cerr << "System Tardiness: " << temp << endl;
	cerr << "Rejection rate: " << (double) rejected / numTasks << endl;
	cerr << "Evictions happened: " << evictions << endl << endl;

	/*cerr << "Kernel delay time: ";
	 for (i = 0; i < numKernels; ++i) {
	 cerr << delay[i] << "\t";
	 }*/

	policy = "replacement." + policy + ".csv";

	fprintf(placer, "%d %d %d\n", t, temp, evictions);
	fclose(placer);

	FILE *file = fopen(policy.c_str(), "w");
	for (i = 0; i < numKernels; ++i) {
		fprintf(file, "%d %d %d %d %d %d\n", place[i][0], place[i][1], place[i][2], place[i][3], place[i][4],
				place[i][5]);
	}
	fclose(file);
}

bool Scheduler::Place(int index, bool replace) {
	int new_x, new_y, new_imp;
	int i, j, k;
	double best_place = numeric_limits<double>::max(), current_place = 0;

	if (replace) {
		vector<int> possible_x, possible_y;
		possible_x.push_back(0);
		possible_y.push_back(0);
		for (i = 0; i < numKernels; ++i) {
			if (i == index)
				continue;
			//if (place[i][1] + place[index][3] < maxW)
			possible_x.push_back(place[i][1]);
			possible_x.push_back(place[i][1] + place[i][3]);
			/*if (place[i][1] + place[i][3] - place[index][3] >= 0)
			 possible_x.push_back(place[i][1] + place[i][3] - place[index][3]);*/
			//if (place[i][2] + place[index][4] < maxH)
			possible_y.push_back(place[i][2]);
			possible_y.push_back(place[i][2] + place[i][4]);

			/*if (place[i][2] + place[i][4] - place[index][4] >= 0)
			 possible_y.push_back(place[i][2] + place[i][4] - place[index][4]);*/

		}
		/*possible_x.push_back(maxW - place[index][3]);
		 possible_x.push_back(maxH - place[index][4]);*/

		sort(possible_x.begin(), possible_x.end());
		possible_x.erase(unique(possible_x.begin(), possible_x.end()), possible_x.end());
		sort(possible_y.begin(), possible_y.end());
		possible_y.erase(unique(possible_y.begin(), possible_y.end()), possible_y.end());

		int imp, w, h, d;
		for (i = 0; i < possible_x.size(); ++i) {
			for (j = 0; j < possible_y.size(); ++j) {
				for (imp = 0; imp < OTreeEncoding::kernels[index].GetNumImplements(); ++imp) {
					OTreeEncoding::kernels[index].SetCurrentImplementation(imp);
					w = OTreeEncoding::kernels[index].GetWidth();
					h = OTreeEncoding::kernels[index].GetHeight();
					d = OTreeEncoding::kernels[index].GetDepth();
					if (possible_x[i] < 0 || possible_x[i] + w >= maxW || possible_y[j] < 0 || possible_y[j] + h >= maxH
							|| d > maxD)
						goto next_choice;
					for (k = 0; k < numKernels; ++k) {
						if (k == index)
							continue;
						if (Overlap(possible_x[i], possible_y[j], w, h, place[k][1], place[k][2], place[k][3],
								place[k][4])) {
							if (kernel_in_use[k])
								goto next_choice;
							else if (offchip[k])
								current_place += delay[k];
							else
								current_place += d * delay[k];
						}
					}
					if (current_place < best_place) {
						best_place = current_place;
						new_x = possible_x[i];
						new_y = possible_y[j];
						new_imp = imp;
					}
					next_choice: current_place = 0;
				}
			}
		}
	} else {
		best_place = 0;
		for (k = 0; k < numKernels; ++k) {
			if (k == index)
				continue;
			if (Overlap(place[index][1], place[index][2], place[index][3], place[index][4], place[k][1], place[k][2],
					place[k][3], place[k][4])) {
				if (kernel_in_use[k]) {
					best_place = numeric_limits<double>::max();
					break;
				}
			}
		}

	}

	if (best_place == numeric_limits<double>::max())
		return false;
	else {
		if (replace) {
			/*cout << "Before Replace: " << endl;
			 PrintDepth();*/
			if (!offchip[index]) {
				for (i = place[index][1]; i < place[index][1] + place[index][3]; ++i) {
					for (j = place[index][2]; j < place[index][2] + place[index][4]; ++j) {
						depth[j][i] -= delay[index];
					}
				}
				/*cout << "Kernel " << index << " been replaced from " << place[index][1] << "\t" << place[index][2] << endl;
				 PrintDepth();*/
			}
			place[index][1] = new_x;
			place[index][2] = new_y;
			OTreeEncoding::kernels[index].SetCurrentImplementation(new_imp);
			place[index][3] = OTreeEncoding::kernels[index].GetWidth();
			place[index][4] = OTreeEncoding::kernels[index].GetHeight();
			place[index][5] = delay[index] = OTreeEncoding::kernels[index].GetDepth();
		}
		for (k = 0; k < numKernels; ++k) {
			if (offchip[k] || k == index)
				continue;
			if (Overlap(place[index][1], place[index][2], place[index][3], place[index][4], place[k][1], place[k][2],
					place[k][3], place[k][4])) {
				assert(!kernel_in_use[k]);
				conflict[index][k] = 1;
				conflict[k][index] = 1;
			}
		}

		for (i = place[index][1]; i < place[index][1] + place[index][3]; ++i) {
			for (j = place[index][2]; j < place[index][2] + place[index][4]; ++j) {
				depth[j][i] += delay[index];
			}
		}
		/*cout << "Kernel " << index << " been replaced to " << place[index][1] << "\t" << place[index][2] << "\t";
		 cout << place[index][3] << "\t" << place[index][4] << "\t" << place[index][5] << endl;
		 PrintDepth();*/

		return true;
	}

}

bool Scheduler::Overlap(int x1, int y1, int w1, int h1, int x2, int y2, int w2, int h2) {

	return ((x1 + w1 > x2 && x2 >= x1) || (x2 + w2 > x1 && x1 >= x2))
			&& ((y1 + h1 > y2 && y2 >= y1) || (y2 + h2 > y1 && y1 >= y2));
}

void Scheduler::PrintDepth() {
	int i, j;
	for (i = maxH - 1; i >= 0; --i) {
		for (j = 0; j < maxW; ++j) {
			cout << depth[i][j] << "\t";
		}
		cout << endl;
	}
}

} /* namespace std */
