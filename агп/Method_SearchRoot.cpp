#include "Method_SearchRoot.hpp"
#include<iostream>
#include <fstream>

void Method_SearchRoot::Solve() {
	Trial current;
	double Rmin;
	size_t Rpos;
	double curr_eps = trials[1].x - trials[0].x;
	std::vector<Trial>::iterator it = trials.begin();
	int itr = 0;
	//curr_eps = zmin;
	//std::ofstream out;
	//out.open("debug.txt", std::ofstream::ios_base::app);  // печать в файл
	double check = (check_method == 0) ? optimum.z : curr_eps;
	while (check > eps && itr <=10000) // optimum.z > eps curr_eps > eps
	{
		Rpos = 1;
		double R = trials[1].z * trials[0].z;
		Rmin = (R >= 0) ? R / (trials[1].x - trials[0].x) : 0;
		
		for (size_t i = 2; i < trials.size(); i++)
		{
			R = trials[i].z * trials[i - 1].z;
			R = (R >= 0) ? R / (trials[i].x - trials[i - 1].x) : 0;
			if (Rmin > R) {
				Rmin = R;
				Rpos = i;
			}
		}
		curr_eps = trials[Rpos].x - trials[Rpos - 1].x;

		std::vector<Trial>::iterator it2 = trials.begin();
		for (it = trials.begin(); it - trials.begin() <= Rpos; it++) it2 = it;

		current.x = (fabs(trials[Rpos].z)*trials[Rpos-1].x + fabs(trials[Rpos - 1].z) * trials[Rpos].x)
			/ (fabs(trials[Rpos].z) + fabs(trials[Rpos - 1].z));
		current.z = Funс(1, key, index_problem, current.x);
		trials.insert(it2, current);
		//out << "Iter = " << itr << "\ncurrent = {" << current.x << ", " << current.z << "}, size = " << trials.size();
		
		if (optimum.z > current.z) {
			//zpos = Rpos;
			best_i = itr;
			optimum = current;
		}
		//out << " opt = {" << optimum.x << ", " << optimum.z << "}\n";

		itr++;
		check = (check_method == 0) ? optimum.z : curr_eps;
	}
	//out.close();
	std::cout << "itr = " << itr << std::endl;
}