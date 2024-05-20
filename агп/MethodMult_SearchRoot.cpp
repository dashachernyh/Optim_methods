#include "MethodMult_SearchRoot.hpp"
#include<iostream>
#include <fstream>

void MethodMult_SearchRoot::SolveMult_SR(double* y) {

	Trial current;
	double Rmin;
	size_t Rpos;
	std::vector<Trial>::iterator it = trials.begin();
	int itr = 0;
	double power = 1 / double(n);
	double curr_eps = pow(trials[1].x - trials[0].x, power);
	std::vector<double> true_opt = GetTrueOpt();
	std::cout << "true_opt " << true_opt[0] <<" " << true_opt[1]<< std::endl;

	std::ofstream out1;
	out1.open("SearchRoot.txt", std::ofstream::ios_base::app);
		//curr_eps = zmin;
	//std::ofstream out;
	//out.open("debug.txt", std::ofstream::ios_base::app);  // печать в файл
	
	//while (curr_eps > eps) z_min > eps curr_eps > eps
	//out1 << true_opt[0] << " " << true_opt[1];
	bool check = (check_method == 0) ? out_optimal[n] > eps : (fabs(true_opt[0] - out_optimal[0]) > eps || fabs(true_opt[1] - out_optimal[1]) > eps);
	while (check && itr < 20000)
	{
		
		Rpos = 1;
		double R = trials[1].z * trials[0].z;
		Rmin = (R >= 0) ? R / pow(trials[1].x - trials[0].x, power) : 0;

		for (size_t i = 2; i < trials.size(); i++)
		{
			R = trials[i].z * trials[i - 1].z;
			R = (R >= 0) ? R / pow(trials[i].x - trials[i - 1].x, power) : 0;
			if (Rmin > R) {
				Rmin = R;
				Rpos = i;
			}
		}
		curr_eps = pow(trials[Rpos].x - trials[Rpos - 1].x, power);

		std::vector<Trial>::iterator it2 = trials.begin();
		for (it = trials.begin(); it - trials.begin() <= Rpos; it++) it2 = it;

		current.x = (fabs(trials[Rpos].z) * trials[Rpos - 1].x + fabs(trials[Rpos - 1].z) * trials[Rpos].x)
			/ (fabs(trials[Rpos].z) + fabs(trials[Rpos - 1].z));
		mapd(current.x, m, y, n, 1);
		InsertScale(y);

		current.z = Funk_mult(y);
		trials.insert(it2, current);

		//out1 << y[0] << " " << y[1] << "\n";

		if (out_optimal[2] > current.z) {
			best_i = itr;
			out_optimal[0] = y[0];
			out_optimal[1] = y[1];
			out_optimal[2] = current.z;
		}
		//out << " opt = {" << optimum.x << ", " << optimum.z << "}\n";
		//std::cout << "itr= " << itr << "\n";
		++itr;
		check = (check_method == 0) ? out_optimal[n] > eps : (fabs(true_opt[0] - out_optimal[0]) > eps || fabs(true_opt[1] - out_optimal[1]) > eps);
	}
	out1.close();
	std::cout << "itr = " << itr << std::endl;
}