#include"MethodMult_DualLipsh.h"

#include <iostream>
#include <fstream>

void MethodMult_DualLipsh::SolveMult_DualLipsh(double* y)
{
	Trial current;  // для подсчета нового испытания
	double M, R_loc, R_glob;
	MaxR Rmax;
	int itr = 0;  // счетчик итераций
	std::vector<Trial>::iterator it = trials.begin();  // для поиска позиции добавления нового испытания в векторе trials

	// минимальное значение функции
	double z_min = trials[0].z;                        
	for (size_t i = 1; i < trials.size(); i++) {
		if (z_min > trials[i].z)
			z_min = trials[i].z;
	}

	double power = 1 / double(n);
	double curr_eps = pow(trials[1].x - trials[0].x, power);
	out_optimal[2] = trials[0].z;

	double ro = (1 - 1 / r_glob) / (1 - 1 / r_loc) * (1 - 1 / r_glob) / (1 - 1 / r_loc);

	// печать в файл
	//std::ofstream out1;
	//out1.open("Grishagin.txt", std::ofstream::ios_base::app);

	std::vector<double> true_opt = GetTrueOpt(task, index_problem);

	while (fabs(true_opt[0] - out_optimal[0]) > eps || fabs(true_opt[1] - out_optimal[1]) > eps)  //while (curr_eps > eps)
	{
		Rmax.pos = 1;

		// начальные значение для вычисления M на интервале 1
		double d_z = fabs(trials[1].z - trials[0].z);
		double d_x = fabs(trials[1].x - trials[0].x);
		d_x = pow(d_x, power);
		M = d_z / d_x;  //M = 50;

		for (size_t i = 2; i < trials.size(); i++)  // поиск со 2 интервала
		{
			double max;
			max = fabs((trials[i].z - trials[i - size_t(1)].z)) / pow(trials[i].x - trials[i - size_t(1)].x, power);
			if (max > M)
				M = max;
			if (z_min > trials[i].z)
				z_min = trials[i].z;
		}
		if (M == 0)
			M = 1;

		// значение на интервале 1

		R_loc = d_x + (pow(d_z / (r_loc * M), 2) / d_x)- 2 * (trials[1].z + trials[0].z - 2 * z_min)/ (r_loc * M);
		R_glob = d_x + (pow(d_z / (r_glob * M), 2) / d_x) - 2 * (trials[1].z + trials[0].z - 2 * z_min) / (r_glob * M);

		if (ro * R_loc > R_glob)
		{
			Rmax.R = ro * R_loc;
			Rmax.r = r_loc;
		}
		else
		{
			Rmax.R = R_glob;
			Rmax.r = r_glob;
		}

		double k, R;

		for (size_t i = 2; i < trials.size(); i++)  // поиск со 2 интервала
		{
			k = pow(trials[i].x - trials[i - size_t(1)].x, power);  // для оптимизации вычисления

			R_loc = k + (pow((trials[i].z - trials[i - size_t(1)].z) / (M * r_loc), 2) / k) -
				2 * (trials[i].z + trials[i - size_t(1)].z - 2 * z_min) / (r_loc * M);
			R_glob = k + (pow((trials[i].z - trials[i - size_t(1)].z) / (M * r_glob), 2) / k) -
				2 * (trials[i].z + trials[i - size_t(1)].z - 2 * z_min) / (r_glob * M);
			
			if (ro * R_loc > R_glob)
			{
				R = ro * R_loc;
				r = r_loc;
			}
			else
			{
				R = R_glob;
				r = r_glob;
			}

			if (R > Rmax.R)
			{
				Rmax.R = R;
				Rmax.pos = i;
				Rmax.r = r;
			}
		}

		//curr_eps = pow(trials[Rmax.pos].x - trials[Rmax.pos - size_t(1)].x, power);

		// поиск поизиции в массиве
		std::vector<Trial>::iterator it2 = trials.begin();
		for (it = trials.begin(); it - trials.begin() <= Rmax.pos; it++) it2 = it;

		double delta_z_pos = trials[Rmax.pos].z - trials[Rmax.pos - size_t(1)].z;

		int sgn = 0;
		if (delta_z_pos < 0)
			sgn = -1;
		if (delta_z_pos > 0)
			sgn = 1;

		current.x = (trials[Rmax.pos].x + trials[Rmax.pos - size_t(1)].x) / 2 - sgn / (2 * Rmax.r) * pow(fabs(delta_z_pos) / M, n);
		mapd(current.x, m, y, n, 1);  // приведение координат
		InsertScale(y);

		current.z = Funk_mult(task, index_problem, y);  // значении функции в точках y
		trials.insert(it2, current);

		//out1 << y[0] << " " << y[1] << std::endl;
		
		if (out_optimal[2] > current.z)
		{
			best_i = itr;
			out_optimal[0] = y[0];
			out_optimal[1] = y[1];
			out_optimal[2] = current.z;
		}

		itr++;
	}

	std::cout << "itr = " << itr << std::endl;
	//out1 <<trials.size() << std::endl;
	//out1.close();
}