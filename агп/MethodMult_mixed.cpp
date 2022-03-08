#include "MethodMult_mixed.h"

#include <iostream>
#include<fstream>

MethodMult_mixed::MethodMult_mixed(int _task, int _index_problem, double* y, double _a, double _b, double _e,
	double _r, int _n, int _m, int _on_step, Mixture _mix, int _alpha)
	:MethodMult(_task, _index_problem, y, _a, _b, _e, _r, _n, _m),
	on_step(_on_step), mix(_mix), alpha(_alpha) {}

//алгоритм для многомерного случая, принимает массив y (координаты), размерности n, значения на [0, 1]
void MethodMult_mixed::SolveMult_mixed(double* y, int a)
{
	Trial current;        // для подсчета нового испытания
	double M, Rmax;
	size_t Rpos;
	int flag = 0, count = 0;
	int itr = 0;           // счетчик итераций
	int change = mix.alg1;
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

	//std::ofstream out1;
	//out1.open("Grishagin.txt", std::ofstream::ios_base::app);  // печать в файл
	std::vector<double> true_opt = GetTrueOpt(task, index_problem);

	while (fabs(true_opt[0] - out_optimal[0]) > eps || fabs(true_opt[1] - out_optimal[1]) > eps)  //while (curr_eps > eps)
	{
		if (itr >= on_step) {
			if (count == change)
			{
				if (!flag && mix.alg2 != 0) {
					flag = 1;
					change = mix.alg2;
				}
				else
				{
					flag = 0;
					change = mix.alg1;
				}
				count = 0;
			}
		}
		Rpos = 1;

		// начальные значение для вычисления M на интервале 1
		double d_z = fabs(trials[1].z - trials[0].z);
		double d_x = fabs(trials[1].x - trials[0].x);
		d_x = pow(d_x, power);
		M = d_z / d_x;

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

		double k,R;

		if (itr >= on_step && flag)  // выполняем локально-адаптивный
		{
			R = d_x + (pow(d_z / (r * M), 2) / d_x)
				- 2 * (trials[1].z + trials[0].z - 2 * z_min) / (r * M);
			double z = (trials[1].z - z_min) * (trials[0].z - z_min);
			Rmax = R / ((pow(z, power) / M) + 1 / pow(1.5, alpha));

			for (size_t i = 2; i < trials.size(); i++)
			{
				// вычисляем характерисику для агп
				k = pow(trials[i].x - trials[i - size_t(1)].x, power);
				R = k + (pow((trials[i].z - trials[i - size_t(1)].z) / (M * r), 2) / k) -
					2 * (trials[i].z + trials[i - size_t(1)].z - 2 * z_min) / (r * M);

				// локально-адаптивная
				z = (trials[i].z - z_min) * (trials[i - size_t(1)].z - z_min);
				double r_alpha = R / (pow(z, power) / M + 1 / pow(1.5, alpha));

				if (r_alpha > Rmax)
				{
					Rmax = r_alpha;
					Rpos = i;
				}
			}

		}
		else  // обычный агп
		{
			// значение на интервале 1
			Rmax = d_x + (pow(d_z / (r * M), 2) / d_x)
				- 2 * (trials[1].z + trials[0].z - 2 * z_min) / (r * M);

			for (size_t i = 2; i < trials.size(); i++) // поиск со 2 интервала
			{
				k = pow(trials[i].x - trials[i - size_t(1)].x, power);   // для оптимизации вычисления

				R = k + (pow((trials[i].z - trials[i - size_t(1)].z) / (M * r), 2) / k) -
					2 * (trials[i].z + trials[i - size_t(1)].z - 2 * z_min) / (r * M);

				if (R > Rmax)
				{
					Rmax = R;
					Rpos = i;
				}
			}
		}

		//curr_eps = pow(trials[Rpos].x - trials[Rpos - size_t(1)].x, power);

		// поиск поизиции в массиве
		std::vector<Trial>::iterator it2 = trials.begin();
		for (it = trials.begin(); size_t(it - trials.begin()) <= Rpos; it++) it2 = it;

		double delta_z_pos = trials[Rpos].z - trials[Rpos - size_t(1)].z;

		int sgn = 0;
		if (delta_z_pos < 0)
			sgn = -1;
		if (delta_z_pos > 0)
			sgn = 1;

		current.x = (trials[Rpos].x + trials[Rpos - size_t(1)].x) / 2 - sgn / (2 * r) * pow(fabs(delta_z_pos) / M, n);
		mapd(current.x, m, y, n, 1);
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
		if (itr >= on_step) count++;
		itr++;
	}

	std::cout << "itr = " << itr << std::endl;
	//out1.open("Grishagin.txt", std::ofstream::ios_base::app);
	//out1 << trials.size() << std::endl;
	//out1.close();
}