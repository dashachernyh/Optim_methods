#include <iostream>
#include<fstream>
#include <random>
#include <ctime>
#include <algorithm>
#include <omp.h>
#include <vector>

#include "MethodMult_mixed.h"

MethodMult_mixed::MethodMult_mixed(int _index_problem, double* y, double _a, double _b, double _e, double _r,
	int _n, int _m, const int _on_step, Mixture _mix, int _alpha) : on_step(_on_step), mix(_mix), alpha(_alpha)
{
	a = _a;
	b = _b;
	eps = _e;
	r = _r;
	n = _n;
	m = _m;

	Trial current, first, second;
	index_problem = _index_problem;
	best_i = 0;
	out_optimal = { 0, 0, 0 };                  // ����������� ������� �������
	std::mt19937 gen;
	gen.seed(static_cast<unsigned int>(time(0)));

	first.x = 0;
	mapd(0, m, y, n, 1);
	ScaleFunc(y[0]);  // ��� ���� ����� ���� ��� �������
	InsertScale(y);
	first.z = Funk_mult(index_problem, y);  // ��������� �������� � ���� ����� 
	trials.push_back(first);

	current.x = 0.64;
	mapd(current.x, m, y, n, 1);
	InsertScale(y);
	current.z = Funk_mult(index_problem, y);
	trials.push_back(current);

	/*double x_1;
	x_1 = gen() % 100 / (100 * 1.0);         // �������� ������������ ����� ������ �� ��������� [0, 1]
	if (x_1 != 1 && x_1 != 0)              // ���� ��������� ����� x1 ����� ������ ��������� [0, 1]
	{
		current.x = x_1;
		mapd(x_1, m, y, n, 1);
		InsertScale(y);
		current.z = Funk_mult(index_problem, y);
		trials.push_back(current);
	}*/

	second.x = 1;
	mapd(1, m, y, n, 1);
	InsertScale(y);
	second.z = Funk_mult(index_problem, y);
	trials.push_back(second);
}

//�������� ��� ������������ ������, ��������� ������ y (����������), ����������� n, �������� �� [0, 1]
void MethodMult_mixed::SolveMult_mixed(double* y)
{
	Trial current;        // ��� �������� ������ ���������
	double M, Rmax;
	size_t Rpos;
	int flag = 0, k = 0;
	int change = mix.alg1;
	std::vector<Trial>::iterator it = trials.begin();  // ��� ������ ������� ���������� ������ ��������� � ������� trials
	double z_min = trials[0].z;       // ����������� �������� �������
	for (size_t i = 1; i < trials.size(); i++) {
		if (z_min > trials[i].z)
			z_min = trials[i].z;
	}
	int itr = 0;           // ������� ��������
	double power = 1 / double(n);
	double curr_eps = pow(trials[1].x - trials[0].x, power);
	out_optimal[2] = trials[0].z;

	std::ofstream out1;
	//out1.open("Grishagin.txt", std::ofstream::ios_base::app);  // ������ � ����
	std::vector<double> true_opt = GetTrueOpt_grish(index_problem);

	while (fabs(true_opt[0] - out_optimal[0]) > eps || fabs(true_opt[1] - out_optimal[1]) > eps)
	//while (curr_eps > eps)
	{
		if (itr >= on_step) {
			if (k == change)
			{
				if (!flag) {
					flag = 1;
					change = mix.alg2;
				}
				else
				{
					flag = 0;
					change = mix.alg1;
				}
				k = 0;
			}
		}
		Rpos = 1;
		
		// ��������� �������� ��� ���������� M �� ��������� 1
		double d_z = fabs(trials[1].z - trials[0].z);
		double d_x = fabs(trials[1].x - trials[0].x);
		d_x = pow(d_x, power);
		M = d_z / d_x;

		for (size_t i = 2; i < trials.size(); i++)      // ����� �� 2 ���������
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

		if (itr >= on_step && flag)  // ��������� ��������-����������
		{
			double R = d_x + (pow(d_z / (r * M), 2) / d_x)
				- 2 * (trials[1].z + trials[0].z - 2 * z_min) / (r * M);
			double z = (trials[1].z - z_min) * (trials[0].z - z_min);
			Rmax = R / ((pow(z, power) / M) + 1 / pow(1.5, alpha));
			for (size_t i = 2; i < trials.size(); i++)
			{
				// ��������� ������������� ��� ���
				double k = pow(trials[i].x - trials[i - size_t(1)].x, power);
				R = k + (pow((trials[i].z - trials[i - size_t(1)].z) / (M * r), 2) / k) -
					2 * (trials[i].z + trials[i - size_t(1)].z - 2 * z_min) / (r * M);
				// ��������-����������
				z = (trials[i].z - z_min) * (trials[i - size_t(1)].z - z_min);
				double r_alpha = R / (pow(z, power) / M + 1 / pow(1.5, alpha));

				if (r_alpha > Rmax)
				{
					Rmax = r_alpha;
					Rpos = i;
				}
			}

		}
		else  // ������� ���
		{
			// �������� �� ��������� 1
			Rmax = d_x + (pow(d_z / (r * M), 2) / d_x)
				- 2 * (trials[1].z + trials[0].z - 2 * z_min) / (r * M);

			for (size_t i = 2; i < trials.size(); i++)       // ����� �� 2 ���������
			{
				double k = pow(trials[i].x - trials[i - size_t(1)].x, power);   // ��� ����������� ����������

				double R = k + (pow((trials[i].z - trials[i - size_t(1)].z) / (M * r), 2) / k) -
					2 * (trials[i].z + trials[i - size_t(1)].z - 2 * z_min) / (r * M);

				if (R > Rmax)
				{
					Rmax = R;
					Rpos = i;
				}
			}
		}

		curr_eps = pow(trials[Rpos].x - trials[Rpos - size_t(1)].x, power);

		// ����� �������� � �������
		std::vector<Trial>::iterator it2 = trials.begin();
		for (it = trials.begin(); it - trials.begin() <= Rpos; it++) it2 = it;

		double delta_z_pos = trials[Rpos].z - trials[Rpos - size_t(1)].z;

		int sgn = 0;
		if (delta_z_pos < 0)
			sgn = -1;
		if (delta_z_pos > 0)
			sgn = 1;

		current.x = (trials[Rpos].x + trials[Rpos - size_t(1)].x) / 2 - sgn / (2 * r) * pow(fabs(delta_z_pos) / M, n);
		// ���������� ���������
		mapd(current.x, m, y, n, 1);
		InsertScale(y);

		current.z = Funk_mult(index_problem, y);  // �������� ������� � ������ y
		trials.insert(it2, current);


		//out1 << y[0] << " " << y[1] << std::endl;

		if (out_optimal[2] > current.z)
		{
			best_i = itr;
			out_optimal[0] = y[0];
			out_optimal[1] = y[1];
			out_optimal[2] = current.z;
		}
		if (itr >= on_step) k++;
		itr++;
		
		/*std::cout << itr << " pos = "<<Rpos<< std::endl;
		for (int i = 0; i < trials.size(); i++)
		{
			std::cout << trials[i].x << " " << trials[i].z << std::endl;
		}*/
		//out1.close();
	}
	//out1.open("Grishagin.txt", std::ofstream::ios_base::app);
	std::cout << "itr = " << itr << std::endl;
	//out1 << trials.size() << std::endl;
	//out1.close();
}