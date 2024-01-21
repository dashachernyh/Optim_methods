#include "MethodMult.h"
#include <iostream>
#include <chrono>
#include <thread>
#include <fstream>

// 0 - Grishagin, 1 - GKLS, 2 - HILL, 3 - Shekel
std::vector<double> MethodMult::GetTrueOpt()
{
	std::vector<double> optimumPoint;
	std::vector<double> point;
	if (task == 0 || task == 1) {
		optimumPoint = grish_fam[index_problem]->GetOptimumPoint();
		point.push_back(optimumPoint[0]);
		point.push_back(optimumPoint[1]);
	}
	else if(task == 2 || task == 3) {
		optimumPoint = gkls_fam[index_problem]->GetOptimumPoint();
		point.push_back(optimumPoint[0]);
		point.push_back(optimumPoint[1]);
	}
	else if  (task == 4) {
		optimumPoint = hill_fam[index_problem]->GetOptimumPoint();
		point.push_back(optimumPoint[0]);
	}
	else {
		optimumPoint = shekel_fam[index_problem]->GetOptimumPoint();
		point.push_back(optimumPoint[0]);
	}

	return point;
}
// 0 - Grishagin, 1 - GKLS, 2 - HILL, 3 - Shekel
void MethodMult::PrintTrueValue()
{
	vector<double> optimumPoint;
	if (task == 0 && task == 1) {
		optimumPoint = grish_fam[index_problem]->GetOptimumPoint();
		std::cout << "x1_min= " << optimumPoint[0] << " x2_min= " << optimumPoint[1];
		std::cout << "  F_min= " << grish_fam[index_problem]->GetOptimumValue() << std::endl;
	}
	else if (task == 2 && task == 3) {
		optimumPoint = gkls_fam[index_problem]->GetOptimumPoint();
		std::cout << "x1_min= " << optimumPoint[0] << " x2_min= " << optimumPoint[1];
		std::cout << "  F_min= " << gkls_fam[index_problem]->GetOptimumValue() << std::endl;
	}
	else if (task == 4) {
		optimumPoint = hill_fam[index_problem]->GetOptimumPoint();
		std::cout << "x_min= " << optimumPoint[0] << "  F_min= " << hill_fam[index_problem]->GetOptimumValue() << "\n";
	}
	else {
		optimumPoint = shekel_fam[index_problem]->GetOptimumPoint();
		std::cout << "x_min= " << optimumPoint[0] << "  F_min= " << shekel_fam[index_problem]->GetOptimumValue() << "\n";
	}
}
// 0 - Grishagin, 1 - GKLS,
double MethodMult::Funk_mult(double* y)
{
	vector<double> val;
	switch (task)
	{
	case 0: {
		val.push_back(y[0]);
		val.push_back(y[1]);
		return grish_fam[index_problem]->ComputeFunction(val);
	}
	case 1: {
		val.push_back(y[0]);
		val.push_back(y[1]);
		return fabs(grish_fam[index_problem]->ComputeFunction(val) - grish_fam[index_problem]->GetOptimumValue());
	}
	case 2: {
		val.push_back(y[0]);
		val.push_back(y[1]);
		return sqrt(fabs(grish_fam[index_problem]->ComputeFunction(val) - grish_fam[index_problem]->GetOptimumValue()));
	}
	case 3: {
		val.push_back(y[0]);
		val.push_back(y[1]);
		return fabs(-(grish_fam[index_problem]->ComputeFunction(val) - grish_fam[index_problem]->GetMaxValue()));
	}
	case 4: {
		val.push_back(y[0]);
		val.push_back(y[1]);
		return gkls_fam[index_problem]->ComputeFunction(val);
	}
	case 5: {
		val.push_back(y[0]);
		val.push_back(y[1]);
		return fabs(gkls_fam[index_problem]->ComputeFunction(val) - gkls_fam[index_problem]->GetOptimumValue());
	}
	case 6: {
		val.push_back(y[0]);
		val.push_back(y[1]);
		return sqrt(fabs(gkls_fam[index_problem]->ComputeFunction(val) - gkls_fam[index_problem]->GetOptimumValue()));
	}
	case 7: {
		val.push_back(y[0]);
		val.push_back(y[1]);
		return fabs(-(gkls_fam[index_problem]->ComputeFunction(val) - gkls_fam[index_problem]->GetOptimumValue()));
	}
	default:
		break;
	}
}

double MethodMult::Funk_multMat(double* y, std::vector<std::vector<int>>& matrix1,
	std::vector<std::vector<int>>& matrix2,
	std::vector<std::vector<std::vector<int>>>& matrix_res,
	int size, int index)
{
	//std::this_thread::sleep_for(std::chrono::milliseconds(100));
	for (int i = 0; i < size; ++i) {
		for (int j = 0; j < size; ++j) {
			matrix_res[index][i][j] = 0;
			for (int k = 0; k < size; ++k) {
				matrix_res[index][i][j] += matrix1[i][k] * matrix2[k][j];
			}
		}
	}

	vector<double> val(2);
	val[0] = y[0];
	val[1] = y[1];
	if (task == 0)
		return grish_fam[index_problem]->ComputeFunction(val);
	else
		return gkls_fam[index_problem]->ComputeFunction(val);
}

void MethodMult::ScaleFunc(double y)
{
	//�����
	if (b - a == 1)
	{
		koef[0] = 1;
		koef[1] = a - y;
	}

	else
	{
		koef[0] = b - a;
		koef[1] = a - koef[0] * y;
	}

	/*koef[0] = b - a;
	koef[1] = fabs(-0.5 - a);*/
}

void MethodMult::InsertScale(double* y)
{
	for (int i = 0; i < n; i++)
	{
		y[i] = y[i] * koef[0] + koef[1];
	}
}

void MethodMult::Init(int _task, int _index_problem, double* y, double _a, double _b, double _e, double _r,
	int _n, int _m) {
	a = _a;
	b = _b;
	eps = _e;
	r = _r;
	n = _n;
	m = _m;
	Trial current, first, second;
	task = _task;
	index_problem = _index_problem;
	best_i = 0;
	out_optimal.resize(n + 1);
	
	first.x = 0;
	mapd(first.x, m, y, n, 1);
	ScaleFunc(y[0]);
	InsertScale(y);
	first.z = Funk_mult(y); 
	trials.push_back(first);
	out_optimal[n] = first.z;
	for (int i = 0; i < n; ++i)
		out_optimal[i] = y[i];

	//std::mt19937 gen;
	//gen.seed(static_cast<unsigned int>(time(0)));
	//double x_1;
	//x_1 = gen() % 100 / (100 * 1.0);
	//if (x_1 != 1 && x_1 != 0)
	//{
	//	current.x = x_1;
	//	if (n != 1) {
	//		mapd(x_1, m, y, n, 1);
	//	} else {
	//		y[0] = current.x;
	//	}
	//	InsertScale(y);
	//	current.z = Funk_mult(y);
	//	trials.push_back(current);
	//}

	current.x = 0.5;
	mapd(current.x, m, y, n, 1);
	InsertScale(y);
	current.z = Funk_mult(y);
	trials.push_back(current);
	if (current.z < out_optimal[n]) {
		out_optimal[n] = current.z;
		for (int i = 0; i < n; ++i)
			out_optimal[i] = y[i];
	}

	second.x = 1;
	mapd(1, m, y, n, 1);
	InsertScale(y);
	second.z = Funk_mult(y);
	trials.push_back(second);

	if(second.z < out_optimal[n]) {
		out_optimal[n] = second.z;
		for (int i = 0; i < n; ++i)
			out_optimal[i] = y[i];
	}
}

void MethodMult::ClearMethod()
{
	trials.clear();
	out_optimal.clear();
	best_i = -1;
	koef[0] = koef[1] = -1;
}

void MethodMult::SolveMult(double* y, std::vector<std::vector<int>>& matrix1,
	std::vector<std::vector<int>>& matrix2,
	std::vector<std::vector<std::vector<int>>>& matrix_res,
	int size)
{
	Trial current;
	double M, Rmax;
	size_t Rpos;
	int itr = 0;
	std::vector<Trial>::iterator it = trials.begin();
	double z_min = out_optimal[n];
	double power = 1 / double(n);
	double curr_eps = pow(trials[1].x - trials[0].x, power);

	//std::ofstream out1;
	//out1.open("Grishagin.txt", std::ofstream::ios_base::app);

	std::vector<double> true_opt = GetTrueOpt();
	//while (curr_eps > eps) (fabs(true_opt[0] - out_optimal[0]) > eps || fabs(true_opt[1] - out_optimal[1]) > eps)
	// fabs(out_optimal[n])
	while ((fabs(true_opt[0] - out_optimal[0]) > eps || fabs(true_opt[1] - out_optimal[1]) > eps) && itr < 12000)
	{
		//out1 << "i = " << itr << " z_min = " << out_optimal[n] << "\n";
		Rpos = 1;

		double d_z = fabs(trials[1].z - trials[0].z);
		double d_x = fabs(trials[1].x - trials[0].x);
		d_x = pow(d_x, power);
		M = d_z / d_x;
		//M = 50;

		for (size_t i = 2; i < trials.size(); i++)
		{
			double max;
			max = fabs((trials[i].z - trials[i - size_t(1)].z)) / pow(trials[i].x - trials[i - size_t(1)].x, power);
			if (max > M)
				M = max;
		}
		if (M == 0)
			M = 1;

		//out1 << "\tM= " << M << "\n";

		Rmax = d_x + (pow(d_z / (r * M), 2) / d_x)
			- 2 * (trials[1].z + trials[0].z - 2 * z_min) / (r * M);

		double k, R;

		for (size_t i = 2; i < trials.size(); ++i)
		{
			k = pow(trials[i].x - trials[i - size_t(1)].x, power);

			R = k + (pow((trials[i].z - trials[i - size_t(1)].z) / (M * r), 2) / k) -
				2 * (trials[i].z + trials[i - size_t(1)].z - 2 * z_min) / (r * M);

			if (R > Rmax)
			{
				Rmax = R;
				Rpos = i;
			}
		}
		//out1 << "\tRmax = " << Rmax << " pos = " << Rpos << "\n";
		//curr_eps = pow(trials[Rpos].x - trials[Rpos - size_t(1)].x, power);

		// ����� �������� � �������
		std::vector<Trial>::iterator it2 = trials.begin();
		it2 += Rpos;

		double delta_z_pos = trials[Rpos].z - trials[Rpos - size_t(1)].z;

		int sgn = 0;
		if (delta_z_pos < 0)
			sgn = -1;
		if (delta_z_pos > 0)
			sgn = 1;

		current.x = (trials[Rpos].x + trials[Rpos - size_t(1)].x) / 2
			- sgn / (2 * r) * pow(fabs(delta_z_pos) / M, n);
		mapd(current.x, m, y, n, 1);
		InsertScale(y);

		current.z = Funk_mult(y);
		trials.insert(it2, current);

		//out1 << "\tx = " << current.x << " z = " << current.z << "\n";

		if (out_optimal[2] > current.z)
		{
			z_min = current.z;
			best_i = itr;
			out_optimal[0] = y[0];
			out_optimal[1] = y[1];
			out_optimal[2] = current.z;
		}

		itr++;
	}
	//std::cout << "itr = " << itr << std::endl;
	//out1.close();
}