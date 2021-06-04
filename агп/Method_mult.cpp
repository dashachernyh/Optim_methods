#include <iostream>
#include<fstream>
#include <random>
#include <ctime>

#include"Trial.h"
#include "Method_mult.h"
#include "Map.h"

std::vector<double> MethodMult::GetTrueOpt_grish(int index_problem)
{	
	std::vector<double> optimumPoint;
	optimumPoint = grishFam[index_problem]->GetOptimumPoint();
	std::vector<double> point{ optimumPoint[0], optimumPoint[1] };
	return point;
}

void MethodMult:: PrintTrueValueGrishagin(int index_problem)
{
	vector<double> optimumPoint;
	optimumPoint = grishFam[index_problem]->GetOptimumPoint();
	std::cout << "x1_min= " << optimumPoint[0] << " x2_min= " << optimumPoint[1] << "  F_min= " << grishFam[index_problem]->GetOptimumValue() << std::endl;
}

// вычисляет значение задачи Гришагина 
double MethodMult::Funk_mult(int index_problem, double* y)
{
	vector<double> val;
	val.push_back(y[0]);
	val.push_back(y[1]);
	return grishFam[index_problem]->ComputeFunction(val);
}

// масштабирует область поиска (куб)
void MethodMult::ScaleFunc(double y)
{
	//сдвиг
	if (b - a == 1)       //длина интервала не изменилась
	{
		koef[0] = 1;
		koef[1] = a - y;
	}
	// масштаб
	else
	{
		koef[0] = b - a;
		koef[1] = a - koef[0] * y;
	}
}

// применение масштабирование к y
void MethodMult::InsertScale(double* y)
{
	// y имеет размерность пространства n
	for (int i = 0; i < n; i++)
	{
		y[i] = y[i] * koef[0] + koef[1];
	}
}

MethodMult::MethodMult(int _index_problem, double *y, double _a, double _b, double _e, double _r,
	double _n, double _m) :a(_a), b(_b), eps(_e), r(_r), n(_n), m(_m)
{
	Trial current;
	index_problem = _index_problem;
	best_i = 0;
	out_optimal = { 0, 0, 0 };                  // оптимальное решение нулевое
	std::mt19937 gen;
	gen.seed(static_cast<unsigned int>(time(0)));
	double x_1;
	x_1 = gen() % 100 / (100 * 1.0);         // выбираем произвольную точку поиска на интервале [0, 1]
	first.x = 0;                                                    
	mapd(0, m, y, n, 1);                        
	ScaleFunc(y[0]);
	InsertScale(y);
	first.z = Funk_mult(index_problem, y);  // вычисляем значение в этой точке 
	trials.push_back(first);

	if (x_1 != 1 && x_1 != 0)              // если выбранная точка x1 лежит внутри интервала [0, 1]
	{
		current.x = x_1;
		mapd(x_1, m, y, n, 1);
		InsertScale(y);
		current.z = Funk_mult(index_problem, y);
		trials.push_back(current);
	}
	second.x = 1;
	mapd(1, m, y, n, 1);
	InsertScale(y);
	second.z = Funk_mult (index_problem, y);
	trials.push_back(second);
}

//алгоритм для многомерного случая, принимает массив y (координаты), размерности n, значения на [0, 1]
void MethodMult::solve_mult(double* y)
{
	Trial current;        // для подсчета нового испытания
	double M, Rmax;
	int Rpos;
	std::vector<Trial>::iterator it = trials.begin();  // для поиска позиции добавления нового испытания в векторе trials
	double z_min = trials[0].z;                        // минимальное значение функции
	for (int i = 1; i < trials.size(); i++) {
		if (z_min > trials[i].z)
			z_min = trials[i].z;
	}
	int itr = 0;           // счетчик итераций
	double power = 1 / n;  
	double curr_eps = pow(trials[1].x - trials[0].x, power);
	out_optimal[2] = first.z;

	//std::ofstream out1;
	while (curr_eps > eps)
	{
		Rpos = 1;
		// начальные значение для вычисления M на интервале 1
		double d_z = fabs(trials[1].z - trials[0].z);   
		double d_x = fabs(trials[1].x - trials[0].x);
		d_x = pow(d_x, power);
		M = d_z / d_x;

		for (int i = 2; i < trials.size(); i++)      // поиск со 2 интервала
		{
			double max;
			max = fabs((trials[i].z - trials[i - 1].z)) / pow(trials[i].x - trials[i - 1].x, power);
			if (max > M)
				M = max;
			if (z_min > trials[i].z)
				z_min = trials[i].z;
		}
		if (M == 0)
			M = 1;

		// значение на интервале 1
		Rmax = d_x + (pow(d_z / (r * M), 2) / d_x)
			- 2* (trials[1].z + trials[0].z - 2 * z_min) / (r * M);

		for (int i = 2; i < trials.size(); i++)       // поиск со 2 интервала
		{
			double k = pow(trials[i].x - trials[i - 1].x, power);   // для оптимизации вычисления

			double R = k + (pow((trials[i].z - trials[i - 1].z) / (M * r), 2) / k) -
				2 * (trials[i].z + trials[i - 1].z - 2 * z_min) / (r * M);

			if (R > Rmax)
			{
				Rmax = R;
				Rpos = i;
			}
		}
		curr_eps = pow(trials[Rpos].x - trials[Rpos - 1].x, power);

		// поиск поизиции в массиве
		std::vector<Trial>::iterator it2 = trials.begin();        
		for (it = trials.begin(); it - trials.begin() <= Rpos; it++) it2 = it;

		double delta_z_pos = trials[Rpos].z - trials[Rpos - 1].z;

		int sgn = 0;
		if (delta_z_pos < 0)
			sgn = -1;
		if (delta_z_pos > 0)
			sgn = 1;

		current.x = (trials[Rpos].x + trials[Rpos - 1].x) / 2 - sgn / (2 * r) * pow(delta_z_pos / M, n) ;
		// приведение координат
		mapd(current.x, m, y, n, 1);
		InsertScale(y);

		current.z = Funk_mult(index_problem, y);  // значении функции в точках y
		trials.insert(it2, current);

		//out1.open("Grishagin.txt", std::ofstream::ios_base::app);  // печать в файл
		//out1 << y[0] << " " << y[1] << std::endl;

		if (out_optimal[2] > current.z)
		{
			best_i = itr;
			out_optimal[0] = y[0];
			out_optimal[1] = y[1];
			out_optimal[2] = current.z;
		}
		itr++;
		//out1.close();
	}
	std::cout << "itr = " << itr << std::endl;
}