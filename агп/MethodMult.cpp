#include "MethodMult.h"

#include <iostream>
#include<fstream>

double MethodMult::Funk_test(int index_problem, double* y)
{
	return y[0]*y[0]+y[1]*y[1] - 2;  // гиперболойд элиптический мин в т (0, 0, -2)
}

std::vector<double> MethodMult::GetTrueOpt(int task, int index_problem)
{	
	std::vector<double> optimumPoint;
	if (task == 0) {
		optimumPoint = grish_fam[index_problem]->GetOptimumPoint();
	}
	else
	{
		optimumPoint = gkls_fam[index_problem]->GetOptimumPoint();
	}
	std::vector<double> point{ optimumPoint[0], optimumPoint[1] };
	return point;
}

void MethodMult::PrintTrueValue(int task, int index_problem)
{
	vector<double> optimumPoint;
	if (task == 0) {
		optimumPoint = grish_fam[index_problem]->GetOptimumPoint();
		std::cout << "x1_min= " << optimumPoint[0] << " x2_min= " << optimumPoint[1];
		std::cout << "  F_min= " << grish_fam[index_problem]->GetOptimumValue() << std::endl;
	}
	else
	{
		optimumPoint = gkls_fam[index_problem]->GetOptimumPoint();
		std::cout << "x1_min= " << optimumPoint[0] << " x2_min= " << optimumPoint[1];
		std::cout << "  F_min= " << gkls_fam[index_problem]->GetOptimumValue() << std::endl;
	}
}

// вычисл€ет значение задачи √ришагина/—ергеева 
double MethodMult::Funk_mult(int task, int index_problem, double* y)
{
	vector<double> val;
	val.push_back(y[0]);
	val.push_back(y[1]);
	if (task == 0)
		return grish_fam[index_problem]->ComputeFunction(val);
	else 
		return gkls_fam[index_problem]->ComputeFunction(val);
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

MethodMult::MethodMult(int _task, int _index_problem, double *y, double _a, double _b, double _e, double _r,
	int _n, int _m) :a(_a), b(_b), eps(_e), r(_r), n(_n), m(_m)
{
	Trial current, first, second;
	task = _task;
	index_problem = _index_problem;
	best_i = 0;
	out_optimal = { 0, 0, 0 };  // оптимальное решение
	std::mt19937 gen;
	gen.seed(static_cast<unsigned int>(time(0)));

	first.x = 0;                                                    
	mapd(0, m, y, n, 1); 
	ScaleFunc(y[0]);  // дл€ всех точек один раз считаем
	InsertScale(y);
	first.z = Funk_mult(task, index_problem, y);  // вычисл€ем значение в этой точке 
	trials.push_back(first);
	
	// 0.64 0.38 0.5 0.63 0. 91
	/*current.x = 0.5;
	mapd(current.x, m, y, n, 1);
	InsertScale(y);
	current.z = Funk_mult(task, index_problem, y);
	trials.push_back(current);*/

	/*current.x = 0.54;
	mapd(current.x, m, y, n, 1);
	InsertScale(y);
	current.z = Funk_mult(task, index_problem, y);
	trials.push_back(current);

	current.x = 0.6;
	mapd(current.x, m, y, n, 1);
	InsertScale(y);
	current.z = Funk_mult(task, index_problem, y);
	trials.push_back(current);*/

	

	double x_1;
	x_1 = gen() % 100 / (100 * 1.0);  // выбираем произвольную точку поиска на интервале [0, 1]
	if (x_1 != 1 && x_1 != 0)         // если выбранна€ точка x1 лежит внутри интервала [0, 1]
	{
		current.x = x_1;
		mapd(x_1, m, y, n, 1);
		InsertScale(y);
		current.z = Funk_mult(task, index_problem, y);
		trials.push_back(current);
	}
	
	second.x = 1;
	mapd(1, m, y, n, 1);
	InsertScale(y);
	second.z = Funk_mult(task, index_problem, y);
	trials.push_back(second);
}

//алгоритм дл€ многомерного случа€, принимает массив y (координаты), размерности n, значени€ на [0, 1]
void MethodMult::SolveMult(double* y)
{
	Trial current;  // дл€ подсчета нового испытани€
	double M, Rmax;
	size_t Rpos;
	int itr = 0;  // счетчик итераций
	std::vector<Trial>::iterator it = trials.begin();  // дл€ поиска позиции добавлени€ нового испытани€ в векторе trials
	double z_min = trials[0].z;                        // минимальное значение функции
	for (size_t i = 1; i < trials.size(); i++) {
		if (z_min > trials[i].z)
			z_min = trials[i].z;
	}
	double power = 1 / double(n);  
	double curr_eps = pow(trials[1].x - trials[0].x, power);
	out_optimal[2] = trials[0].z;

	// печать в файл
	std::ofstream out1;
	out1.open("Grishagin.txt", std::ofstream::ios_base::app);
	

	std::vector<double> true_opt = GetTrueOpt(task, index_problem);

	 while (fabs(true_opt[0] - out_optimal[0]) > eps || fabs(true_opt[1] - out_optimal[1]) > eps)//while (curr_eps > eps)
	{
		Rpos = 1;
		
		// начальные значение дл€ вычислени€ M на интервале 1
		double d_z = fabs(trials[1].z - trials[0].z);   
		double d_x = fabs(trials[1].x - trials[0].x);
		d_x = pow(d_x, power);
		M = d_z / d_x;
		//M = 50;

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
		Rmax = d_x + (pow(d_z / (r * M), 2) / d_x)
			- 2* (trials[1].z + trials[0].z - 2 * z_min) / (r * M);
	
		double k, R;

		for (size_t i = 2; i < trials.size(); i++)  // поиск со 2 интервала
		{
			k = pow(trials[i].x - trials[i - size_t(1)].x, power);  // дл€ оптимизации вычислени€

			R = k + (pow((trials[i].z - trials[i - size_t(1)].z) / (M * r), 2) / k) -
				2 * (trials[i].z + trials[i - size_t(1)].z - 2 * z_min) / (r * M);

			if (R > Rmax)
			{
				Rmax = R;
				Rpos = i;
			}
		}

		//curr_eps = pow(trials[Rpos].x - trials[Rpos - size_t(1)].x, power);
	
		// поиск поизиции в массиве
		std::vector<Trial>::iterator it2 = trials.begin();        
		for (it = trials.begin(); it - trials.begin() <= Rpos; it++) it2 = it;

		double delta_z_pos = trials[Rpos].z - trials[Rpos - size_t(1)].z;

		int sgn = 0;
		if (delta_z_pos < 0)
			sgn = -1;
		if (delta_z_pos > 0)
			sgn = 1;

		current.x = (trials[Rpos].x + trials[Rpos - size_t(1)].x) / 2 - sgn / (2 * r) * pow(fabs(delta_z_pos) / M, n) ;
		mapd(current.x, m, y, n, 1);  // приведение координат
		InsertScale(y);

		current.z = Funk_mult(task, index_problem, y);  // значении функции в точках y
		trials.insert(it2, current);

		out1 << y[0] << " " << y[1] << std::endl;
		//out1 << "itr= " << itr <<  " y0 = "<<y[0]<<" y1 = "<<y[1]<<" z= " << current.z << std::endl;

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
	out1.close();
}