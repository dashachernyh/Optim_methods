#pragma once
#include "GrishaginOption.hpp"
#include "GrishaginProblemFamily.hpp"
#include "GKLSOption.hpp"
#include "GKLSProblemFamily.hpp"
#include "Trial.h"
#include "Map.h"

#include <random>
#include <ctime>
#include <algorithm>

class MethodMult
{
protected:
	std::vector<Trial> trials;        // вектор испытаний
	std::vector<double> out_optimal;  // значение out_optimal[0] = y*, out_optimal[1] = z* 
	int index_problem;
	int best_i;
	int n, m;
	double eps, r;
	double a, b;                       //интервал куба
	double koef[2];                    // коэффициенты масштабирования koef[0] для сдвига, koef[1] для увеличения масштаба
	TGrishaginProblemFamily grish_fam;  //Семейство задач
	TGKLSProblemFamily gkls_fam;
	int task;                          // 0 - Гришагин, 1 - Сергеев
public:                         
	MethodMult()
	{
		a = 0;
		b = 0;
		n = 0;
		m = 0;
		eps = 0;
		r = 0;
		koef[0] = koef[1] = 0;
		index_problem = best_i = task = 0;
	}
	MethodMult(int _task, int _index_problem, double* y, double _a,  double _b, double _e, double _r, int _n, int _m);
	std::vector<double> GetOpt() { return out_optimal; }  // возвращает оптимальное значение
	int GetBestIndex() { return best_i; }
	void SolveMult(double * y);
	void ScaleFunc(double y);                             // масштабирует область поиска (куб)
	void InsertScale(double* y);                          // применение масштабирование к y
	void PrintTrueValue(int task, int index_problem); 
	double Funk_mult(int task, int index_problem, double* y);
	std::vector<double> GetTrueOpt(int task, int index_problem);
	double Funk_test(int index_problem, double* y);
};
