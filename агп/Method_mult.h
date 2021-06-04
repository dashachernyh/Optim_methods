#pragma once
#include "GrishaginOption.hpp"
#include "GrishaginProblemFamily.hpp"
#include <vector>

class MethodMult
{
	std::vector<Trial> trials;        // вектор испытаний
	Trial first;                      // испытание в первой точке
	Trial second;                     // во второй 
	std::vector<double> out_optimal;  // значение out_optimal[0] = y*, out_optimal[1] = z* 
	int index_problem;
	int best_i;
	double n, m;
	double eps, r;
	double a, b;                       //интервал куба
	double koef[2];                    // коэффициенты масштабирования koef[0] для сдвига, koef[1] для увеличения масштаба
	TGrishaginProblemFamily grishFam;  //Семейство задач
public:                         
	std::vector<double> GetOpt() { return out_optimal; }  // возвращает оптимальное значение
	int GetBestIndex() { return best_i; }                 
	MethodMult(int _index_problem, double* y, double _a,  double _b, double _e, double _r, double _n, double _m);
	void solve_mult(double * y);
	void ScaleFunc(double y);                             // масштабирует область поиска (куб)
	void InsertScale(double* y);                          // применение масштабирование к y
	void PrintTrueValueGrishagin(int index_problem); 
	double Funk_mult(int index_problem, double* y);
	std::vector<double> GetTrueOpt_grish(int index_problem);
};
