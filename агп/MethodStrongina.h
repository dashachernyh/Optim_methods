#pragma once
#include <algorithm>
#include<vector>
#include <math.h>

#include"Trial.h"
#include "HansenProblem.hpp"
#include "HansenProblemFamily.hpp"
#include "HillProblem.hpp"
#include "HillProblemFamily.hpp"

class Method
{
protected:
	std::vector<Trial> trials;
	Trial optimum;
	Trial first;
	Trial second;
	int index_problem;
	int best_i;
	int key;
	double eps, r;
	THillProblemFamily hillFam;
	THansenProblemFamily hansFam;
public:
	Method() {
		Trial zero;
		zero.x = 0;
		zero.z = 0;
		optimum = first = second = zero;
		key = 0;
		index_problem = best_i= 0;
		r = eps = 0;
	};

	Method(int _key, int index_problem, std::vector<double> x_0, std::vector<double> x_n, double _e, double _r);
	Trial GetOpt() { return optimum; }
	int GetBestIndex() { return best_i; }
	void InitIntervalHans(int index, vector<double> &begin, vector<double> &end);
	void Solve();
	void  PrintTrueValueHans(int index_problem);
	void PrintTrueValueHill(int index_problem);
	double GetTrueOpt_hans(int index_problem);
	double GetTrueOpt_hill(int index_problem);
	double Funk(int key, int index_problem, double x);
};
