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
	double eps, r;
	THillProblemFamily hillFam;
	THansenProblemFamily hansFam;
public:
	int key;
	Method() {};
	Method(int _key, int index_problem, double x_0, double x_n, double _e, double _r);
	Trial GetOpt() { return optimum; }
	int GetBestIndex() { return best_i; }
	void InitIntervalHans(int index, vector<double> &begin, vector<double> &end);
	void solve();
	void  PrintTrueValueHans(int index_problem);
	void PrintTrueValueHill(int index_problem);
protected:
	double Funk(int key, int index_problem, double x);
};
