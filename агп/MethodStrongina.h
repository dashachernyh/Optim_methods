#pragma once
#include <algorithm>
#include<vector>
#include <math.h>

#include"Trial.h"
#include "HansenProblem.hpp"
#include "HansenProblemFamily.hpp"
#include "HillProblem.hpp"
#include "HillProblemFamily.hpp"
#include "ShekelProblem.hpp"
#include "ShekelProblemFamily.hpp"

class Method
{
protected:
	std::vector<Trial> trials;
	Trial optimum;
	int index_problem;
	int best_i;
	int key;
	int check_method;
	double eps, r;
	THillProblemFamily hillFam;
	THansenProblemFamily hansFam;
	TShekelProblemFamily shekelFam;
	//double zpos;

public:
	Method() {
		Trial zero;
		zero.x = 0;
		zero.z = 0;
		optimum = zero;
		key = 0;
		index_problem = best_i= -1;
		r = eps = 0;
	};
	void Clear() {
		trials.clear();
		best_i = -1;
		optimum.x = NULL;
		optimum.z = NULL;
	}
	//Method(int _key, int index_problem, std::vector<double> x_0, std::vector<double> x_n, double _e, double _r);
	void Init(int method, int _check_method, int _key, int index_problem, std::vector<double> x_0,
		std::vector<double> x_n, double _e, double _r = 0);
	Trial GetOpt() { return optimum; }
	int GetBestIndex() { return best_i; }
	void InitIntervalHans(int index, vector<double> &begin, vector<double> &end);
	void Solve();
	void  PrintTrueValue(int index_problem, int task);
	double GetTrueOpt(int index_problem, int task);
	double Funñ (int method, int key, int index_problem, double x);
};
