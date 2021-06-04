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

	Method() {
		Trial n;
		n.x = 0;
		n.z = 0;
		trials.push_back(n);
		first = second= optimum = n;
		index_problem = best_i = eps = r = 0;
	};

	Method(int _key, int index_problem, std::vector<double> x_0, std::vector<double> x_n, double _e, double _r);
	Trial GetOpt() { return optimum; }
	int GetBestIndex() { return best_i; }
	void InitIntervalHans(int index, vector<double> &begin, vector<double> &end);
	void solve();
	void  PrintTrueValueHans(int index_problem);
	void PrintTrueValueHill(int index_problem);
	double GetTrueOpt_hans(int index_problem);
	double GetTrueOpt_hill(int index_problem);
protected:
	double Funk(int key, int index_problem, double x);
};
