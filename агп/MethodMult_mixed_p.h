#pragma once
#include "MethodMult_p.h"
#include "Mixture.h"

class MethodMult_mixed_p :public MethodMult_p
{
	int alpha;               // параметр локализации
	int on_step;             // шаг включения смеси
	Mixture mix;
public:
	MethodMult_mixed_p(int _task, int _index_problem, double* y, double _a, double _b, double _e, double _r,
		int _n, int _m, const int _on_step, Mixture _mix, int _alpha, int _p);
	void SolveMult_mixed_p(double* y);
};

