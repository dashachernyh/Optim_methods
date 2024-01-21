#pragma once
#include "MethodMult.h"
#include "Mixture.h"


class MethodMult_mixed :public MethodMult
{
	int alpha;               // параметр локализации
	int on_step;             // шаг включения смеси
	Mixture mix;
public:
	MethodMult_mixed() {
		alpha = on_step = -1;
	}
	void Init_Mix(int _task, int _index_problem, double* y, double _a, double _b, double _e,
		double _r, int _n, int _m, int _on_step, Mixture _mix, int _alpha);
	void SolveMult_mixed(double* y);
};
