#pragma once
#include "MethodMult.h"
#include "MaxR.h"

class MethodMult_DualLipsh :public MethodMult
{
protected:
	double r_loc, r_glob;
public:
	MethodMult_DualLipsh(int _task, int _index_problem, double* y, double _a, double _b,
		double _e, double r_l, double r_g, int _n, int _m)
		:MethodMult(_task,_index_problem, y, _a, _b, _e, r_g, _n, _m), r_loc(r_l), r_glob(r_g) {}

	void SolveMult_DualLipsh(double* y);
};