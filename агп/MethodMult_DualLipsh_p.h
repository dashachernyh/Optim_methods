#pragma once
#include "MethodMult_p.h"
#include "MaxR.h"


struct Characteristic_dual {
	double R, r;
	int pos;

	Characteristic_dual(double _R = 0, double _r = 0, int _pos = 0) {
		R = _R;
		r = _r;
		pos = _pos;
	}
	bool operator <  (const Characteristic_dual& ch)
	{
		return R < ch.R;
	}

};

class MethodMult_DualLipsh_p:public MethodMult_p
{
	double r_loc, r_glob;

public:
	MethodMult_DualLipsh_p(int _task, int _index_problem, double* y, double _a, double _b, double _e,
		 double _r_l, double _r_g, int _n, int _m, const int p)
		:MethodMult_p(_task,_index_problem, y, _a,_b, _e, _r_g, _n, _m, p),
		r_loc(_r_l), r_glob(_r_g){}
	void SolveMult_DualLipsh_p(double* y);
};