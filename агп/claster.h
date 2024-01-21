#pragma once
#include "MethodMult_p.h"
#include "MaxR.h"
#include <chrono>
#include <thread>
#include <random>

struct Characteristic_dual {
	double R, r;
	int pos;
	Trial_thread right;

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

class MethodMult_DualLipsh_p :public MethodMult_p
{
	double r_loc, r_glob;

public:
	MethodMult_DualLipsh_p(int _p) {
		p = _p;
		r_loc = 0;
		r_glob = 0;
	}
	void Init_Dual_p(int _task, int _index_problem, double* y, double _a, double _b, double _e,
		double _r_l, double _r_g, int _n, int _m, const int p) {
		Init_p(_task, _index_problem, y, _a, _b, _e, _r_g, _n, _m, p);
		r_loc = _r_l;
		r_glob = _r_g;
	}
	void SolveMult_DualLipsh_p(double* y, std::vector<std::vector<int>> matrix1,
		std::vector<std::vector<int>> matrix2, std::vector<std::vector<std::vector<int>>> matrix_res,
		int size);
	double Funk_multMat(double* y, std::vector<std::vector<int>> matrix1,
		std::vector<std::vector<int>> matrix2, std::vector<std::vector<std::vector<int>>> matrix_res,
		int size, int index);
};