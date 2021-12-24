#pragma once
#include "MethodMult.h"
#include "Mixture.h"


class MethodMult_mixed :public MethodMult
{
	int alpha;               // �������� �����������
	int on_step;             // ��� ��������� �����
	Mixture mix;
public:
	MethodMult_mixed(int _task, int _index_problem, double* y, double _a, double _b, double _e,
		double _r, int _n, int _m, const int _on_step, Mixture _mix, int _alpha);
	void SolveMult_mixed(double* y);
};
