#pragma once
#include "MethodMult.h"

#include <omp.h>

struct Trial_thread
{
	double x, z;
	int thread;
	Trial_thread()
	{
		x = 0;
		z = 0;
		thread = 0;
	}
	Trial_thread& operator = (const Trial_thread& tr)
	{
		x = tr.x;
		z = tr.z;
		thread = tr.thread;
		return *this;
	}
	bool operator == (const Trial_thread& tr)
	{
		if (x == tr.x && z == tr.z) return true;
		else
			return false;
	}
	bool operator < (const Trial_thread& tr) // ��� ������� sort
	{
		if (x < tr.x)return true;
		else
			return false;
	}
};

struct Characteristic {
	double R;
	int pos;

	Characteristic(double _R = 0, int _pos = 0) {
		R = _R;
		pos = _pos;
	}
	bool operator <  (const Characteristic& ch)
	{
		return R < ch.R;
	}

};

class MethodMult_p :public MethodMult
{
protected:
	std::vector<Trial_thread> trials_thread;
	int p;                       // �������� ��� ������������ ����������, ���-�� �������
public:
	MethodMult_p()
	{
		p = 0;
	}
	void Init_p(int _task, int _index_problem, double* y, double _a, double _b, double _e,
		double _r, int _n, int _m, const int p);
	void ClearMethod_p();
	void SolveMult_p(double* y);
};
