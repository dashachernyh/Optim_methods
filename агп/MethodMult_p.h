#pragma once
#include "MethodMult.h"

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
	bool operator < (const Trial_thread& tr) // для функции sort
	{
		if (x < tr.x)return true;
		else
			return false;
	}
};

struct Characteristic {
	double r;
	int pos;

	Characteristic(double _r = 0, int _pos = 0) {
		r = _r;
		pos = _pos;
	}
	bool operator <  (const Characteristic& ch2)
	{
		return r < ch2.r;
	}

};

class MethodMult_p :public MethodMult
{
protected:
	std::vector<Trial_thread> trials_thread;
	int p;                       // параметр при параллельной реализации, кол-во потоков
public:
	MethodMult_p()
	{
		p = 0;

	}
	MethodMult_p(int _index_problem, double* y, double _a, double _b, double _e,
		double _r, int _n, int _m, const int p);
	void SolveMult_p(double* y);
};
