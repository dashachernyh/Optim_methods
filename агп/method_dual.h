#pragma once
#include "MethodStrongina.h"

class MethodDual:public Method
{
public:
	MethodDual(int _key, int index_problem, std::vector<double> x_0, std::vector<double> x_n, double _e, double _r);
	void SolveDual();
};
