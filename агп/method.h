#pragma once
#include <algorithm>
#include<vector>
#include <math.h>

#include"Trial.h"
#include "HansenProblem.hpp"
#include "HansenProblemFamily.hpp"
#include "HillProblem.hpp"
#include "HillProblemFamily.hpp"

// инициализурем поисковый интервал 
void InitInterval(int index, vector<double> &begin, vector<double> &end) 
{
	THansenProblemFamily hansFam;
		hansFam[index]->GetBounds(begin, end);
}

//вывод истинных значений для задачи Hans
void PrintTrueValue(int index_hans)
{
	THansenProblemFamily hansFam;	
	vector<double> optimumPoint;
	optimumPoint = hansFam[index_hans]->GetOptimumPoint();
		std::cout<<"x_min= "<<optimumPoint[0]<<"  F_min= "<< hansFam[index_hans]->GetOptimumValue()<<std::endl;
}

void PrintTrueValueHill(int index_problem) {

	std::cout << "x_min = " << minHill[index_problem][1] << " F_min = " << minHill[index_problem][0] << std::endl;
}


// вычисляет значение функции, key - выбор метода или задачи
double Funk(int key, int index_problem,double x)
{
	if (key == 1 || key == 2) 
		key = 1;
	else 
		key = 2;
	switch (key)
	{
	case 1:
	{
		THansenProblemFamily hansFam;
		return hansFam[index_problem]->ComputeFunction({ x });
	}
	case 2:
	{
		double x_0 = 0.5 + double(index_problem) / 2000.0;
		THillProblemFamily hillFam;
		return hillFam[index_problem]->ComputeFunction({ x });
	}
	}
}

class Method
{
	std::vector<Trial> trials;
	Trial optimum;
	Trial first;
	Trial second;
	int index_problem;
	int best_i;
	double eps, r;
public:
	int key;
	Trial GetOpt() { return optimum; }
	int GetBestIndex() { return best_i; }
	Method(int _key, int index_problem, double x_0, double x_n, double _e, double _r);
	void solve();
};
Method::Method(int _key, int _index_problem,double x_0, double x_n, double _e, double _r)
{
	index_problem = _index_problem;
	key = _key;
	first.x = x_0;
	first.z = Funk(key, index_problem,x_0);
	second.x= x_n;
	second.z = Funk(key, index_problem,x_n);
	trials.push_back(first);
	trials.push_back(second);
	best_i = 0;
	eps = _e;
	r = _r;	
}
void Method:: solve()
{
	Trial current;
	double M,Rmax,Rpos;
	double curr_eps = second.x - first.x;
	std::vector<Trial>::iterator it = trials.begin();
	int itr=0;
	optimum.z = first.z;
	while(curr_eps > eps)
	{
		Rpos = 1;
		M = fabs((trials[1].z - trials[0].z) / (trials[1].x - trials[0].x));
		for (int i = 2; i < trials.size(); i++)
		{
			double max;
			max= fabs((trials.at(i).z - trials.at(i - 1).z) / (trials.at(i).x - trials.at(i - 1).x));
			if (max>M)
				M = max;
		}

		if (M == 0)
			M = 1;
		else
			M = r * M;
		Rmax = M * (trials[1].x-trials[0].x) + (pow((trials[1].z-trials[0].z), 2) / (M * (trials[1].x-trials[0].x))) - 2 * (trials[1].z+trials[0].z);
		for (int i=2; i<trials.size();i++) 
		{
			double k = M * (trials[i].x - trials[i-1].x);
			double R = k + (pow((trials[i].z - trials[i - 1].z), 2) / k) - 2 * (trials[i].z + trials[i - 1].z);	
			if (R > Rmax)
			{
				Rmax = R;
				Rpos = i;
			}
		}
		curr_eps = trials.at(Rpos).x - trials.at(Rpos - 1).x;

		std::vector<Trial>::iterator it2 = trials.begin();
		for (it = trials.begin(); it - trials.begin() <= Rpos; it++) it2 = it;

		current.x = (trials[Rpos].x + trials[Rpos - 1].x) / 2 - (trials[Rpos].z - trials[Rpos - 1].z) / (2 * M);
		current.z = Funk(key, index_problem,current.x);
		trials.insert(it2,current);

		if (optimum.z > current.z)
		{
			best_i = itr;
			optimum = current;
		}
		itr++;
	}
	std::cout << "itr = " << itr << std::endl;
}