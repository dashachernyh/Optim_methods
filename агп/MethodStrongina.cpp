#include"MethodStrongina.h"
#include<iostream>
#include<fstream>

// инициализурем поисковый интервал 
void Method::InitIntervalHans(int index, vector<double>& begin, vector<double>& end)
{
	hansFam[index]->GetBounds(begin, end);
}

double Method::GetTrueOpt(int index_problem, int task)
{
	if (task == 0 || task == 1) {
		return hansFam[index_problem]->GetOptimumPoint()[0];
	}
	else if (task == 2) {
		return hansFam[index_problem]->GetMaxPoint()[0];
	}
	else if (task == 3 || task == 4 ) {
		return hillFam[index_problem]->GetOptimumPoint()[0];
	}
	else if (task == 5) {
		return hillFam[index_problem]->GetMaxPoint()[0];
	}
	else if (task == 6 || task == 7) {
		return shekelFam[index_problem]->GetOptimumPoint()[0];
	}
	else {
		return shekelFam[index_problem]->GetMaxPoint()[0];
	}
}

//вывод истинных значений для задачи Hans
void  Method::PrintTrueValue(int index_problem, int task)
{
	if (task >= 0 and task <=2) {
		vector<double> optimumPoint;
		optimumPoint = hansFam[index_problem]->GetOptimumPoint();
		std::cout << "x_min= " << optimumPoint[0] << "  F_min= " << hansFam[index_problem]->GetOptimumValue() << std::endl;
	}
	else if (task >= 3 and task <= 5){
		std::cout << "x_min = " << minHill[index_problem][1] << " F_min = " << minHill[index_problem][0] << std::endl;
		std::cout << "x_max= " << hillFam[index_problem]->GetMaxPoint()[0] << " F_max = " << hillFam[index_problem]->GetMaxValue() << std::endl;
	}
	else {
		std::cout << "x_min= " << shekelFam[index_problem]->GetOptimumPoint()[0] << "  F_min= " << shekelFam[index_problem]->GetOptimumValue() << std::endl;
		std::cout << "x_max= " << shekelFam[index_problem]->GetMaxPoint()[0] << " F_max = " << shekelFam[index_problem]->GetMaxValue() << std::endl;
	}
}

// вычисляет значение функции, key - выбор метода или задачи
double Method::Funс(int method, int key, int index_problem, double x)
{
	if (method) {
		switch (key) {
		case(0): // module_hans
		{
			return fabs(hansFam[index_problem]->ComputeFunction({ x }) - hansFam[index_problem]->GetOptimumValue());
		}
		case(1): // sqrt(module_hans_min)
		{
			return sqrt(fabs(hansFam[index_problem]->ComputeFunction({ x }) - hansFam[index_problem]->GetOptimumValue()));
		}
		case(2): // module(_h_max)
		{
			return fabs(-1 * (hansFam[index_problem]->ComputeFunction({ x }) - hansFam[index_problem]->GetMaxValue()));
		}
		case(3): // module_hill
		{
			return fabs(hillFam[index_problem]->ComputeFunction({ x }) - hillFam[index_problem]->GetOptimumValue());
		}
		case(4): //sqrt(module_hill_min)
		{
			return sqrt(fabs(hillFam[index_problem]->ComputeFunction({ x }) - hillFam[index_problem]->GetOptimumValue()));
		}
		case(5): // module(_hill_max)
		{
			return fabs(-1 * (hillFam[index_problem]->ComputeFunction({ x }) - hillFam[index_problem]->GetMaxValue()));
		}
		case(6):
		{
			return fabs(shekelFam[index_problem]->ComputeFunction({ x }) - shekelFam[index_problem]->GetOptimumValue());
		}
		case(7):
		{
			return sqrt(fabs(shekelFam[index_problem]->ComputeFunction({ x }) - shekelFam[index_problem]->GetOptimumValue()));
		}
		case(8):
		{
			return fabs(-1 * (shekelFam[index_problem]->ComputeFunction({ x }) - shekelFam[index_problem]->GetMaxValue()));
		}
		}
	}
	else {
		switch (key) {
		case(0):
			return hansFam[index_problem]->ComputeFunction({ x });
		case(1):
			return hillFam[index_problem]->ComputeFunction({ x });
		case(2):
			return shekelFam[index_problem]->ComputeFunction({ x });
		}
	}
}

void Method::Init(int method, int _key, int _index_problem, std::vector<double> x_0,
	std::vector<double> x_n, double _e, double _r)
{
	std::ofstream f_out;
	f_out.open("graph265.txt", std::ofstream::ios_base::app);
	f_out << GetTrueOpt(_index_problem, _key) <<" " << Funс(method, _key, _index_problem, GetTrueOpt(_index_problem, _key))<< "\n";
	Trial first, second;
	index_problem = _index_problem;
	key = _key;
	if (key == 0)
		InitIntervalHans(index_problem, x_0, x_n);
	first.x = x_0[0];
	first.z = Funс(method, key, index_problem, first.x);
	second.x = x_n[0];
	second.z = Funс(method, key, index_problem, second.x);
	trials.push_back(first);
	trials.push_back(second);
	best_i = 0;
	eps = _e;
	r = _r;

	optimum = first.z < second.z ? first : second;
	f_out << first.x << " " << first.z << "\n";
	f_out << second.x << " " << second.z << "\n";
	f_out.close();
	//zpos = 1;
}

void Method::Solve()
{
	Trial current;
	double M, Rmax, max, k, R;
	size_t Rpos;
	double curr_eps = trials[1].x - trials[0].x;
	int itr = 0;
	std::ofstream f_out;
	f_out.open("graph265.txt", std::ofstream::ios_base::app);
	while (curr_eps > eps && itr <= 10000) //(curr_eps > eps) optimum.z
	{
		Rpos = 1;
		M = fabs((trials[1].z - trials[0].z) / (trials[1].x - trials[0].x));
		for (size_t i = 2; i < trials.size(); i++)
		{
			max = fabs((trials[i].z - trials[i - size_t(1)].z)
				/ (trials[i].x - trials[i - size_t(1)].x));
			if (max > M)
				M = max;
		}
		M = (M == 0) ? 1 : r * M;

		Rmax = M * (trials[1].x - trials[0].x) + (pow((trials[1].z - trials[0].z), 2)
			/ (M * (trials[1].x - trials[0].x))) - 2 * (trials[1].z + trials[0].z);
		for (int i = 2; i < trials.size(); i++)
		{
			k = M * (trials[i].x - trials[i - size_t(1)].x);
			R = k + (pow((trials[i].z - trials[i - size_t(1)].z), 2) / k)
				- 2 * (trials[i].z + trials[i - size_t(1)].z);
			if (R > Rmax)
			{
				Rmax = R;
				Rpos = i;
			}
		}
		curr_eps = trials[Rpos].x - trials[Rpos - size_t(1)].x;

		std::vector<Trial>::iterator it2 = trials.begin();
		for (auto it = trials.begin(); it - trials.begin() <= Rpos; it++) it2 = it;

		current.x = (trials[Rpos].x + trials[Rpos - 1].x) / 2 - (trials[Rpos].z
			- trials[Rpos - size_t(1)].z) / (2 * M);
		current.z = Funс(1, key, index_problem, current.x);
		trials.insert(it2, current);
		f_out << current.x << " " << current.z << "\n";
		if (current.x == 0.750128) {
			double mp = 0;
		}
		if (optimum.z > current.z)
		{
			best_i = itr;
			optimum = current;
		}
		itr++;
	}
	std::cout << "itr = " << itr << std::endl;
	f_out.close();
}