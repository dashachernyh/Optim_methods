#include "GrishaginOption.hpp"
#include "GrishaginProblemFamily.hpp"
#include "GKLSOption.hpp"
#include "GKLSProblemFamily.hpp"
#include "ShekelProblem.hpp"
#include "ShekelProblemFamily.hpp"
#include "HansenProblem.hpp"
#include "HansenProblemFamily.hpp"
#include "HillProblem.hpp"
#include "HillProblemFamily.hpp"
#include <vector>
#include <fstream>

THansenProblemFamily hansFam;
THillProblemFamily hillFam;
TGrishaginProblemFamily grishFam;  //Семейство задач
TGKLSProblemFamily gklsFam;
TShekelProblemFamily shekelFam;

void grid_build(int t_class, std::vector<double>& x, std::vector<double>& z,
	const double a, const double b, int index_problem, const int h)
{
	std::ofstream out2;
	out2.open("Python.txt", std::ofstream::ios_base::app);
	const int s = (b - a) / h;
	double point;
	switch (t_class)
	{
	case(0): {
		// 2-numerics 
		for (int i = 0; i < s; i++)
		{
			point = a + i * h;
			x.push_back(point);

		}
		
		for (int i = 0; i < x.size(); i++)
		{
			for (int j = 0; j < x.size(); j++)
			{
				z.push_back(gklsFam[index_problem]->ComputeFunction({ x[j], x[i] }));
				out2 << z[i * x.size() + j] << " ";
			}
			out2 << std::endl;
		}
		out2.close();

	}
	case(2): {
		// 2-numerics 
		for (int i = 0; i < s; i++)
		{
			point = a + i * h;
			x.push_back(point);
		}

		for (int i = 0; i < x.size(); i++)
		{
			z.push_back(shekelFam[index_problem]->ComputeFunction({ x[i] }));
			out2 << z[i] << " ";
		}
		out2.close();
	}
	}
}

//вычисляет значение функции, key - выбор метода или задачи
double funс(int key, int index_problem, double x)
{
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

vector<double> get_opt(int key, int index_problem) {
	vector<double> optimumPoint;
	
	if (key == 0 || key == 1) {
		optimumPoint = hansFam[index_problem]->GetOptimumPoint();
		optimumPoint.push_back( hansFam[index_problem]->GetOptimumValue());
		return optimumPoint;
	}
	else if (key == 2) // sqrt(module_hans_min)
	{
		optimumPoint = hansFam[index_problem]->GetMaxPoint();
		optimumPoint.push_back(hansFam[index_problem]->GetMaxValue());
		return optimumPoint;
	}
	else if (key ==3 || key == 4) // module(_h_max)
	{
		optimumPoint = hillFam[index_problem]->GetOptimumPoint();
		optimumPoint.push_back(hillFam[index_problem]->GetOptimumValue());
		return optimumPoint;
	}
	else if (key == 5) // module_hill
	{
		optimumPoint = hillFam[index_problem]->GetMaxPoint();
		optimumPoint.push_back(hillFam[index_problem]->GetMaxValue());
		return optimumPoint;
	}
	else if (key ==6 || key ==7) //sqrt(module_hill_min)
	{
		optimumPoint = shekelFam[index_problem]->GetOptimumPoint();
		optimumPoint.push_back(shekelFam[index_problem]->GetOptimumValue());
		return optimumPoint;
	}
	else // module(_hill_max)
	{
		optimumPoint = shekelFam[index_problem]->GetMaxPoint();
		optimumPoint.push_back(shekelFam[index_problem]->GetMaxValue());
		return optimumPoint;
	}
}

void grid_build_txt(int t_class, const double a, const double b, int index_problem, const int s,
	int key)
{
	std::ofstream out2;
	out2.open("Python_txt.txt", std::ofstream::ios_base::app);
	double  h = (b - a) / s;
	std::vector<double> x(s);
	double point;
	switch (t_class)
	{
	case(0): {

		// 2-numerics 
		for (int i = 0; i < s; i++)
		{
			point = a + i * h;
			x[i]  =point;

		}

		for (int i = 0; i < x.size(); i++)
		{
			for (int j = 0; j < x.size(); j++)
			{
				out2 << gklsFam[index_problem]->ComputeFunction({ x[j], x[i] }) << " ";
				//out2 << grish_fam[index_problem]->ComputeFunction({ x[j], x[i] }) << " ";
			}
			out2 << std::endl;
		}
		out2.close();

	}
	case(1): {
		vector<double> opt_val = get_opt(key, index_problem);
		out2 << opt_val[0] << " "<< opt_val[1]<<"\n";
		for (int i = 0; i < s; i++)
		{
			point = a + i * h;
			x[i] = point;
			out2 << point << " ";
		}
		out2 << "\n";
		for (int i = 0; i < x.size(); i++)
		{
			out2 << funс(key, index_problem, x[i]) << " ";
		}
		out2 << "\n";
		out2.close();
	}
	}
}
