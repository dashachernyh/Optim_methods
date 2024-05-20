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
TGKLSProblemFamily gklsFam{ 2, GKLSClass::Hard, GKLSFuncionType::TD };
TShekelProblemFamily shekelFam;


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

double funс_mult(int key, int index_problem, std::vector<double> y)
{
	switch (key) {
		case 0: {
			return grishFam[index_problem]->ComputeFunction(y);
		}
		case 1: {
			return fabs(grishFam[index_problem]->ComputeFunction(y) - grishFam[index_problem]->GetOptimumValue());
		}
		case 2: {
			return sqrt(fabs(grishFam[index_problem]->ComputeFunction(y) - grishFam[index_problem]->GetOptimumValue()));
		}
		case 3: {

			return fabs(-(grishFam[index_problem]->ComputeFunction(y) - grishFam[index_problem]->GetMaxValue()));
		}
		case 4: {
			return gklsFam[index_problem]->ComputeFunction(y);
		}
		case 5: {
			return fabs(gklsFam[index_problem]->ComputeFunction(y) - gklsFam[index_problem]->GetOptimumValue());
		}
		case 6: {
			return sqrt(fabs(gklsFam[index_problem]->ComputeFunction(y) - gklsFam[index_problem]->GetOptimumValue()));
		}
		case 7: {
			return fabs(-gklsFam[index_problem]->ComputeFunction(y) + gklsFam[index_problem]->GetMaxValue());
		}
	}
}

vector<double> get_opt(int key, int index_problem) {
	vector<double> optimumPoint;

	if (key == 0 || key == 1) {
		optimumPoint = hansFam[index_problem]->GetOptimumPoint();
		optimumPoint.push_back(hansFam[index_problem]->GetOptimumValue());
		return optimumPoint;
	}
	else if (key == 2) // sqrt(module_hans_min)
	{
		optimumPoint = hansFam[index_problem]->GetMaxPoint();
		optimumPoint.push_back(hansFam[index_problem]->GetMaxValue());
		return optimumPoint;
	}
	else if (key == 3 || key == 4) // module(_h_max)
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
	else if (key == 6 || key == 7) //sqrt(module_hill_min)
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

vector<double> get_opt_mult(int key, int index_problem) {
	vector<double> optimumPoint;

	if (key >= 0 and key <= 2) {
		optimumPoint = grishFam[index_problem]->GetOptimumPoint();
		//optimumPoint.push_back(grishFam[index_problem]->GetOptimumValue());
		
	}
	else if (key == 3) 
	{
		optimumPoint = hansFam[index_problem]->GetMaxPoint();
		//optimumPoint.push_back(hansFam[index_problem]->GetMaxValue());
	}
	else if (key >= 4 and key <= 6) // module(_h_max)
	{
		optimumPoint = gklsFam[index_problem]->GetOptimumPoint();
		//optimumPoint.push_back(hillFam[index_problem]->GetOptimumValue());
	}
	else if (key == 7) // module_hill
	{
		optimumPoint = gklsFam[index_problem]->GetMaxPoint();
		//optimumPoint.push_back(hillFam[index_problem]->GetMaxValue());
	}
	return optimumPoint;
}

void grid_build(int t_class, int key,
	const double a, const double b, int index_problem, const int n)
{
	std::ofstream out2;
	out2.open("Python50.txt", std::ofstream::ios_base::app);
	const double h = double(b - a) / double(n);
	double point;
	std::vector<double> x(n+1);
	out2 << index_problem << "\n";
	switch (t_class)
	{
	case(0): {
		vector<double> opt_val = get_opt(key, index_problem);
		out2 << opt_val[0] << " " << opt_val[1] << "\n";
		for (int i = 0; i < n + 1; i++)
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
		break;
	}
	case(1): {
		// 2-numerics 
		vector<double> opt_val = get_opt_mult(key, index_problem);
		out2 << opt_val[0] << " " << opt_val[1] << "\n";
		for (int i = 0; i < n; i++)
		{
			x[i] = a + i * h;

		}

		for (int i = 0; i < x.size(); i++)
		{
			for (int j = 0; j < x.size(); j++)
			{
				//out2 << x[i] << " " << x[j] << " ";
				out2 << funс_mult(key, index_problem, { x[j], x[i] }) << " ";
			}
			out2 << "\n";
		}
		out2.close();
		break;
	}
	}
}
