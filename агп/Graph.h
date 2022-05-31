#include "GrishaginOption.hpp"
#include "GrishaginProblemFamily.hpp"
#include "GKLSOption.hpp"
#include "GKLSProblemFamily.hpp"
#include <vector>
#include <fstream>


TGrishaginProblemFamily grish_fam;  //Семейство задач
TGKLSProblemFamily gkls_fam;

void grid_build(std::vector<double>& x, std::vector<double>& y, std::vector<double>& z,
	double a, double b, int index_problem, double h)
{
	std::ofstream out2;
	out2.open("Python_gkls_0.txt", std::ofstream::ios_base::app);

	int s = (b - a) / h;
	for (int i = 0; i < s ; i ++)
	{
		double point = a + i * h;
		x.push_back(point);
		y.push_back(point);
		
	}
	for (int i = 0; i < y.size(); i++)
	{
		for (int j = 0; j < x.size(); j++)
		{
			z.push_back (gkls_fam[index_problem]->ComputeFunction({ x[j], y[i] }));
			// z.push_back(x[j] + y[i]);
			out2 << z[i* x.size() + j]<<" ";
		}
		out2 << std::endl;
	}
	out2.close();
}
