#pragma once 
#include <iostream>
#include<fstream>
#include <conio.h>

#include "Method_dual.h"  // алгоритм двойственный
#include "MethodMult_mixed_p.h"
#include "MethodMult_mixed.h"

bool Checked_method(double val_meth, double val_true, double eps)
{
	if (fabs(val_true - val_meth) <= eps)
		return true;
	else
		return false;
}

bool Checked_method_grish(std::vector<double> val_meth, std::vector<double> val_true, double eps)
{
	if (fabs(val_true[0] - val_meth[0]) <= eps && fabs(val_true[1] - val_meth[1]) <= eps)
		return true;
	else
		return false;
}
 
int main()
{
	double E = 0.001, r = 3.2;
	int n = 2, m = 10;
	vector<double> begin(20,0), end(20,0);  // значения начала и конца поискового интервала  для задач Hans, Hill
	int k = 0;
	double count_true = 0;

	std::ofstream out;
	out.open("Graph.txt",  std::ofstream::ios_base::app);
	while (k != 9) {
		std::cout << " 1- hans, 2 - hill, 3 - hans mdual, 4 - hill mdual, 5 - mult_method, 6 - mult_method_p, 7 - mult_method_mixed, 8 - mult_method_mixed_p, 9 - exit " << std::endl;
		std::cin >> k;
		switch (k) {
		case 1:
		{
			count_true = 0;
			out << "Hans" << std::endl;
			for (int index = 0; index < 20; index++)
			{
				
				// инициализируем метод
				Method met(k, index, begin, end, E, r);

				// вызов метода
				met.Solve();

				std::cout << "HansProblem[" << index << "]" << std::endl;
				// вывод истинных значений x,y задачи index
				met.PrintTrueValueHans(index);

				// вывод итерации, на которой было получено лучшее значение
				std::cout << "the best value on " << met.GetBestIndex() << " iterator" << std::endl;

				// вывод полученных  значений x,y задачи index
				std::cout << "x*=" << met.GetOpt().x << "  " << "F(x*)=" << met.GetOpt().z << std::endl;
				std::cout << std::endl;
				if (Checked_method(met.GetOpt().x, met.GetTrueOpt_hans(index), E))
				{
					count_true ++;
					out << met.GetBestIndex() << std::endl;
				}
				else
					out << "wrong" << std::endl;

			}
			out << "count_true " << count_true << std::endl;
			out << "Hans_finish" << std::endl;
			break;
		}
		case 2:
		{
			count_true = 0;
			out << "Hill" << std::endl;
			for (int index = 0; index < 50; index++) {
				// все задачи решаются на итервале [ 0; 1 ]
				std::vector<double> a{0.0};
				std::vector<double> b{ 1.0 };
				Method met(k, index, a, b, E, r);
				met.Solve();
				
				std::cout << "HillProblem[" << index << "]" << std::endl;
				met.PrintTrueValueHill(index);
				std::cout << "the best value on " << met.GetBestIndex() << " iterator" << std::endl;
				std::cout << "x*=" << met.GetOpt().x << "  " << "F(x*)=" << met.GetOpt().z << std::endl;
				std::cout << std::endl;
				if (Checked_method(met.GetOpt().x, met.GetTrueOpt_hill(index), E))
				{
					count_true ++;
					out << met.GetBestIndex() << std::endl;
				}
				else
					out << "wrong" << std::endl;
			}
			out << "count_true " << count_true << std::endl;
			out << "Hill_finish" << std::endl;
			break;
		}
		case 3:
		{
			count_true = 0;
			out << "Hans_dual" << std::endl;
			for (int index = 0; index < 20; index ++)
			{
				MethodDual met(k, index, begin, end, E, r);
				met.SolveDual();
				
				std::cout << "HansProblem[" << index << "]" << std::endl;
				met.PrintTrueValueHans(index);
				std::cout << "the best value on " << met.GetBestIndex() << " iterator" << std::endl;
				std::cout << "x*=" << met.GetOpt().x << "  " << "F(x*)=" << met.GetOpt().z << std::endl;
				std::cout << std::endl;
				if (Checked_method(met.GetOpt().x, met.GetTrueOpt_hans(index), E))
				{
					count_true ++;
					out << met.GetBestIndex() << std::endl;
				}
				else
					out << "wrong" << std::endl;
			}
			out << "count_true " << count_true << std::endl;
			out << "Hans_finish_dual" << std::endl;
			break;
		}
		case 4:
		{
			count_true = 0;
			out << "Hill_dual" << std::endl;
			for (int index = 0; index < 50; index++) {
				std::vector<double> a{ 0.0 };
				std::vector<double> b{ 1.0 };
				MethodDual met(k, index, a, b, E, r);
				met.SolveDual();
				std::cout << "HillProblem[" << index << "]" << std::endl;
				met.PrintTrueValueHill(index);
				std::cout << "the best value on " << met.GetBestIndex() << " iterator" << std::endl;
				std::cout << "x*=" << met.GetOpt().x << "  " << "F(x*)=" << met.GetOpt().z << std::endl;
				std::cout << std::endl;
				if (Checked_method(met.GetOpt().x, met.GetTrueOpt_hill(index), E))
				{
					count_true ++;
					out << met.GetBestIndex() << std::endl;
				}
				else
					out << "wrong" << std::endl;
			}
			out << "count_true " << count_true << std::endl;
			out << "Hill_finish_dual" << std::endl;
			break;

		}
		case 5:
		{
			count_true = 0;
			out << "Mult" << std::endl;
			int index = 0;
			double* y = new double[n];
			//for (int index = 0; index < 100; index++)
			//{
				MethodMult met(index, y, 0, 1, E, r, n, m);
				met.SolveMult(y);
				std::cout << "GrishaginProblem[" << index << "]" << std::endl;
				met.PrintTrueValueGrishagin(index);
				std::cout << "the best value on " << met.GetBestIndex() << " iterator" << std::endl;
				std::cout << "x*= " << met.GetOpt()[0] << " y*= " << met.GetOpt()[1] << "  " << " F(x*)= " << met.GetOpt()[2] << std::endl;
				std::cout << std::endl;
				if (Checked_method_grish(met.GetOpt(), met.GetTrueOpt_grish(index), E))
				{
					count_true ++;
					out << met.GetBestIndex() << std::endl;
				}
				else
					out << "wrong" << std::endl;
			//}
			out << "count_true " << count_true << std::endl;
			out << "Mult_finish" << std::endl;
			break;
		}
		case 6:
		{
			count_true = 0;
			out << "Mult_p" << std::endl;
			double* y = new double[n];
			int index = 0;
			int p = 2;
			//for (int index = 0; index < 100; index++)
			//{
				MethodMult_p met(index, y, 0, 1, E, r, n, m, p);
				met.SolveMult_p(y);
				std::cout << "GrishaginProblem[" << index << "]" << std::endl;
				met.PrintTrueValueGrishagin(index);
				std::cout << "the best value on " << met.GetBestIndex() << " iterator" << std::endl;
				std::cout << "x*= " << met.GetOpt()[0] << " y*= " << met.GetOpt()[1] << "  " << " F(x*)= " << met.GetOpt()[2] << std::endl;
				std::cout << std::endl;
				if (Checked_method_grish(met.GetOpt(), met.GetTrueOpt_grish(index), E))
				{
					count_true++;
					out << met.GetBestIndex() << std::endl;
				}
				else
					out << "wrong" << std::endl;
			//}
			out << "count_true " << count_true << std::endl;
			out << "Mult_p_finish" << std::endl;
			break;
		}
		case 7:
		{
			int _step = 100;
			int _alpha = 15;
			Mixture _mixed(3, 1);
			count_true = 0;
			out << "Mult_mix" << std::endl;
			int index = 0;
			double* y = new double[n];
			for (int index = 0; index < 100; index++)
			{
				MethodMult_mixed met(index, y, 0, 1, E, r, n, m, _step, _mixed, _alpha);
				met.SolveMult_mixed(y);
				std::cout << "GrishaginProblem[" << index << "]" << std::endl;
				met.PrintTrueValueGrishagin(index);
				std::cout << "the best value on " << met.GetBestIndex() << " iterator" << std::endl;
				std::cout << "x*= " << met.GetOpt()[0] << " y*= " << met.GetOpt()[1] << "  " << " F(x*)= " << met.GetOpt()[2] << std::endl;
				std::cout << std::endl;
				if (Checked_method_grish(met.GetOpt(), met.GetTrueOpt_grish(index), E))
				{
					count_true++;
					out << met.GetBestIndex() << std::endl;
				}
				else
					out << "wrong" << std::endl;
			}
			out << "count_true " << count_true << std::endl;
			out << "Mult_mixed_finish" << std::endl;
			break;
		}
		case 8:
		{
			int _step = 50;
			int _alpha = 15;
			int p = 2;
			Mixture _mixed(1, 3);
			count_true = 0;
			out << "Mult_mix_p" << std::endl;
			int index = 0;
			double* y = new double[n];
			for (int index = 0; index < 100; index++)
			{
				MethodMult_mixed_p met(index, y, 0, 1, E, r, n, m, _step, _mixed, _alpha, p);
				met.SolveMult_mixed_p(y);
				std::cout << "GrishaginProblem[" << index << "]" << std::endl;
				met.PrintTrueValueGrishagin(index);
				std::cout << "the best value on " << met.GetBestIndex() << " iterator" << std::endl;
				std::cout << "x*= " << met.GetOpt()[0] << " y*= " << met.GetOpt()[1] << "  " << " F(x*)= " << met.GetOpt()[2] << std::endl;
				std::cout << std::endl;
				if (Checked_method_grish(met.GetOpt(), met.GetTrueOpt_grish(index), E))
				{
					count_true++;
					out << met.GetBestIndex() << std::endl;
				}
				else
					out << "wrong" << std::endl;
			}
			out << "count_true " << count_true << std::endl;
			out << "Mult_mixed_p_finish" << std::endl;
			break;
		}
		}
	}
	out.close();
	_getch();
}
