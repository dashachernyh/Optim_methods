#include <iostream>
#include <conio.h>
#include <fstream>

#include "Method_Strongina.h" // алгоритм Стронгина
#include "method_dual.h"  // алгоритм двойственный
#include "Map.h"
#include "Method_mult.h"  // алгоритм для ф-ии нескольких переменных

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
	double E = 0.0001, r = 2.5;
	int n = 2, m = 10;
	vector<double> begin(20,0), end(20,0);  // значения начала и конца поискового интервала  для задач Hans, Hill
	int k = 0;
	double count_true = 0;

	std::ofstream out;
	out.open("Graph.txt",  std::ofstream::ios_base::app);
	while (k != 6) {
		std::cout << " 1- hans, 2 - hill, 3 - hans mdual, 4 - hill mdual, 5 - mult_method, 6 - exit " << std::endl;
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
				met.solve();

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
				met.solve();
				
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
				met.solve_dual();
				
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
				met.solve_dual();
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
			double* y = new double[n];
			int key = 1;
			int index = 1;
			for (int index = 0; index < 50; index++)
			{
				MethodMult met(index, y, 0, 1, E, r, n, m);
				met.solve_mult(y);
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
			}
			out << "count_true " << count_true << std::endl;
			out << "Mult_finish" << std::endl;
			break;
		}
		}
	}
	out.close();
	_getch();
}
