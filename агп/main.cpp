#include <iostream>
#include <conio.h>
#include <fstream>

#include "Method_Strongina.h" // �������� ���������
#include "method_dual.h"  // �������� ������������
#include "Map.h"
#include "Method_mult.h"  // �������� ��� �-�� ���������� ����������


int main()
{
	double E = 0.001, r = 3;
	int n = 2, m = 10;
	vector<double> begin(20,0), end(20,0);  // �������� ������ � ����� ���������� ���������  ��� ����� Hans, Hill

	int k = 0;
	std::ofstream out;
	out.open("Graph.txt",  std::ofstream::ios_base::app);
	while (k != 6) {
		std::cout << " 1- hans, 2 - hill, 3 - hans mdual, 4 - hill mdual, 5 - mult_method, 6 - exit " << std::endl;
		std::cin >> k;
		switch (k) {
		case 1:
		{
			out << "Hans" << std::endl;
			for (int index = 0; index < 20; index++)
			{
				// �������������� �����
				Method met(k, index, begin[0], end[0], E, r);

				// �������������� ��������� ��������� ��� ������ � ������� index
				met.InitIntervalHans(index, begin, end);

				// ����� ������
				met.solve();

				out << met.GetBestIndex() << std::endl;
				std::cout << "HansProblem[" << index << "]" << std::endl;
				// ����� �������� �������� x,y ������ index
				met.PrintTrueValueHans(index);

				// ����� ��������, �� ������� ���� �������� ������ ��������
				std::cout << "the best value on " << met.GetBestIndex() << " iterator" << std::endl;

				// ����� ����������  �������� x,y ������ index
				std::cout << "x*=" << met.GetOpt().x << "  " << "F(x*)=" << met.GetOpt().z << std::endl;
				std::cout << std::endl;
			}
			out << "Hans_finish" << std::endl;
			break;
		}
		case 2:
		{
			out << "Hill" << std::endl;
			for (int index = 0; index < 20; index++) {
				// ��� ������ �������� �� �������� [ 0; 1 ]
				Method met(k, index, 0.0, 1.0, E, r);
				met.solve();
				out << met.GetBestIndex() << std::endl;
				std::cout << "HillProblem[" << index << "]" << std::endl;
				met.PrintTrueValueHill(index);
				std::cout << "the best value on " << met.GetBestIndex() << " iterator" << std::endl;
				std::cout << "x*=" << met.GetOpt().x << "  " << "F(x*)=" << met.GetOpt().z << std::endl;
				std::cout << std::endl;
			}
			out << "Hill_finish" << std::endl;
			break;
		}
		case 3:
		{
			out << "Hans_dual" << std::endl;
			for (int index = 0; index < 20; index++)
			{
				MethodDual met(k, index, begin[0], end[0], E, r);
				met.InitIntervalHans(index, begin, end);
				met.solve_dual();
				out << met.GetBestIndex() << std::endl;
				std::cout << "HansProblem[" << index << "]" << std::endl;
				met.PrintTrueValueHans(index);
				std::cout << "the best value on " << met.GetBestIndex() << " iterator" << std::endl;
				std::cout << "x*=" << met.GetOpt().x << "  " << "F(x*)=" << met.GetOpt().z << std::endl;
				std::cout << std::endl;
			}
			out << "Hans_finish_dual" << std::endl;
			break;
		}
		case 4:
		{
			out << "Hill_dual" << std::endl;
			for (int index = 0; index < 20; index++) {
				MethodDual met(k, index, 0.0, 1.0, E, r);
				met.solve_dual();
				out << met.GetBestIndex() << std::endl;
				std::cout << "HillProblem[" << index << "]" << std::endl;
				met.PrintTrueValueHill(index);
				std::cout << "the best value on " << met.GetBestIndex() << " iterator" << std::endl;
				std::cout << "x*=" << met.GetOpt().x << "  " << "F(x*)=" << met.GetOpt().z << std::endl;
				std::cout << std::endl;
			}
			out << "Hill_finish_dual" << std::endl;
			break;
		}
		case 5:
		{
			double* y = new double[n];
			int key = 1;
			int index = 1;
			for (int index = 0; index < 10; index++)
			{
				MethodMult met(index, y, 0, 1, E, r, n, m);
				met.solve_mult(y);
				std::cout << "GrishaginProblem[" << index << "]" << std::endl;
				met.PrintTrueValueGrishagin(index);
				std::cout << "the best value on " << met.GetBestIndex() << " iterator" << std::endl;
				std::cout << "x*= " << met.GetOpt()[0] << " y*= " << met.GetOpt()[1] << "  " << " F(x*)= " << met.GetOpt()[2] << std::endl;
				std::cout << std::endl;
			}
		}
		}
	}
	out.close();
	_getch();
}
