#include<iostream>
#include<conio.h>
#include"method.h"
#include<fstream>


int main()
{
	double max, E=0.0001, r=3;
	vector<double> begin(20,0), end(20,0);
	std::ofstream out;
	out.open("Graph.txt",  std::ofstream::ios_base::app);
	std::cout << " 1- hans, 2 - hill " << std::endl;
	int k = 0;
	std::cin >> k;
	switch (k) {
	case 1:
	{
		out << "Hans" << std::endl;
		InitInterval(begin, end);
		for (int index = 0; index < 20; index++)
		{
			Method met(k, index, begin[index], end[index], E, r);
			met.solve();
			out << met.GetBestIndex() << std::endl;
			std::cout << "HansProblem[" << index << "]" << std::endl;
			PrintTrueValue(index);
			std::cout << "the best value on " << met.GetBestIndex() << " iterator" << std::endl;
			std::cout << "x*=" << met.GetCurr().x << "  " << "F(x*)=" << met.GetCurr().z << " pos =" << met.GetPos() << std::endl;
			std::cout << std::endl;
		}
		out << "Hans_finish" << std::endl;
		break;
	}
	case 2:
	{
		out << "Hill" << std::endl;
		for (int index = 0; index < 20; index++) {
			Method met(k, index, 0.0, 1.0, E, r);
			met.solve();
			out << met.GetBestIndex() << std::endl;
			std::cout << "HillProblem[" << index << "]" << std::endl;
			PrintTrueValueHill(index);
			std::cout << "the best value on " << met.GetBestIndex() << " iterator" << std::endl;
			std::cout << "x*=" << met.GetCurr().x << "  " << "F(x*)=" << met.GetCurr().z << " pos =" << met.GetPos() << std::endl;
			std::cout << std::endl;
		}
		out << "Hill_finish" << std::endl;
		break;
	}
	}
	out.close();
	_getch();
}



