#include<math.h>
#include<iostream>
#include<vector>
#include <algorithm>
#include<conio.h>
#include"method.h"

using namespace std;

int main()
{
	
	double x_0 , x_k , max, E=0.00001, r=3.01;

	std::cout << "enter interval X" << std::endl;
	std::cin >> x_0 >> x_k;
	/*std::cout << "enter E (stop condition)" << std::endl;
	std::cin >> E;
	std::cout << "enter coefficient r" << std::endl;
	std::cin >> r;*/
	
	Method met(x_0,x_k);
	met.solve(r, E);
	
	std::cout << "pos =" << met.GetPos() << "  " << "x*=" << met.GetCurr().x << "  " << "F(x*)=" <<met.GetCurr().z << std::endl;
	_getch();
}



