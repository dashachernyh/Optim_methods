#include<math.h>
#include<iostream>
#include"List.h"

double Funk(double x)
{
	return fabs(x - 1);
}


void rec(TList<double> kor, TList<double> fun,double max,double r,double E, int &Rpos)
{ 
	double curr, R;
	kor.Reset();
	fun.Reset();
	for (int i = 0; i < kor.GetSize(); i++)
	{
		curr = fabs((fun.GetCurr() - fun.GetPrev()) / (kor.GetCurr() - kor.GetPrev()));
		fun.GoNext();
		kor.GoNext();
		if (curr > max)
			max = curr;
		
	}

	if (max == 0) max = 1;
	if (max > 0) max = r * max;

	kor.Reset();
	fun.Reset();
	double Rmax;
	Rmax = max * (kor.GetCurr() - kor.GetPrev()) + (pow((fun.GetCurr() - fun.GetPrev()), 2) / (max*(kor.GetCurr() - kor.GetPrev()))) - 2 * (fun.GetCurr() + fun.GetPrev());
	for (int i = 0; i < kor.GetSize(); i++)
	{
		double k = max * (kor.GetCurr() - kor.GetPrev());
		R = k + (pow((fun.GetCurr() - fun.GetPrev()), 2) / k) - 2 * (fun.GetCurr() + fun.GetPrev());
		if (R > Rmax)
		{
			Rmax = R;
			Rpos =i;
		}
	}

	kor.SetPos(Rpos);
	double Xk;
	fun.SetPos(Rpos);

	Xk = (kor.GetCurr() + kor.GetPrev()) / 2 - (fun.GetCurr() - fun.GetPrev()) / (2 * max);
	kor.InsOrder(Xk);
	fun.SetPos(kor.GetPos());
	if ((kor.GetCurr() - kor.GetPrev()) >= E)
	{
		rec(kor, fun, max, r, E, Rpos);
	}
}

int main()
{
	TList<double> kor, fun;
	double x_0, x_k, max, E, r;
	int *pos;
	pos = 0;

	std::cout << "enter interval X" << std::endl;
	std::cin >> x_0 >> x_k;
	std::cout << "enter E (stop condition)" << std::endl;
	std::cin >> E;
	std::cout << "enter coefficient r" << std::endl;
	std::cin >> r;
	//создаем два списка: один с координатами точек, другой с их значениями (по порядку)

	kor.InsFirst(x_0);
	kor.InsLast(x_k);

	fun.InsFirst(Funk(x_0));
	fun.InsLast(Funk(x_k));

	max = fabs((Funk(x_k) - Funk(x_0)) / (x_k - x_0)); //для самого начала
	rec(kor, fun, max, r, E, *pos);
	kor.SetPos(*pos);
	fun.SetPos(*pos);
	double x_min = kor.GetCurr(), f_min=fun.GetCurr();
	std::cout << "pos =" << pos << "  " << "x*=" << x_min << "  " << "F(x*)=" << f_min << std::endl;
}