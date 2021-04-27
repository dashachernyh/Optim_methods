#pragma once
#include <vector>

class MethodMult
{
	std::vector<Trial> trials;        // ������ ���������
	Trial first;                      // ��������� � ������ �����
	Trial second;                     // �� ������ 
	std::vector<double> out_optimal;  // �������� out_optimal[0] = y*, out_optimal[1] = z* 
	int index_problem;
	int best_i;
	double n, m;
	double eps, r;
	double a, b;                       //�������� ����
	double koef[2];                    // ������������ ��������������� koef[0] ��� ������, koef[1] ��� ���������� ��������
public:
	int key;                           
	std::vector<double> GetOpt() { return out_optimal; }  // ���������� ����������� ��������
	int GetBestIndex() { return best_i; }                 
	MethodMult(int _index_problem, double* y, double _a,  double _b, double _e, double _r, double _n, double _m);
	void solve_mult(double * y);
	void ScaleFunc(double y);                             // ������������ ������� ������ (���)
	void InsertScale(double* y);                          // ���������� ��������������� � y
	void PrintTrueValueGrishagin(int index_problem);        
};
