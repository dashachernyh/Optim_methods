#pragma once 
//#include <iostream>
#include<fstream>
#include <conio.h>
#include <iomanip>
#include <chrono>

#include "Method_dual.h"  // include MethodStrongina
#include "MethodMult_mixed.h"  // include MethodMult
#include "MethodMult_mixed_p.h"  // include MethodMult_p // include MethodMult
#include "MethodMult_DualLipsh.h"  // include MethodMult
#include "MethodMult_DualLipsh_p.h" // include MethodMult_p
#include "Method_SearchRoot.hpp"
#include "Method_AGP_ASR.hpp"
#include "MethodMult_AGP_ASR.hpp"
#include "MethodMult_SearchRoot.hpp"
#include "Graph.h"



bool Checked_method(double val_meth, double val_true, double eps)
{
	if (fabs(val_true - val_meth) <= eps)
		return true;
	else
		return false;
}

bool Checked_method_mult(int n, std::vector<double> val_meth, std::vector<double> val_true, double eps)
{
	if (n == 2) {
		if (fabs(val_true[0] - val_meth[0]) <= eps && fabs(val_true[1] - val_meth[1]) <= eps)
			return true;
		else
			return false;
	}
	else {
		if (fabs(val_true[0] - val_meth[0]) <= eps)
			return true;
		else
			return false;
	}
}

bool Checked_value(double val_meth, double val_true, double eps)
{
	if (fabs(val_true - val_meth) <= eps)
		return true;
	else
		return false;
}
 
int main()
{
	double E = 0.005, r =6.0;
	int n = 1, m = 10;
	int key = 0;
	double count_true = 0;
	double average_time = 0;
	double aver = 0;
	std::ofstream out;
	out.open("Graph.txt",  std::ofstream::ios_base::app);
	while (key != 14) {
		std::cout << " 0- grid, 1- mstrong, 2- mdual, 3 - mult_method, 4 - mult_method_p,";
		std::cout << "\n5 - mult_method_mixed, 6 - mult_method_mixed_p, 7 - mult_metgod_dualLipsh, 9- mult_metgod_dualLipsh_p,";
		std::cout << "\n10 - searchRoot, 11 - AGP_ASR, 12 - mult AGP_ASR, 13 - mult_ASR" << std::endl;
		std::cin >> key;
		/*std::cout << "Enter parameters eps, r "<<std::endl;
		std::cin >> E >> r;*/
		switch (key) {
		case(0): {

			int task_dim = 1;
			int index_task = 0;
			int len = 1; //  1 - HILL, 10 - Shekel, 2- GKLS, 1- Grish
			std::cout << "0- Hill, Shekel, 1 - Gkls, Grish\n";
			std::cin >> task_dim;
			std::cout << "[0-2]-Hans, [3-5]-Hill, [6-8]-Shekel";
			std::cout << "[0-3] - Grish, [4-7] - GKLS\n";
			std::cin >> index_task;

			count_true = 0;
			// значения начала и конца поискового интервала  для задач Hans
			double a, b;
			int index_problem = 50;
			if (task_dim == 0 and index_task >= 3 and index_task <= 5) {
				// Hill на [0; 1]
				a = 0.0;
				b = 1.0;
			}
			else if (task_dim == 0 and index_task >= 6 and index_task <= 8) {
				// Shekel на [0; 10]
				len = 10;
				a = 0.0;
				b = 10.0;
			}
			else if (task_dim == 1 and index_task >= 0 and index_task <= 3) {
				len = 1;
				a = 0.0;
				b = 1.0;
			}
			else {
				len = 2;
				a = -1.0;
				b = 1.0;
			}
			grid_build(task_dim, index_task, a, b, index_problem, 1000);
			break;
		}
		case 1:
		{
			int task = 0;
			int len = 1; //  1 - Grishagin, 2 - GKLS, 1 - HILL, 10 - Shekel
			std::cout << "[0-2] - Hans, [3-5] - Hill, [6-8] - Shekel\n";
			std::cin >> task;
			int check_method = 0;
			std::cout << "check_method: 0- z_min, 1 - curr_eps\n";
			std::cin >> check_method;
			count_true = 0;
			// значения начала и конца поискового интервала  для задач Hans
			vector<double> begin(20, 0), end(20, 0);
			int maxTask = 20;
			if (task >= 3 and task <= 5) {
				// Hill на [0; 1]
				begin = { 0.0 };
				end = { 1.0 };
				maxTask = 1000; // 1000
			}
			else if (task >= 6 and task <= 8) {
				// Shekel на [0; 10]
				len = 10;
				begin = { 0.0 };
				end = { 10.0 };
				maxTask = 1000;
				
			}
			Method met;
			aver = 0;
			out <<"["<<key<< "] Task" << task << " e " << len * E<< " r " << r<< " check "<< check_method << std::endl;
			for (int index = 0; index < maxTask; index++)
			{
				// инициализируем метод 0 - usual, 1 - fabs
				met.Init(1, check_method,  task, index, begin, end, len * E, r);
				// вызов метода
				met.Solve();

				std::cout << "Problem[" << index << "]" << std::endl;
				// вывод истинных значений x,y задачи index
				met.PrintTrueValue(index, task);

				// вывод итерации, на которой было получено лучшее значение
				std::cout << "the best value on " << met.GetBestIndex() << " iterator" << std::endl;
				// вывод полученных  значений x,y задачи index
				std::cout << "x*=" << met.GetOpt().x << "  " << "F(x*)=" << met.GetOpt().z << std::endl;
				std::cout << std::endl;
				bool check = (check_method == 0) ? met.GetOpt().z <= len * E : Checked_method(met.GetOpt().x, met.GetTrueOpt(index, task), len * E);
				if (check) //)met.GetOpt().z <= len * E Checked_method(met.GetOpt().x, met.GetTrueOpt(index, task), E)
				{
					count_true++;
					aver += met.GetBestIndex();
					out << met.GetBestIndex() << std::endl;
				}
				else {
					//Checked_method(met.GetOpt().x, met.GetTrueOpt(index, task), E);
					out << "wrong " << index << std::endl;
				}
				met.Clear();
			}
			out << "count_true " << count_true << " aver= " << aver / count_true << std::endl;
			std::cout << "count_true " << count_true << " aver= " << aver / count_true << std::endl;
			out << "Task_finish" << std::endl;
			break;
		}
		case 2:
		{
			int task = 0;    // 0 - Hans, 1 - Hill
			int len = 1; //  1 - Grishagin, 2 - GKLS, 1 - HILL, 10 - Shekel
			std::cout << "[0-2] - Hans, [3-5] - Hill, [6-8] - Shekel\n";
			std::cin >> task;

			int check_method = 0;
			std::cout << "check_method: 0- z_min, 1 - curr_eps\n";
			std::cin >> check_method;

			count_true = 0;
			// значения начала и конца поискового интервала  для задач Hans
			vector<double> begin(20, 0), end(20, 0);
			int maxTask = 20;
			if (task == 1) {
				// Hill на [0; 1]
				begin = { 0.0 };
				end = { 1.0 };
				maxTask = 1000; // 1000
			}
			else if (task == 2) {
				// Shekel на [0; 10]
				len = 10;
				begin = { 0.0 };
				end = { 10.0 };
				maxTask = 1000;
			}
			MethodDual met;
			out << "[" << key << "] Task" << task << " e " << len * E << " r " << r << " check " << check_method << std::endl;
			for (int index = 0; index < maxTask; index++) {

				met.Init(0, check_method, task, index, begin, end, len * E, r);
				met.Solve();
				
				std::cout << "Problem[" << index << "]" << std::endl;
				met.PrintTrueValue(index, task);
				std::cout << "the best value on " << met.GetBestIndex() << " iterator" << std::endl;
				std::cout << "x*=" << met.GetOpt().x << "  " << "F(x*)=" << met.GetOpt().z << std::endl;
				std::cout << std::endl;

				bool check = (check_method == 0) ? met.GetOpt().z <= len * E : Checked_method(met.GetOpt().x, met.GetTrueOpt(index, task), len * E);
				if (check)
				{
					count_true ++;
					out << met.GetBestIndex() << std::endl;
				}
				else
					out << "wrong" << std::endl;
				met.Clear();
			}
			out << "count_true " << count_true << std::endl;
			out << "Test_finish" << std::endl;
			break;
		}
		case 3:
		{
			n = 2;
			count_true = 0;
			MethodMult met;
			double a = 0.0, b = 1.0;
			int task = 1; // 0 - Grishagin, 1 - GKLS
			int len = 1; //  1 - Grishagin, 2 - GKLS, 1 - HILL, 10 - Shekel
			std::cout << "[0, 3] - Grishagin, [4, 7] - GKLS\n";
			std::cin >> task;
			int check_method = 0;
			std::cout << "check_method: 0- z_min, 1 - curr_eps\n";
			std::cin >> check_method;

			if (task >= 4 && task <= 7) {
				a = -1.0;
				len = 2;
			}
			double* y = new double[n];
			std::vector<std::vector<int>> matrix1 = { { 1 } };
			std::vector<std::vector<int>> matrix2 = { {2} };
			std::vector<std::vector<std::vector<int>>> matrix_res = { {{0}} };
			out << "[" << key << "] Task" << task << " e " << len * E << " r " << r << " check " << check_method << std::endl;
			double aver = 0;
			for (int index = 50; index < 51; index++)
			{
				met.Init(task, check_method, index, y, a, b, len * E, r, n, m);
				std::chrono::system_clock::time_point start = std::chrono::system_clock::now();
				met.SolveMult(y, matrix1, matrix2, matrix_res, 1);
				std::chrono::system_clock::time_point end = std::chrono::system_clock::now();
				auto sec = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
				average_time += sec.count();
				std::cout << "GrishaginProblem[" << index << "]" << std::endl;
				met.PrintTrueValue();
				std::cout << "the best value on " << met.GetBestIndex() << " iterator" << std::endl;
				std::cout << "x*= " << met.GetOpt()[0] << " y*= " << met.GetOpt()[1] << "  " << " F(x*)= " << met.GetOpt()[2] << std::endl;
				std::cout << std::endl;
				bool check = (check_method == 0) ? met.GetOpt()[2] <= len * E : Checked_method_mult(n, met.GetOpt(), met.GetTrueOpt(), len * E);
				if (check) //Checked_method_mult(n, met.GetOpt(), met.GetTrueOpt(), len * E) met.GetOpt()[2] <= len * E
				{
					count_true ++;
					aver += met.GetBestIndex();
					out << met.GetBestIndex() << std::endl;
				}
				else
					out << "wrong " << index << std::endl;
				met.ClearMethod();
			}
			std::cout << " Average Time = " << average_time / 100  << " count_true " << count_true << "\n";
			aver /= count_true;
			out << "aver = " << aver << "\n";
			out << "count_true " << count_true << std::endl;
			std::cout << "aver = " << aver << "\n";
			std::cout << "count_true " << count_true << std::endl;
			//out << "Mult_finish" << std::endl;
			delete[]y;
			break;
		}
		case 4:
		{
			n = 2;
			int p = 4;
			count_true = 0;
			MethodMult_p met;
			//out << "Mult_p" << std::endl;
			double* y = new double[n];
			int task = 1; // 0 - Grishagin, 1 - GKLS
			int len = 1; //  1 - Grishagin, 2 - GKLS, 1 - HILL, 10 - Shekel
			double a = 0.0, b = 1.0;
			std::cout << "0 - Grishagin, 1 - GKLS\n";
			std::cin >> task;

			if (task) {
				a = -1.0;
				len = 2;
			}
			out << "[" << key << "] Task" << task << " e " << len * E << " r " << r << std::endl;
			for (int index = 0; index < 1; index++)
			{
				met.Init_p(task, index, y, a, b, len * E, r, n, m, p);
				met.SolveMult_p(y);
				std::cout << "GrishaginProblem[" << index << "]" << std::endl;
				met.PrintTrueValue();
				std::cout << "the best value on " << met.GetBestIndex() << " iterator" << std::endl;
				std::cout << "x*= " << met.GetOpt()[0] << " y*= " << met.GetOpt()[1] << "  " << " F(x*)= " << met.GetOpt()[2] << std::endl;
				std::cout << std::endl;
				if (Checked_method_mult(n, met.GetOpt(), met.GetTrueOpt(), len * E))
				{
					count_true++;
					out << met.GetBestIndex() << std::endl;
				}
				else
					out << "wrong" << std::endl;
				met.ClearMethod();
			}
			//out << "count_true " << count_true << std::endl;
			//out << "Mult_p_finish" << std::endl;
			delete[]y;
			break;
		}
		case 5:
		{
			n = 2;
			int _step =250;
			int _alpha = 15;
			Mixture _mixed(0,0);
			MethodMult_mixed method;
			std::cin >> _mixed.alg1 >> _mixed.alg2;
			std::cout << _mixed.alg1 << " " << _mixed.alg2 << std::endl;
			count_true = 0;
			//out << "Mult_mix" << std::endl;
			int index = 0;
			int task = 1; // 0 - Grishagin, 1 - GKLS
			int len = 1; //  1 - Grishagin, 2 - GKLS, 1 - HILL, 10 - Shekel
			double a = 0.0, b = 1.0;
			std::cout << "0 - Grishagin, 1 - GKLS\n";
			std::cin >> task;
			if (task) {
				a = -1.0;
				len = 2;
			}
			double* y = new double[n];
			out << "[" << key << "] Task" << task << " e " << len * E << " r " << r << std::endl;
			for (int index = 0; index < 1000; index++)
			{
				//auto  metod = std::make_shared<MethodMult_mixed>(task, index, y, -1, 1, E, r, n, m, _step, _mixed, _alpha);
				method.Init_Mix(task, index, y, a, b, len * E, r, n, m, _step, _mixed, _alpha);
				method.SolveMult_mixed(y);
				std::cout << "GrishaginProblem[" << index << "]" << std::endl;
				method.PrintTrueValue();
				std::cout << "the best value on " << method.GetBestIndex() << " iterator" << std::endl;
				std::cout << "x*= " << method.GetOpt()[0] << " y*= " << method.GetOpt()[1] << "  " << " F(x*)= " << method.GetOpt()[2] << std::endl;
				std::cout << std::endl;
				if (Checked_method_mult(n, method.GetOpt(), method.GetTrueOpt(), len * E))
				{
					count_true++;
					out << method.GetBestIndex() << std::endl;
				}
				else
					out << "wrong" << std::endl;
			}
			//met.CleaMethod();
			//out << "count_true " << count_true << std::endl;
			//out << "Mult_mixed_finish" << std::endl;
			delete[]y;
			break;
		}
		case 6:
		{
			n = 2;
			int _step = 250;
			int _alpha = 15;
			int p = 2;
			MethodMult_mixed_p met;
			Mixture _mixed(1, 2);
			count_true = 0;
			out << "Mult_mix_p" << std::endl;
			int index = 0;
			int task = 1; // 0 - Grishagin, 1 - GKLS
			int len = 1; //  1 - Grishagin, 2 - GKLS, 1 - HILL, 10 - Shekel
			std::cout << "0 - Grishagin, 1 - GKLS\n";
			std::cin >> task;
			double a = 0.0, b = 1.0;
			if (task) {
				a = -1.0;
				len = 2;
			}
			double* y = new double[n];
			out << "[" << key << "] Task" << task << " e " << len * E << " r " << r << std::endl;
			for (int index = 0; index < 1000; index++)
			{
				met.Init_Mix_p(task, index, y, a, b, len * E, r, n, m, _step, _mixed, _alpha, p);
				met.SolveMult_mixed_p(y);
				std::cout << "GrishaginProblem[" << index << "]" << std::endl;
				met.PrintTrueValue();
				std::cout << "the best value on " << met.GetBestIndex() << " iterator" << std::endl;
				std::cout << "x*= " << met.GetOpt()[0] << " y*= " << met.GetOpt()[1] << "  " << " F(x*)= " << met.GetOpt()[2] << std::endl;
				std::cout << std::endl;
				if (Checked_method_mult(n, met.GetOpt(), met.GetTrueOpt(), len * E))
				{
					count_true++;
					out << met.GetBestIndex() << std::endl;
				}
				else
					out << "wrong" << std::endl;
				met.ClearMethod();
			}
			out << "count_true " << count_true << std::endl;
			out << "Mult_mixed_p_finish" << std::endl;
			delete[]y;
			break;
		}case 7:
		{
			n = 2;
			count_true = 0;
			MethodMult_DualLipsh met;
			//out << "Mult_DualLipsh_p" << std::endl;
			int index = 0;
			int task = 1; // 0 - Grishagin, 1 - GKLS
			int len = 1; //  1 - Grishagin, 2 - GKLS, 1 - HILL, 10 - Shekel
			double a = 0.0, b = 1.0;
			std::cout << "0 - Grishagin, 1 - GKLS\n";
			std::cin >> task;
			if (task == 1) {
				a = -1.0;
				len = 2;
			}
			double r_loc = 1.48;
			/*std::cout << "r_loc = " << std::endl;
			std::cin >> r_loc;*/
			double* y = new double[n];
			aver = 0;
			out << "[" << key << "] Task" << task << " e " << len * E << " r " << r << std::endl;
			for (int index = 0; index < 100; index++)
			{
				met.Init_Dual(task, index, y, a, b, len * E, r_loc, r, n, m);
				met.SolveMult_DualLipsh(y);
				std::cout << "GrishaginProblem[" << index << "]" << std::endl;
				met.PrintTrueValue();
				std::cout << "the best value on " << met.GetBestIndex() << " iterator" << std::endl;
				std::cout << "x*= " << met.GetOpt()[0] << " y*= " << met.GetOpt()[1] << "  " << " F(x*)= " << met.GetOpt()[2] << std::endl;
				std::cout << std::endl;
				if (Checked_method_mult(n, met.GetOpt(), met.GetTrueOpt(), len * E))
				{
					count_true++;
					aver += met.GetBestIndex();
					out << met.GetBestIndex() << std::endl;
				}
				else
					out << "wrong" << std::endl;
				met.ClearMethod();
			}
			out << "count_true " << count_true << " aver= " << aver / count_true << std::endl;
			//out << "Mult_DualLipsh_p_finish" << std::endl;
			delete[]y;
			break;
		}
		case 8:
		{
			//out << "Mult_DualLipsh" << std::endl;
			int task = 3; // 0 - Grishagin, 1 - GKLS, 2 - HILL, 3 - Shekel
			std::cout << "0 - Grishagin, 1 - GKLS, 2 - HILL, 3 - Shekel\n";
			std::cin >> task;
			int len = 1; //  1 - Grishagin, 2 - GKLS, 1 - HILL, 10 - Shekel
			double a = 0.0, b = 1.0;
			if (task == 1) {
				a = -1.0;
				len = 2;
			}
			else if (task == 3) {
				b = 10.0;
				len = 10;
			}
			double r_loc;
			double r_glob = r;
			double eps = 0.01;
			std::vector<double> z_iter;
			double* y = new double[n];
			int k = 10;
			int p = 10;
			
			std::ofstream out_log;
			out_log.open("LogFile.txt", std::ofstream::ios_base::app);
			
			out_log << "k = " << k << " p = " << p << std::endl;
			out_log << "r_glob |" << " r_loc   |" << " task" << std::endl;
		
			double h_glob = r_glob / k;
			double h_loc = (1 - 1 / r_glob - 2 * eps) / p;
			double alpha = 1 / r_glob + eps;
			double a_0 = alpha;
			r_loc = alpha * r_glob;

			MethodMult_DualLipsh met;

			/*std::ofstream out2;
			out2.open("3D_graph.txt", std::ofstream::ios_base::app);*/

			for (int j = 0; j <= p; j++)  //p+1/2
			{
				alpha = a_0 + j * h_loc;
				for (int i = 0; i <= k; i++)														
				{																					
					r_glob = r + i * h_glob;														
					r_loc = alpha * r_glob;															
					double iter_count = 0;															
					count_true = 0;																	
					out_log << std::fixed << std::setprecision(3) << r_glob << "  | ";				
					out_log << std::fixed << std::setprecision(3) << r_loc << "   | ";	
					std::cout << i << std::endl;
																									
					for (int index = 0; index < 100; index++)										
					{																				
						met.Init_Dual(task, index, y, a, b, len * E, r_loc, r_glob, n, m); //GKLS [-1; 1] 
						met.SolveMult_DualLipsh(y);
						std::cout << index << std::endl;

						if (Checked_method_mult(n, met.GetOpt(), met.GetTrueOpt(), len * E))
						{
							count_true++;
							iter_count += met.GetBestIndex();
						}
						else {
							//out << "wrong" << std::endl;
						}
						met.ClearMethod();
					}
					out_log << count_true<< std::endl;
					z_iter.push_back(iter_count / count_true);
					//out2 << iter_count / count_true<<" ";
				}
				//out2 << std::endl;
			}

			std::ofstream out2;
			out2.open("3DD_graph.txt", std::ofstream::ios_base::app);
			for (int j = 0; j <= p; j++) {
				for (int i = 0; i <= k; i++)
					out2 << z_iter[j * (k + 1) + i] << " ";
				out2 << std::endl;
			}
			out2.close();
			out_log.close();
			delete[]y;
			break;
		}
		case 9:
		{
			//count_true = 0;
			//MethodMult_DualLipsh_p met;
			////out << "Mult_DualLipsh_p" << std::endl;
			//int index = 0;
			//int task = 0; // 0 - Grishagin, 1 - GKLS
			//double a = 0.0, b = 1.0;
			//if (task)
			//	a = -1.0;
			//double r_loc = 1.7;
			//int p = 4;
			///*std::cout << "r_loc = " << std::endl;
			//std::cin >> r_loc;*/
			//double* y = new double[n];
			//for (int index = 0; index < 1; index++)
			//{
			//	met.Init_Dual_p(task, index, y, a, b, E, r_loc, r, n, m, p);
			//	met.SolveMult_DualLipsh_p(y);
			//	std::cout << "GrishaginProblem[" << index << "]" << std::endl;
			//	met.PrintTrueValue();
			//	std::cout << "the best value on " << met.GetBestIndex() << " iterator" << std::endl;
			//	std::cout << "x*= " << met.GetOpt()[0] << " y*= " << met.GetOpt()[1] << "  " << " F(x*)= " << met.GetOpt()[2] << std::endl;
			//	std::cout << std::endl;
			//	if (Checked_method_mult(n, met.GetOpt(), met.GetTrueOpt(), E))
			//	{
			//		count_true++;
			//		out << met.GetBestIndex() << std::endl;
			//	}
			//	else
			//		out << "wrong" << std::endl;
			//	met.ClearMethod();
			//}
			////out << "count_true " << count_true << std::endl;
			////out << "Mult_DualLipsh_p_finish" << std::endl;
			//delete[]y;
			//break;
		}
		case 10: {
			count_true = 0;
			aver = 0;
			Method_SearchRoot met;
			int task = 2;         // 0- Hans, 1 - Hill, 2 - Shekel
			int len = 1; //  1 - Grishagin, 2 - GKLS, 1 - HILL, 10 - Shekel
			std::cout << "[0-2] - Hans, [3-5] - Hill, [6-8] - Shekel\n";
			std::cin >> task;

			int check_method = 0;
			std::cout << "check_method: 0- z_min, 1 - curr_eps\n";
			std::cin >> check_method;

			vector<double> begin(20, 0), end(20, 0);
			int maxTask = 20;
			if (task >= 3 and task <= 5) {
				// Hill на [0; 1]
				begin = { 0.0 };
				end = { 1.0 };
				maxTask = 1000; // 1000
			}
			else if (task >= 6 and task <= 8) {
				// Shekel на [0; 10]
				len = 10;
				begin = { 0.0 };
				end = { 10.0 };
				maxTask = 1000;
			}

			out << "[" << key << "] Task" << task << " e " << len * E << " r " << r << " check " << check_method << std::endl;
			for (int index = 0; index < maxTask; index++) {

				met.Init(1,check_method, task, index, begin, end, len * E);
				met.Solve();
				met.PrintTrueValue(index, task);

				std::cout << "Problem[" << index << "]" << std::endl;
				// вывод истинных значений x,y задачи index
				//met.PrintTrueValue(index, task);

				// вывод итерации, на которой было получено лучшее значение
				std::cout << "the best value on " << met.GetBestIndex() << " iterator" << std::endl;

				// вывод полученных  значений x,y задачи index
				std::cout << "x*=" << met.GetOpt().x << "  " << "F(x*)=" << met.GetOpt().z << std::endl;
				std::cout << std::endl;
				bool check = (check_method == 0) ? met.GetOpt().z <= len * E : Checked_method(met.GetOpt().x, met.GetTrueOpt(index, task), len * E);
				if (check) // met.GetOpt().z <= len*E Checked_method(met.GetOpt().x, met.GetTrueOpt(index, task), E)
				{
					count_true++;
					aver += met.GetBestIndex();
					out << met.GetBestIndex() << std::endl;
				}
				else {
					out << "wrong " << index << std::endl;
				}
				met.Clear();
			}
			out << "count_true " << count_true << " aver= "<< aver/count_true<< std::endl;
			std::cout << "count_true " << count_true << " aver= " << aver / count_true << std::endl;
			out << "Test_func_finish" << std::endl;
			break;
		}
		case 11: {
			count_true = 0;
			aver = 0;
			Method_AGP_ASR met;
			int task = 1;
			std::cout << "[0-2] - Hans, [3-5] - Hill, [6-8] - Shekel\n";
			int len = 1; //  1 - Grishagin, 2 - GKLS, 1 - HILL, 10 - Shekel
			std::cin >> task;

			int check_method = 0;
			std::cout << "check_method: 0- z_min, 1 - curr_eps\n";
			std::cin >> check_method;

			vector<double> begin(20, 0), end(20, 0);
			int maxTask = 20;
			if (task >= 3 and task <= 5) {
				// Hill на [0; 1]
				begin = { 0.0 };
				end = { 1.0 };
				maxTask = 1000; // 1000
			}
			else if (task >= 6 and task <= 8) {
				// Shekel на [0; 10]
				len =  10;
				begin = { 0.0 };
				end = { 10.0 };
				maxTask = 1000;
			}

			out << "[" << key << "] Task" << task << " e " << len * E << " r " << r << " check " << check_method << std::endl;
			for (int index = 0; index < maxTask; index++) {

				met.PrintTrueValue(index, task);
				met.Init(1, check_method, task, index, begin, end, len * E, r);
				met.Solve();

				std::cout << "Problem[" << index << "]" << std::endl;
				// вывод истинных значений x,y задачи index
				
				met.PrintTrueValue(index, task);
				// вывод итерации, на которой было получено лучшее значение
				std::cout << "the best value on " << met.GetBestIndex() << " iterator" << std::endl;

				// вывод полученных  значений x,y задачи index
				std::cout << "x*=" << met.GetOpt().x << "  " << "F(x*)=" << met.GetOpt().z << std::endl;
				std::cout << std::endl;


				bool check = (check_method == 0) ? met.GetOpt().z <= len * E : Checked_method(met.GetOpt().x, met.GetTrueOpt(index, task), len * E);
				if (check)  // met.GetOpt().z <= len *E Checked_method(met.GetOpt().x, met.GetTrueOpt(index, task), E)
				{
					count_true++;
					aver += met.GetBestIndex();
					out << met.GetBestIndex() << std::endl;
				}
				else {
					out << "wrong " << index << std::endl;
				}
				met.Clear();
			}
			out << "count_true " << count_true << " aver= " << aver / count_true << std::endl;
			std::cout << "count_true " << count_true << " aver= " << aver / count_true << std::endl;
			out << "Test_func_finish" << std::endl;
			break;
		}
		case(12): {
			n = 2;
			count_true = 0;
			MethodMult_AGP_ASR met;
			double a = 0.0, b = 1.0;
			int task = 1; // 0 - Grishagin, 1 - GKLS
			int len = 1; //  1 - Grishagin, 2 - GKLS, 1 - HILL, 10 - Shekel
			std::cout << "[0, 3] - Grishagin, [4, 7] - GKLS\n";
			std::cin >> task;

			int check_method = 0;
			std::cout << "check_method: 0- z_min, 1 - curr_eps\n";
			std::cin >> check_method;

			if (task >= 4 and task <= 7) { // ?
				a = -1.0;
				len = 2;
			}
			double* y = new double[n];

			out << "[" << key << "] Task" << task << " e " << len * E << " r " << r << " check " << check_method << std::endl;
			double aver = 0;
			for (int index = 50; index < 51; index++)
			{
				met.Init(task,check_method, index, y, a, b, len * E, r, n, m);
				std::chrono::system_clock::time_point start = std::chrono::system_clock::now();
				met.SolveMult_AGP_ASR(y);
				std::chrono::system_clock::time_point end = std::chrono::system_clock::now();
				auto sec = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
				average_time += sec.count();
				std::cout << "Problem[" << index << "]" << std::endl;
				met.PrintTrueValueMult();
				std::cout << "the best value on " << met.GetBestIndex() << " iterator" << std::endl;
				std::cout << "x*= " << met.GetOpt()[0] << " y*= " << met.GetOpt()[1] << "  " << " F(x*)= " << met.GetOpt()[2] << std::endl;
				std::cout << std::endl;

				bool check = (check_method == 0) ? met.GetOpt()[2] <= len * E : Checked_method_mult(n, met.GetOpt(), met.GetTrueOpt(), len * E);
				if (check) // Checked_method_mult(n, met.GetOpt(), met.GetTrueOpt(), len * E) met.GetOpt()[2] <= len * E
				{
					count_true++;
					aver += met.GetBestIndex();
					out << met.GetBestIndex() << std::endl;
				}
				else
					out << "wrong " << index << std::endl;
				met.ClearMethod();
			}
			std::cout << "Average Time = " << average_time / 100 << " count_true " << count_true << "\n";
			aver /= count_true;
			out << "aver = " << aver << "\n";
			out << "count_true " << count_true << std::endl;
			std::cout << "aver = " << aver << "\n";
			std::cout << "count_true " << count_true << std::endl;
			//out << "Mult_finish" << std::endl;
			delete[]y;
			break;
		}
		case 13: {
			n = 2;
			count_true = 0;
			aver = 0;
			MethodMult_SearchRoot met;
			double a = 0.0, b = 1.0;
			int task = 2;
			std::cout << "[0, 3] - Grishagin, [4, 7] - GKLS\n";
			int len = 1; //  1 - Grishagin, 2 - GKLS
			double* y = new double[n];
			std::cin >> task;

			int check_method = 0;
			std::cout << "check_method: 0- z_min, 1 - curr_eps\n";
			std::cin >> check_method;

			if (task > 3) {
				a = -1.0;
				len = 2;
			}

			out << "[" << key << "] Task" << task << " e " << len * E << " r " << r << " check " << check_method << std::endl;
			for (int index = 0; index < 100; index++) {

				met.Init(task,check_method, index, y, a, b, len* E, r, n, m);
				met.SolveMult_SR(y);

				std::cout << "Problem[" << index << "]" << std::endl;
				// вывод истинных значений x,y задачи index
				met.PrintTrueValueMult();

				// вывод итерации, на которой было получено лучшее значение
				std::cout << "the best value on " << met.GetBestIndex() << " iterator" << std::endl;

				std::cout << "x*= " << met.GetOpt()[0] << " y*= " << met.GetOpt()[1] << "  " << " F(x*)= " << met.GetOpt()[2] << std::endl;
				std::cout << std::endl;

				bool check = (check_method == 0) ? met.GetOpt()[2] <= len * E : Checked_method_mult(n, met.GetOpt(), met.GetTrueOpt(), len * E);
				if (check)// met.GetOpt()[2] <= len*E Checked_method_mult(n, met.GetOpt(), met.GetTrueOpt(), len * E)
				{
					count_true++;
					aver += met.GetBestIndex();
					out << met.GetBestIndex() << std::endl;
				}
				else {
					out << "wrong " << index << std::endl;
				}
				met.ClearMethod();
			}
			out << "count_true " << count_true << " aver= " << aver / count_true << std::endl;
			std::cout << "count_true " << count_true << " aver= " << aver / count_true << std::endl;
			out << "Test_func_finish" << std::endl;
			break;
		}
		}
	}
	out.close();
// int count = 400;
//double a = -1;
//double b = 1;
//double h = (b - a) / count;
//std::vector<double> x, y, z;
//grid_build(x, y, z, a, b, 0, h);
	_getch();
}
