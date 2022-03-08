#include "MethodMult_DualLipsh_p.h"

#include <iostream>
#include <fstream>

void MethodMult_DualLipsh_p::SolveMult_DualLipsh_p(double* y)
{

	Trial current;  // для подсчета нового испытания
	double M, R_loc, R_glob;
	MaxR Rmax;
	int itr = 0; // счетчик итераций
	std::vector<Characteristic_dual> charact(p);
	std::vector<Trial_thread>::iterator it = trials_thread.begin();  // для поиска позиции добавления нового испытания в векторе trials
	
	// минимальное значение функции
	double z_min = trials_thread[0].z;
	for (size_t i = 1; i < trials_thread.size(); i++) {
		if (z_min > trials_thread[i].z)
			z_min = trials_thread[i].z;
	}
	
	double power = 1 / double(n);
	double curr_eps = pow(trials_thread[1].x - trials_thread[0].x, power);
	out_optimal[2] = trials_thread[0].z;

	double ro = (1 - 1 / r_glob) / (1 - 1 / r_loc) * (1 - 1 / r_glob) / (1 - 1 / r_loc);

	//std::ofstream out1;
	//out1.open("Grishagin.txt", std::ofstream::ios_base::app);
	std::vector<double> true_opt = GetTrueOpt(task, index_problem);

	while (fabs(true_opt[0] - out_optimal[0]) > eps || fabs(true_opt[1] - out_optimal[1]) > eps)  //while (curr_eps > eps)
	{
		// начальные значение для вычисления M на интервале 1
		double d_z = fabs(trials_thread[1].z - trials_thread[0].z);
		double d_x = fabs(trials_thread[1].x - trials_thread[0].x);
		d_x = pow(d_x, power);
		M = d_z / d_x;  //M = 50;

		for (size_t i = 2; i < trials_thread.size(); i++)  // поиск со 2 интервала
		{
			double max;
			max = fabs((trials_thread[i].z - trials_thread[i - size_t(1)].z))
				/ pow(trials_thread[i].x - trials_thread[i - size_t(1)].x, power);
			if (max > M)
				M = max;
			if (z_min > trials_thread[i].z)
				z_min = trials_thread[i].z;
		}
		if (M == 0)
			M = 1;

		double R, k;
		for (size_t i = 1; i < trials_thread.size(); i++)  // поиск с 1 интервала
		{
			k = pow(trials_thread[i].x - trials_thread[i - 1].x, power);  // для оптимизации вычисления

			R_loc = k + (pow((trials_thread[i].z - trials_thread[i - size_t(1)].z) / (M * r_loc), 2) / k) -
				2 * (trials_thread[i].z + trials_thread[i - size_t(1)].z - 2 * z_min) / (r_loc * M);
			R_glob = k + (pow((trials_thread[i].z - trials_thread[i - size_t(1)].z) / (M * r_glob), 2) / k) -
				2 * (trials_thread[i].z + trials_thread[i - size_t(1)].z - 2 * z_min) / (r_glob * M);

			if (ro * R_loc > R_glob)
			{
				R = ro * R_loc;
				r = r_loc;
			}
			else
			{
				R = R_glob;
				r = r_glob;
			}

			// обновляем каждый раз до кол-ва потоков, после выбираем лучшие если есть из чего
			if (i <= p)
			{
				charact[i - 1].R = R;
				charact[i - 1].r = r;
				charact[i - 1].pos = i;

				if (i == p)
				{
					sort(charact.begin(), charact.end());
				}
			}
			else {
				if (R > charact[0].R)
				{
					charact[0].R = R;
					charact[0].r = r;
					charact[0].pos = i;
					sort(charact.begin(), charact.end());
				}
			}
		}

		// выбираем минимальный curr_eps
		/*curr_eps = pow(trials_thread[charact[0].pos].x - trials_thread[charact[0].pos - size_t(1)].x, power);
		for (int i = 1; i < p; i++)
		{
			double new_curr = pow(trials_thread[charact[i].pos].x - trials_thread[charact[i].pos - size_t(1)].x, power);
			if (curr_eps > new_curr)
				curr_eps = new_curr;
		}*/

		//заполняем массив значениями правых границ 
		//интервала с лучшими р характеристиками 
		//для поиска этих элементов после добавления новых точек
		std::vector<Trial_thread> elem_of_ch(p);
		for (int i = 0; i < p; i++) {
			elem_of_ch[i] = trials_thread[charact[i].pos];
		}
		//можно параллелить
		double delta_z_pos, sgn;
		std::vector<Trial> vect_current(p);
		int pos = 0;

#pragma omp parallel num_threads(p) shared(vect_current) private(delta_z_pos, sgn, current, pos)
		{
#pragma omp for schedule(static)
			for (pos = 0; pos < p; pos++) {

				// вектор испытаний не изменяли, все концы интервалов находятся на своем месте

				delta_z_pos = trials_thread[charact[pos].pos].z - trials_thread[charact[pos].pos - size_t(1)].z;

				sgn = 0;
				if (delta_z_pos < 0)
					sgn = -1;
				if (delta_z_pos > 0)
					sgn = 1;

				current.x = (trials_thread[charact[pos].pos].x + trials_thread[charact[pos].pos - size_t(1)].x) / 2
					- sgn / (2 * charact[pos].r) * pow(delta_z_pos / M, n);

				// каждый поток сохраняет полученную точку в общие данные
				vect_current[pos].x = current.x;
				vect_current[pos].z = omp_get_thread_num();
			}
		}

		// передаем управление главному потоку 0
		for (int j = 0; j < p; j++) {
			Trial_thread new_trial;
			new_trial.x = vect_current[j].x;
			new_trial.thread = (int)vect_current[j].z;
			mapd(new_trial.x, m, y, n, 1);
			InsertScale(y);
			new_trial.z = Funk_mult(task, index_problem, y);
			std::vector<Trial_thread>::iterator it2 = find(trials_thread.begin(), trials_thread.end(), elem_of_ch[j]);
			size_t pos_elem_of_ch = std::distance(trials_thread.begin(), it2);
			trials_thread.insert(it2, new_trial);

			//out1 << y[0] << " " << y[1] <<" "<< new_trial.thread<<std::endl;

			if (out_optimal[2] > new_trial.z)
			{
				best_i = itr;
				out_optimal[0] = y[0];
				out_optimal[1] = y[1];
				out_optimal[2] = new_trial.z;
			}
		}

		itr++;
	}

	std::cout << "itr = " << itr << std::endl;
	//out1 << trials_thread.size() << std::endl;
	//out1.close();
}