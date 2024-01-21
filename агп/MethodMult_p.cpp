#include"MethodMult_p.h"

#include <iostream>
#include<fstream>


void MethodMult_p::Init_p(int _task, int _index_problem, double* y, double _a, double _b,
	double _e, double _r, int _n, int _m, int _p){
	a = _a;
	b = _b;
	p = _p;
	eps = _e;
	r = _r;
	n = _n;
	m = _m;
	task = _task;
	index_problem = _index_problem;

	Trial_thread current, first, second;
	
	best_i = 0;
	out_optimal = { 0, 0, 0 };  // оптимальное решение
	std::mt19937 gen;
	gen.seed(static_cast<unsigned int>(time(0)));
	double x_next;

	first.x = 0;
	first.thread = 0;
	mapd(0, m, y, n, 1);
	ScaleFunc(y[0]);
	InsertScale(y);
	// вычисляем значение в этой точке
	first.z = Funk_mult(y);
	trials_thread.push_back(first);
	
	/*current.x = 0.5;
	mapd(current.x, m, y, n, 1);
	InsertScale(y);
	current.z = Funk_mult(y);
	trials_thread.push_back(current);

	current.x = 0.6;
	mapd(current.x, m, y, n, 1);
	InsertScale(y);
	current.z = Funk_mult(y);
	trials_thread.push_back(current);

	current.x = 0.8;
	mapd(current.x, m, y, n, 1);
	InsertScale(y);
	current.z = Funk_mult(y);
	trials_thread.push_back(current);*/

	//потоковая реализация, необходимо p интервалов (есть 1 либо 2 - если х_1 отлична от 0 и 1) нужно отсортировать точки 
	std::vector<double> next_point(p - size_t(1));
	for (int j = 0; j < p - 1; j++) {
		x_next = gen() % 100 / (100 * 1.0);
		while (x_next == 1 || x_next == 0) 
		{  // если выбранная точка x_next - концы интервала
			x_next = gen() % 100 / (100 * 1.0);
		}
		next_point[j] = x_next;
	}
	sort(next_point.begin(), next_point.end());

	
	for (int j = 0; j < p - 1; j++) {
		current.x = next_point[j];  //(j + 1) / double(p);
		current.thread = 0;
		mapd(current.x, m, y, n, 1);
		InsertScale(y);
		current.z = Funk_mult(y);
		trials_thread.push_back(current);
	}
	
	second.x = 1;
	second.thread = 0;
	mapd(1, m, y, n, 1);
	InsertScale(y);
	second.z = Funk_mult(y);
	trials_thread.push_back(second);
}

void MethodMult_p::ClearMethod_p() {
	trials_thread.clear();
	out_optimal.clear();
}

void MethodMult_p::SolveMult_p(double* y)
{
	Trial current;  // для подсчета нового испытания
	double M;
	int itr = 0;  // счетчик итераций
	std::vector<Characteristic> charact(p);
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

	std::ofstream out1;
	out1.open("Grishagin_p.txt", std::ofstream::ios_base::app);

	std::vector<double> true_opt = GetTrueOpt();

	while (fabs(true_opt[0] - out_optimal[0]) > eps || fabs(true_opt[1] - out_optimal[1]) > eps)  //while (curr_eps > eps)
	{
		// начальные значение для вычисления M на интервале 1
		double d_z = fabs(trials_thread[1].z - trials_thread[0].z);
		double d_x = fabs(trials_thread[1].x - trials_thread[0].x);
		d_x = pow(d_x, power);
		M = d_z / d_x;
		//M = 50;

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

		double k, R;

		for (size_t i = 1; i < trials_thread.size(); i++)  // поиск со 2 интервала
		{
			k = pow(trials_thread[i].x - trials_thread[i - 1].x, power);  // для оптимизации вычисления

			R = k + (pow((trials_thread[i].z - trials_thread[i - 1].z) / (M * r), 2) / k) -
				2 * (trials_thread[i].z + trials_thread[i - 1].z - 2 * z_min) / (r * M);
			// обновляем каждый раз до кол-ва потоков, после выбираем лучшие если есть из чего
			if (i <= p)
			{
				charact[i - 1].R = R;
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

				current.x = (trials_thread[charact[pos].pos].x + trials_thread[charact[pos].pos - size_t(1)].x) / 2 - sgn / (2 * r) * pow(delta_z_pos / M, n);

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
			new_trial.z = Funk_mult(y);

			std::vector<Trial_thread>::iterator it2 = find(trials_thread.begin(), trials_thread.end(), elem_of_ch[j]);
			size_t pos_elem_of_ch = std::distance(trials_thread.begin(), it2);
			trials_thread.insert(it2, new_trial);

			out1 << y[0] << " " << y[1] << " " << new_trial.thread << std::endl;

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
	out1.close();
}
