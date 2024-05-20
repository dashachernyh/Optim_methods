#include "Method_AGP_ASR.hpp"
#include <fstream>

double func(double x) {
	return fabs(x - 0.345);
}

// f(x) >= 0
void Method_AGP_ASR::Solve() {
	Trial current;
	double M, Rmax, max, k, R, R0;
	double ro;
	size_t Rpos, zpos = 1;
	int flag = 1; // для первой итерации
	// M Strongina -> [x0, x1], [a, b] != [0, 1]

	double curr_eps = trials[1].x - trials[0].x;
	int itr = 0;
	std::ofstream out1;
	out1.open("Log_AGP_ASR.txt", std::ofstream::ios_base::app);

	double mult_z = trials[1].z * trials[0].z;
	double sum_z = trials[1].z + trials[0].z;
	double dif_x = trials[1].x - trials[0].x;

	double check = (check_method == 0) ? optimum.z : curr_eps;
	while (check > eps && itr <= 10000) //optimum.z > eps curr_eps > eps
	{
		Rpos = 1;
		M = fabs((trials[1].z - trials[0].z) / (trials[1].x - trials[0].x));
		for (size_t i = 2; i < trials.size(); i++)
		{
			max = fabs((trials[i].z - trials[i - 1].z) / (trials[i].x - trials[i - 1].x));
			if (max > M)
				M = max;
		}

		M = M == 0 ? 1 : r * M;

		if (flag) {
			double m_dt = M * (trials[zpos].x - trials[zpos - 1].x);
			// было без + 1 ?????
			ro = m_dt * pow(1 - (trials[zpos].z + trials[zpos - 1].z) / m_dt, 2);
			flag = 0;
		}

		out1 << "ro = " << ro << "\n";
		
		k = M * (trials[1].x - trials[0].x);
		R0 = -4 * (trials[1].z * trials[0].z) / k;
		R = k * pow(1 - (trials[1].z + trials[0].z) / k, 2) + R0;
		
		// было R0 + ro ??????
		Rmax = std::max(R, R0 + ro);
		//Rmax = R > R0 * ro ? R : R0 * ro;

		for (int i = 2; i < trials.size(); i++)
		{
			k = M * (trials[i].x - trials[i - 1].x);
			R0 = -4 * (trials[i].z * trials[i - 1].z) / k;
			R = k*pow((1 - (trials[i].z  + trials[i - 1].z)/k), 2) + R0;

			// было без *
			R0 += ro;

			if (R >= R0 && R > Rmax)
			{
				Rmax = R;
				Rpos = i;
			}
			else if (R0 > R && R0 > Rmax) {
				Rmax = R0;
				Rpos = i;
			}

		}
		
		//out1 << "Rmax = " << Rmax << " Rpos =" << Rpos << "\n";

		curr_eps = trials[Rpos].x - trials[Rpos - size_t(1)].x;

		std::vector<Trial>::iterator it2 = trials.begin();
		it2 += Rpos;
	
		current.x = (trials[Rpos].x + trials[Rpos - 1].x) * 0.5 - (trials[Rpos].z
			- trials[Rpos - size_t(1)].z) / (2 * M);
		current.z = Funс(1, key, index_problem, current.x);
		trials.insert(it2, current);

		//out1 << "x = " << current.x << " z = " << current.z << "\n";

		if (optimum.z > current.z)
		{
			zpos = Rpos;
			best_i = itr;
			flag = 1;
			optimum = current;
		}
		else if (zpos == Rpos) {
			flag = 1;
			++zpos;
		}

		//out1 << "zmin = " << z_min << " zpos = " << zpos << "\n";
		itr++;
		check = (check_method == 0) ? optimum.z : curr_eps;
	}
	//out1.close();
	std::cout << "itr = " << itr << std::endl;
}