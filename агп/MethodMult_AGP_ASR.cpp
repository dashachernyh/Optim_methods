#include "MethodMult_AGP_ASR.hpp"
#include <fstream>

// f(x) >= 0
void MethodMult_AGP_ASR::SolveMult_AGP_ASR(double* y) {
	Trial current;
	double M, Rmax, max, R, R0;
	double ro;
	size_t Rpos, zpos = 1;
	int flag = 1; // для первой итерации
	// M Strongina -> [x0, x1], [a, b] != [0, 1]

	double z_min = out_optimal[n];
	double power = 1 / double(n);
	double curr_eps = pow(trials[1].x - trials[0].x, power);
	int itr = 0;
	std::vector<double> true_opt = GetTrueOpt();
	std::ofstream out1;
	out1.open("Log_AGP_ASR.txt", std::ofstream::ios_base::app);

	//optimum.z > eps curr_eps > eps
	while (fabs(true_opt[0] - out_optimal[0]) > eps || fabs(true_opt[1] - out_optimal[1]) > eps)
	{
		//out1 << "i = " << itr << " z_min = " << out_optimal[n] <<"\n";
		Rpos = 1;

		double d_z = fabs(trials[1].z - trials[0].z);
		double d_x = fabs(trials[1].x - trials[0].x);
		d_x = pow(d_x, power);
		M = d_z / d_x;

		for (size_t i = 2; i < trials.size(); i++)
		{
			max = fabs((trials[i].z - trials[i - 1].z)) / pow(trials[i].x - trials[i - 1].x, power);
			if (max > M)
				M = max;
		}

		if (M == 0)
			M = 1;
		//out1 << "\tM= " << M << "\n";
		if (flag) {
			double dt = pow(trials[zpos].x - trials[zpos - 1].x, power);
			ro = dt * pow(1 - (trials[zpos].z + trials[zpos - 1].z) / (r * M *dt), 2) + 4*z_min/(r* M);
			flag = 0;
		}

		//out1 << "\tro = " << ro << "\n";

		d_x = pow(trials[1].x - trials[0].x, power);
		R0 = -4 * (trials[1].z * trials[0].z) / (pow(r * M, 2) * d_x);
		R = d_x * pow(1 - (trials[1].z + trials[0].z) / (r * M * d_x), 2)
			+ R0 + 4 * z_min / (r * M);

		Rmax = std::max(R, R0 + ro);

		for (int i = 2; i < trials.size(); i++)
		{
			d_x = pow(trials[i].x - trials[i - 1].x, power);
			R0 = -4 * (trials[i].z * trials[i - 1].z) /(pow(r * M, 2) * d_x);
			R = d_x * pow(1 - (trials[i].z + trials[i - 1].z) / (r * M * d_x), 2)
				+ R0 + 4 * z_min / (r * M);

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
		//out1 << "\tRmax = " << Rmax << " pos = " << Rpos << "\n";
		curr_eps = trials[Rpos].x - trials[Rpos - size_t(1)].x;

		std::vector<Trial>::iterator it2 = trials.begin();
		it2 += Rpos;

		double delta_z_pos = trials[Rpos].z - trials[Rpos - size_t(1)].z;

		int sgn = 0;
		if (delta_z_pos < 0)
			sgn = -1;
		if (delta_z_pos > 0)
			sgn = 1;

		current.x = (trials[Rpos].x + trials[Rpos - size_t(1)].x) / 2
			- sgn / (2 * r) * pow(fabs(delta_z_pos) / M, n);
		mapd(current.x, m, y, n, 1);
		InsertScale(y);

		current.z = Funk_mult(y);
		trials.insert(it2, current);

		//out1 << "\tx = " << current.x << " z = " << current.z << "\n";

		if (z_min > current.z)
		{
			z_min = current.z;
			zpos = Rpos;
			best_i = itr;
			flag = 1;
			out_optimal[0] = y[0];
			out_optimal[1] = y[1];
			out_optimal[2] = z_min;
		}
		else if (zpos == Rpos) {
			flag = 1;
			++zpos;
		}

		//out1 << "zmin = " << z_min << " zpos = " << zpos << "\n";
		itr++;
	}
	//out1.close();
	std::cout << "itr = " << itr << std::endl;
}