#pragma once

class MethodDual
{
	std::vector<Trial> trials;
	Trial optimum;
	Trial first;
	Trial second;
	int index_problem;
	int best_i;
	double eps, r;
public:
	int key;
	Trial GetOpt() { return optimum; }
	int GetBestIndex() { return best_i; }
	MethodDual(int _key, int index_problem, double x_0, double x_n, double _e, double _r);
	void solve_dual();
};

MethodDual::MethodDual(int _key, int _index_problem, double x_0, double x_n, double _e, double _r)
{
	index_problem = _index_problem;
	key = _key;
	first.x = x_0;
	first.z = Funk(key, index_problem, x_0);
	second.x = x_n;
	second.z = Funk(key, index_problem, x_n);
	trials.push_back(first);
	trials.push_back(second);
	best_i = 0;
	eps = _e;
	r = _r;

}
void MethodDual::solve_dual()
{
	Trial current;
	double M, Rmax, Rpos;
	double curr_eps = second.x - first.x;
	std::vector<Trial>::iterator it = trials.begin();
	double z_min = trials[0].z;
	if (trials[1].z < trials[0].z) z_min = trials[1].z;
	int itr = 0;
	optimum.z = first.z;
	while (curr_eps > eps)
	{
		Rpos = 1;
		double d_z = fabs(trials[1].z - trials[0].z);
		double d_x = trials[1].x - trials[0].x;
		M = d_z / d_x;
		for (int i = 2; i < trials.size(); i++)
		{
			double max;
			max = fabs((trials.at(i).z - trials.at(i - 1).z) / (trials.at(i).x - trials.at(i - 1).x));
			if (max > M)
				M = max;
			if (z_min > trials[i].z)
				z_min = trials[i].z;
		}

		if (M == 0)
			M = 1;
		else
			M = r * M;
		Rmax = d_x + (pow(d_z / (r * M), 2) / d_x) - 2 * (trials[1].z + trials[0].z - 2 * z_min) / (r * M);
		for (int i = 2; i < trials.size(); i++)
		{

			double k = trials[i].x - trials[i - 1].x;
			double R = k + (pow((trials[i].z - trials[i - 1].z) / (M * r), 2) / k) -
				2 * (trials[i].z + trials[i - 1].z - 2 * z_min) / (r * M);

			if (R > Rmax)
			{
				Rmax = R;
				Rpos = i;
			}
		}
		curr_eps = trials.at(Rpos).x - trials.at(Rpos - 1).x;

		std::vector<Trial>::iterator it2 = trials.begin();
		for (it = trials.begin(); it - trials.begin() <= Rpos; it++) it2 = it;

		int sgn = 0;
		if (trials[Rpos].z - trials[Rpos - 1].z < 0)
			sgn = -1;
		if (trials[Rpos].z - trials[Rpos - 1].z > 0)
			sgn = 1;

		current.x = (trials[Rpos].x + trials[Rpos - 1].x) / 2 - sgn * (trials[Rpos].z - trials[Rpos - 1].z) / (2 * M * r);
		current.z = Funk(key, index_problem, current.x);
		trials.insert(it2, current);

		int pos = find(trials.begin(), trials.end(), current) - trials.begin();
		if (optimum.z > current.z)
		{
			best_i = itr;
			optimum = current;
		}
		itr++;
	}
	std::cout << "itr = " << itr << std::endl;
}