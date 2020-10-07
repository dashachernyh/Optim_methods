#pragma once
#include<vector>
double Funk(double x)
{
	//return fabs(x - 1);
	return (x*x);
}

struct Trial
{
	double x, z;
	Trial& operator = (const Trial &tr)
	{
		x = tr.x;
		z = tr.z;
		return *this;
	}
	bool operator == (const Trial &tr)
	{
		if (x == tr.x && z == tr.z) return true;
		else
			return false;
	}
};

class Method
{

	std::vector<Trial> trials;
	Trial current;
	Trial first;
	Trial second;
public:
	Trial GetCurr() { return current; }
	int GetPos() { return (find(trials.begin(), trials.end(), current) - trials.begin()); }
	Method(double x_0, double x_n);
	void solve(double r, double eps);
};
Method::Method(double x_0, double x_n)
{
	first.x = x_0;
	first.z = Funk(x_0);
	second.x= x_n;
	second.z = Funk(x_n);
	trials.push_back(first);
	trials.push_back(second);
	
}
void Method:: solve(double r, double eps)
{
	double M,Rmax,Rpos=1;
	double curr_eps = second.x - first.x;
	std::vector<Trial>::iterator it=trials.begin();

	M = fabs((trials[1].z - trials[0].z) / (trials[1].x - trials[0].x));

	while(curr_eps > eps)
	{
		for (int i = 1; i < trials.size(); i++)
		{
			double max;
			max= fabs((trials.at(i).z - trials.at(i - 1).z) / (trials.at(i).x - trials.at(i - 1).x));
			if (M > max)
				M = max;
		}

		if (M == 0)
			M = 1;
		else
			M = r * M;

		Rmax = M * (trials[1].x-trials[0].x) + (pow((trials[1].z-trials[0].z), 2) / (M * (trials[1].x-trials[0].x))) - 2 * (trials[1].z+trials[0].z);

		for (int i=1; i<trials.size();i++)
		{
			
			double k = M * (trials.at(i).x - trials.at(i-1).x);
			double R = k + (pow((trials.at(i).z - trials.at(i - 1).z), 2) / k) - 2 * (trials.at(i).z + trials.at(i - 1).z);
			
			if (R > Rmax)
			{
				Rmax = R;
				Rpos = i;
			}
		}
		curr_eps = trials.at(Rpos).x - trials.at(Rpos - 1).x;

		std::vector<Trial>::iterator it2 = trials.begin();
		for (it = trials.begin(); it - trials.begin() <= Rpos; it++) it2 = it;

		current.x = (trials.at(Rpos).x + trials.at(Rpos - 1).x) / 2 - (trials.at(Rpos).z - trials.at(Rpos - 1).z) / (2 * M);
		current.z = Funk(current.x);


		trials.insert(it2,current);
		int pos = find(trials.begin(), trials.end(), current) - trials.begin();
		std::cout <<"pos="<< pos<<" x="<<trials.at(pos).x<<" z="<<trials.at(pos).z<<std::endl;
	
	}


}