#pragma once
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
	bool operator < (const Trial& tr) // для функции sort
	{
		if (x < tr.x)return true;
		else
			return false;
	}

	friend Trial min(const Trial& tr, const Trial& tr_other) {
		if (tr.z <= tr_other.z)
			return tr;
		else
			return tr_other;
	}
};
