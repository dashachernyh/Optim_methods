#pragma once

struct MaxR
{
	double R, r;
	int pos;

	MaxR(double _R = 0, double _r = 0, int _pos = 0)
		:R(_R), r(_r), pos(_pos) {}
};