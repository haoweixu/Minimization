#ifndef _GoldenSearch_
#define _GoldenSearch_

#include "bracket.h"

struct GoldenSearch : public BracketMethod
{
public:

	double xmin, fmin;
	double tol;

	GoldenSearch(double _tol = 1e-8);
	~GoldenSearch();

	template <class Func>
	double minimize(Func & func);
	
};


template<class Func>
double GoldenSearch::minimize(Func & func)
{
	const double R = 0.61803399, C = 1.0 - R;
	double x1, x2;
	double x0 = ax;
	double x3 = cx;
	if (abs(cx - bx) > abs(bx - ax)) {
		x1 = bx;
		x2 = bx + C * (cx - bx);
	}
	else {
		x2 = bx;
		x1 = bx - C * (bx - ax);
	}

	double f1 = func(x1);
	double f2 = func(x2);

	while (abs(x3 - x0) > tol*(abs(x1) + abs(x2))) {
		if (f2 < f1) {
			shft3(x0, x1, x2, R*x2 + C * x3);
			shft2(f1, f2, func(x2));
		}
		else {
			shft3(x3, x2, x1, R*x1 + C * x0);
			shft2(f2, f1, func(x1));
		}
	}
	if (f1<f2) {
		xmin = x1;
		fmin = f1;
	}
	else {
		xmin = x2;
		fmin = f2;
	}
	return xmin;
}


#endif // !_GoldenSearch_

