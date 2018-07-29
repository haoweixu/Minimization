#ifndef __POWELL__
#define __POWELL__

#include "LineMethod.h"
#include "macros.h"
#include <vector>

template <class T>
class Powell :
	public LineMethod<T>
{
public:
	int iter;
	double fret;
	using LineMethod<T> ::func;
	using LineMethod<T> ::linmin;
	using LineMethod<T> ::p;
	using LineMethod<T> ::xi;

	const double ftol;
	Powell(T &func, const double _ftol = 3.0e-8) :
		LineMethod<T>(func), ftol(_ftol) {};

	vector<double> minimize(vector<double> &pp);
	vector<double> minimize(vector<double> &pp, Matrix &ximat);

	Powell() {};
	~Powell() {};
};

template <class T>
vector<double> Powell<T>::minimize(vector<double> &pp) {
	int n = pp.size();
	Matrix ximat(n, n);
	for (size_t i = 0; i < n; i++) {
		for (size_t j = 0; j < n; j++)
			ximat(i, j) = 0;
		ximat(i, i) = 1.0;
	}
	return minimize(pp, ximat);
}

template <class T>
vector<double> Powell<T>::minimize(vector<double> &pp, Matrix &ximat) {
	const int ITMAX = 200;
	const double TINY = 1.0e-25;
	double fptt;
	int n = pp.size();
	p = pp;
	vector<double> pt(n), ptt(n);
	xi.resize(n);
	fret = func(p);
	for (int j = 0; j < n; j++)
		pt[j] = p[j];
	for (iter = 0; ; ++iter) {
		double fp = fret;
		int ibig = 0;
		double del = 0.0;
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++)
				xi[j] = ximat(j,i);
			fptt = fret;
			fret = linmin();
			if (fptt - fret > del) {
				del = fptt - fret;
				ibig = i + 1;
			}
		}

		if (2.0*(fp - fret) <= ftol * (abs(fp) + abs(fret)) + TINY) {
			return p;
		}
		if (iter == ITMAX)
			throw("Powell exceeding maximum iterations");
		for (int j = 0; j < n; j++) {
			ptt[j] = 2.0*p[j] - pt[j];
			xi[j] = p[j] - pt[j];
			pt[j] = p[j];
		}
		fptt = func(ptt);
		if (fptt < fp) {
			double t = 2.0*(fp - 2.0*fret + fptt)*pow(fp - fret - del, 2) - del * pow(fp - fptt, 2);
			if (t < 0) {
				fret = linmin();
				for (int j = 0; j < n; j++) {
					ximat(j, ibig - 1) = ximat(j, n - 1);
					ximat(j, n - 1) = xi[j];
				}
			}
		}
	}
}

#endif // !__POWELL__


