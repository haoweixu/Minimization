#ifndef __AMOEBA__
#define __AMOEBA__
#include "bracket.h"
#include "macros.h"
#include "vector"

using namespace std;

class Amoeba : public LineSearch
{
public:

	const double ftol;
	int nfunc; // number of function evaluations;
	int mpts;
	int ndim;
	double fmin;
	vector<double> y; // function values at the vertices of the simplex;
	Matrix p; // current simplex;

	Amoeba(double _ftol = 1.0e-8) : ftol(_ftol) {};
	~Amoeba() {};

	template <class Func>
	vector<double> minimize(vector<double> &point, double del, Func func);

	template <class Func>
	vector<double> minimize(vector<double> &point, vector<double> &dels, Func func);

	template <class Func>
	vector<double> minimize(Matrix &pp, Func func);

	template <class Func>
	double amotry(Matrix &p, vector<double> &y, vector<double> &psum, const int ihi, const double fac, Func func);

	inline void get_psum(Matrix &p, vector<double> &psum);

};

template <class Func>
vector<double> Amoeba::minimize(vector<double> &point, double del, Func func) {
	vector<double> dels(point.size(), del);
//	point[0] *= 2;
	return minimize(point, dels, func);
}

template <class Func>
vector<double> Amoeba::minimize(vector<double> &point, vector<double> &dels, Func func) {
	ndim = point.size();
	Matrix pp(ndim + 1, ndim);
	for (size_t i = 0; i < ndim + 1; i++) {
		for (size_t j = 0; j < ndim; j++) 
			pp(i, j) = point[j];
		if (i != 0) pp(i, i - 1) += dels[i - 1];
	}
	return minimize(pp, func);
}

template <class Func>
vector<double> Amoeba::minimize(Matrix &pp, Func func) {
	const int NMAX = 50000;
	const double TINY = 1.0e-10;
	int ihi, ilo, inhi;
	mpts = pp.nrows();
	ndim = pp.ncols();
	vector<double> psum(ndim), pmin(ndim), x(ndim);
	p = pp;
	y.resize(mpts);
	for (size_t i = 0; i < mpts; i++) {
		for (size_t j = 0; j < ndim; j++)
			x[j] = p(i, j);
		y[i] = func(x);
	}
	nfunc = 0;

	get_psum(p, psum);
	for (;;) {
		ilo = 0;
		ihi = y[0] > y[1] ? (inhi = 1, 0) : (inhi = 0, 1);
		// determine which point is the highest, next-highest, and lowest, by looping over the points in the simplex;
		for (rsize_t i = 0; i < mpts; i++) {
			if (y[i] <= y[ilo]) ilo = i;
			if (y[i] > y[ihi]) {
				inhi = ihi;
				ihi = i;
			}
			else if (y[i] > y[inhi] && i != ihi) inhi = i;
		}

		double rtol = 2.0*abs(y[ihi] - y[ilo]) / (abs(y[ihi] + y[ilo]) + TINY);
		if (rtol < ftol) {
			SWAP(y[0], y[ilo]);
			for (size_t i = 0; i < ndim; i++) {
				SWAP(p(0, i), p(ilo, i));
				pmin[i] = p(0, i);
			}
			fmin = y[0];
			return pmin;
		}

		if (nfunc >= NMAX) throw("NMAX exceeded");
		nfunc += 2;

		// new iteration. First extrapolate by a factor -1 through the face of simplx across from the high point, i.e., reflect the simplex from the high point.
		double ytry = amotry(p, y, psum, ihi, -1.0, func);
		if (ytry <= y[ilo])
			// reflection gives a better result than the best point, so try an additional extrapolation by a factor of 2.
			ytry = amotry(p, y, psum, ihi, 2.0, func);
		else if (ytry >= y[ihi]) {
			// reflection gives a worse result thant the best point, so look for an intermediate lower point, i.e., do a one-dimensional contraction.
			double ysave = y[ihi];
			ytry = amotry(p, y, psum, ihi, 0.5, func);
			if (ytry > ysave) { // can't seem to get rid of the high point.
				for (size_t i = 0; i < mpts; i++) {
					if (i != ilo) {
						for (size_t j = 0; j < ndim; j++)
							p(i, j) = psum[j] = 0.5*(p(i, j) + p(ilo, j));
						y[i] = func(psum);
					}
				}
				nfunc += ndim;
				get_psum(p, psum);
			}
		}
		else --nfunc;
	}
}

inline void Amoeba::get_psum(Matrix &p, vector<double> &psum) {
	for (size_t j = 0; j < ndim; j++) {
		double sum = 0.0;
		for (size_t i = 0; i < mpts; i++)
			sum += p(i, j);
		psum[j] = sum;
	}
}

template <class Func>
double Amoeba::amotry(Matrix &p, vector<double> &y, vector<double> &psum, const int ihi, const double fac, Func func) {
	vector<double> ptry(ndim);
	double fac1 = (1.0 - fac) / ndim;
	double fac2 = fac1 - fac;
	for (size_t j = 0; j < ndim; j++)
		ptry[j] = psum[j] * fac1 - p(ihi, j)*fac2;
	double ytry = func(ptry);
	if (ytry < y[ihi]) {
		y[ihi] = ytry;
		for (size_t j = 0; j < ndim; j++) {
			psum[j] += ptry[j] - p(ihi, j);
			p(ihi, j) = ptry[j];
		}
	}
	return ytry;
}


#endif // !__AMOEBA__




