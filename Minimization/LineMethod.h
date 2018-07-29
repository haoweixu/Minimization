#ifndef __LINEMETHOD__
#define __LINEMETHOD__

#include<vector>
#include"Brent.h"

template <class T>
class LineMethod
{
public:
	vector<double> p;
	vector<double> xi;
	T &func;
	int n;

	LineMethod(T &_func) : func(_func) {};
	~LineMethod() {};
	double linmin();

};

template <class T>
double LineMethod<T>::linmin() {
	double ax, xx, xmin;
	n = p.size();
	F1dim<T> f1dim(p, xi, func);
	ax = 0.0;
	xx = 1.0;
	Brent brent;
	brent.bracket(ax, xx, f1dim);
	xmin = brent.minimize(f1dim);
	//xmin = 0;
	for (int j = 0; j < n; j++) {
		xi[j] *= xmin;
		p[j] += xi[j];
	}
	return brent.fmin;
}

template <class T>
class F1dim {
public:

	const vector<double> &p;
	const vector<double> &xi;
	int n;
	T &func;
	vector<double> xt;

	F1dim(vector<double> & _p, vector<double> & _xi, T & _func) :
		p(_p), n(_p.size()), xi(_xi), func(_func), xt(n) {};
	double operator () (const double x) {
		for (int j = 0; j < n; j++) 
			xt[j] = p[j] + x * xi[j];
		return func(xt);
	}

};



#endif // !__LINEMETHOD__




