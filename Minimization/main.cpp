#include<iostream>
#include <cmath>
#include "macros.h"
#include "GoldenSearch.h"
#include "Brent.h"
#include "bracket.h"
#include "Amoeba.h"

using namespace std;

double func(vector<double> p) {
	//return pow(x, 3) + 5 * pow(x, 2) + 3 * x + 5;
	double y = 0;
	for (int i = 0; i < p.size(); i++)
		y += cos(p[i])*pow(p[i], 2);
	return y;
}

int main() {
	Amoeba *gs = new Amoeba();
	vector<double> point = {0.0, 0.0};
	//double a, b;
	//cin >> a >> b;
	////cout << a << b;
	////cout << func(a);
	//gs->bracket(a, b, func);
	double tol = 1.0;
	vector<double> xmin = gs->minimize(point, 2.5, func);
	double fmin = gs->fmin;
	double fsave = 0;
	while (tol >= 1.0e-8) {
		fsave = fmin;
		xmin = gs->minimize(xmin, 2.5, func);
		fmin = gs->fmin;
		tol = abs(fmin - fsave);	
	}	
	for (int i=0; i<xmin.size(); i++)
		cout << xmin[i] << "\n";

	cout << fmin;
	return 0;
}
