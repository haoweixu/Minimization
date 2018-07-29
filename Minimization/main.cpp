#include<iostream>
#include <cmath>
#include "macros.h"
#include "GoldenSearch.h"
#include "Brent.h"
#include "bracket.h"
#include "Amoeba.h"
#include "Powell.h"

using namespace std;

double func(vector<double> p) {
	
	double y = 0;
	for (int i = 0; i < p.size(); i++)
		y += cos(p[i])*pow(p[i], 2);
	return y;
}

double func1(double x) {
	return pow(x, 3) + 5 * pow(x, 2) + 3 * x + 5;
}

int main() {
	Amoeba *gs = new Amoeba();
	vector<double> point = {4.0, -2.0};
	double tol = 1.0;
	Powell<double (vector<double>)> powell(func);
	point = powell.minimize(point);

	for (int i = 0; i < point.size(); i++)
		cout << point[i] << "\t";
	cout << "\n" << powell.fret;



	/*vector<double> xmin = gs->minimize(point, 2.5, func);
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
	cout << fmin;*/

	/*double a, b;
	cin >> a >> b;
	LineSearch *gs;

	int flag;
	cin >> flag;

	if (flag == 1)
		gs = new Brent();
	else if (flag == 2)
		gs = new GoldenSearch();
	else
		throw("undefined linesearch method.");

	gs->bracket(a, b, func1);
	GS->minimize(func1);
	cout << GS->fmin;*/

	return 0;
}
