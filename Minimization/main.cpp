#include<iostream>
#include <cmath>
#include "macros.h"
#include "GoldenSearch.h"
#include "Brent.h"
#include "bracket.h"

using namespace std;

double func(double x) {
	//return pow(x, 3) + 5 * pow(x, 2) + 3 * x + 5;
	return cos(x);// +5 * sin(x);
}

int main() {
	Brent *gs = new Brent();
	double a = 1, b = 2;
	//cout << func(a);
	gs->bracket(a, b, func);
	double xmin = gs->minimize(func);
	cout << gs->ax << "\n" << gs->bx << "\n" << gs->cx << "\n" << xmin;
	return 0;
}
