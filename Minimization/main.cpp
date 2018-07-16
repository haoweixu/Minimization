#include<iostream>
#include <cmath>
#include "macros.h"
#include "GoldenSearch.h"
#include "bracket.h"

using namespace std;

double func(double x) {
	//return pow(x, 3) + 5 * pow(x, 2) + 3 * x + 5;
	return cos(x) + 5 * sin(x);
}

int main() {
	GoldenSearch *gs = new GoldenSearch();
	double a = 20, b = 100;
	//cout << func(a);
	gs->bracket(a, b, func);
	double xmin = gs->minimize(func);
	cout << gs->ax << "\n" << gs->bx << "\n" << gs->cx << "\n" << xmin;
	return 0;

	aaaaaaa
}