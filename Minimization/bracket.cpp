#include <iostream>
#include <cmath>
#include "macros.h"
#include "bracket.h"

using namespace std;

LineSearch::LineSearch()
{
	ax = bx = cx = 0;
	fa = fb = fc = 0;
}

LineSearch::~LineSearch()
{
}
