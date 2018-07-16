#ifndef _Bracket_
#define _Bracket_

struct BracketMethod
{
public:
	double ax, bx, cx, fa, fb, fc;

	BracketMethod();
	~BracketMethod();

	template <typename T, typename Func>
	void bracket(T, T, Func);

	inline void shft2(double &, double &, const double);
	inline void shft3(double &a, double &b, double &c,
		const double d);
	inline void mov3(double &a, double &b, double &c,
		const double d, const double e, const double f);

};

template <typename T, typename Func>
void BracketMethod::bracket(T a, T b, Func func) {
	const double GOLD = 1, GLIMIT = 100.0, TINY = 1.e-20;
	ax = a; bx = b;
	double fu;
	fa = func(ax);
	fb = func(bx);
	if (fb > fa) {
		SWAP(ax, bx);
		SWAP(fa, fb);
	}

	cx = bx + GOLD * (bx - ax);
	fc = func(cx);

	while (fb > fc) {
		double r = (bx - ax)*(fb - fc);
		double q = (bx - cx)*(fb - fa);
		double u = bx - ((bx - cx)*q - (bx - ax)*r) /
			(2.0*SIGN(MAX(abs(q - r), TINY), q - r));
		double ulim = bx + GLIMIT * (cx - bx);
		if ((bx - u) * (u - cx) > 0.0) {
			fu = func(u);
			if (fu < fc) {
				ax = bx;
				bx = u;
				fa = fb;
				fb = fu;
				return;
			}
			else if (fu > fb) {
				cx = u;
				fc = fu;
				return;
			}
			u = cx + GOLD * (cx - bx);
			fu = func(u);
		}
		else if ((cx - u)*(u - ulim) > 0.0) {
			fu = func(u);
			if (fu < fc) {
				shft3(bx, cx, u, u + GOLD * (u - cx));
				shft3(fb, fc, fu, func(u));
			}
		}
		else if ((u - ulim)*(ulim - cx) >= 0.0) {
			u = ulim;
			fu = func(u);
		}
		else {
			u = cx + GOLD * (cx - bx);
			fu = func(u);
		}
		shft3(ax, bx, cx, u);
		shft3(fa, fb, fc, fu);
	}
}

inline void BracketMethod::shft2(double &a, double &b, const double c) {
	a = b;
	b = c;
}

inline void BracketMethod::shft3(double &a, double &b, double &c, const double d) {
	a = b;
	b = c;
	c = d;
}

inline void BracketMethod::mov3(double &a, double &b, double &c,
	const double d, const double e, const double f) {
	a = d;
	b = e;
	c = f;
}


#endif