#ifndef _MACROS_
#define _MACROS_

#include<vector>

using namespace std;

template<class T>
inline T SQR(const T a) { return a * a; }

template<class T>
inline const T &MAX(const T &a, const T &b)
{
	return b > a ? (b) : (a);
}

template<class T>
inline const T &MIN(const T &a, const T &b)
{
	return b < a ? (b) : (a);
}

template<class T>
inline T SIGN(const T &a, const T &b)
{
	return b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a);
}

template<class T>
inline void SWAP(T &a, T &b)
{
	T dum = a; a = b; b = dum;
}


class Matrix
{
public:
	Matrix();
	Matrix(size_t rows, size_t cols);
	double& operator()(size_t i, size_t j);
	double operator()(size_t i, size_t j) const;
	int ncols();
	int nrows();

private:
	size_t mRows;
	size_t mCols;
	std::vector<double> mData;
};



#endif // !_MACROS_
