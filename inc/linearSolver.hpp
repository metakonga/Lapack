#ifndef LINEAR_SOLVER_HPP
#define LINEAR_SOLVER_HPP

#include "lapacke.h"

class linearSolver
{
public:
	linearSolver(
		int _major = LAPACK_COL_MAJOR,
		int _n = 0,
		int _nrhs = 0,
		int _lda = 0,
		int _ldb = 0)
		: major(_major)
		, n(_n)
		, nrhs(_nrhs)
		, lda(_lda)
		, ldb(_ldb)
	{
		ipiv = new lapack_int[n];
	}
	~linearSolver()
	{
		delete[] ipiv; ipiv = NULL;
	}

	lapack_int solve(double* a, double* b)
	{
		lapack_int info;
		info = LAPACKE_dgesv(major, n, nrhs, a, lda, ipiv, b, ldb);
		return info;
	}

	void setMajor(int _major) { major = _major; }

	lapack_int major;
	lapack_int n;
	lapack_int nrhs;
	lapack_int lda;
	lapack_int ldb;
	lapack_int *ipiv;
};

#endif