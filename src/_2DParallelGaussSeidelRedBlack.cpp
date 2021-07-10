#include "_2DParallelGaussSeidelRedBlack.hpp"
#include <functional>
#include "opt.hpp"
#include <stdio.h>

using namespace std;
using namespace ff;

_2DParallelGaussSeidelRedBlack::_2DParallelGaussSeidelRedBlack(double tol,
															   int it,
															   int nw,
															   int chk,
															   int i,
															   int j) : _2DIterativePoissonSolver(tol, it)
{

	numberOfWorkers = nw;
	chunk = chk;
	ci = i;
	cj = j;
}

inline int _2DParallelGaussSeidelRedBlack::getNumberOfWorkers() { return numberOfWorkers; }

inline int _2DParallelGaussSeidelRedBlack::getChunk() { return chunk; }

inline void _2DParallelGaussSeidelRedBlack::operator()(_2DPoissonEquation *eq)
{

	_2DGrid *grid = eq->getGrid();

	long int _n = grid->getXlen(),
			 _m = grid->getYlen();

	double *__restrict__ __attribute__((aligned(ALIGNMENT))) _unew = grid->getU(),
															 *__restrict__ _uold;
	//* __restrict__ __attribute__((aligned(ALIGNMENT))) _uold = (double *)malloc(sizeof(double) * _n *_m);

	posix_memalign((void **)&_uold, ALIGNMENT, _n * _m * sizeof(double));

	double *_f = eq->getF();

	double hhxx = grid->getXspacing() * grid->getXspacing(),
		   hhyy = grid->getYspacing() * grid->getYspacing();

	double Error = 10 * getTolerance();

	long int CI = ci,
			 CJ = cj;

	ParallelForReduce<double> pfr(getNumberOfWorkers(), (getNumberOfWorkers() < ff_numCores()));

	auto Fsum = [](double &v, const double elem)
	{ v += elem; };

	auto Fcopy = [&_uold, &_unew, _m](const long i)
	{
		double *__restrict__ __attribute__((aligned(ALIGNMENT))) __uold = _uold;
		double *__restrict__ __attribute__((aligned(ALIGNMENT))) __unew = _unew;

#pragma ivdep
		for (int j = 0; j < _m; j++)
			_uold[_m * i + j] = _unew[_m * i + j];
	};

	auto FreduceRedBlocked = [&_uold, &_unew, _f, _m, _n, hhxx, hhyy, CI, CJ](const long ii, double &Error)
	{
		double *__restrict__ __uold = (double *)__builtin_assume_aligned(_uold, ALIGNMENT);
		double *__restrict__ __unew = (double *)__builtin_assume_aligned(_unew, ALIGNMENT);
		double *__restrict__ __f = (double *)__builtin_assume_aligned(_f, ALIGNMENT);

#pragma ivdep
		for (long jj = 1; jj < _m; jj += CJ)
		{
			for (long i = ii; i < min(ii + CI, _n - 1); i++)
			{
				for (long j = jj; j < min(jj + CJ, _m - 1); j++)
				{
					if ((i + j) % 2 == 1)
					{

						const long double val = 0.5 * (-hhxx * hhyy * _f[i * _m + j] + hhyy * (__uold[i * _m + j - 1] + __uold[i * _m + j + 1]) + hhxx * (__uold[(i - 1) * _m + j] + __uold[(i + 1) * _m + j])) / (hhyy + hhxx);

						__unew[i * _m + j] = val;
						Error += (__uold[i * _m + j] - val) * (__uold[i * _m + j] - val);
					}
				}
			}
		}
	};

	auto FreduceBlackBlocked = [&_uold, &_unew, _f, _m, _n, hhxx, hhyy, CI, CJ](const long ii, double &Error)
	{
		double *__restrict__ __uold = (double *)__builtin_assume_aligned(_uold, ALIGNMENT);
		double *__restrict__ __unew = (double *)__builtin_assume_aligned(_unew, ALIGNMENT);
		double *__restrict__ __f = (double *)__builtin_assume_aligned(_f, ALIGNMENT);

#pragma ivdep
		for (long jj = 1; jj < _m; jj += CJ)
		{
			for (long i = ii; i < min(ii + CI, _n - 1); i++)
			{
				for (long j = jj; j < min(jj + CJ, _m - 1); j++)
				{
					if ((i + j) % 2 == 0)
					{

						const long double val = 0.5 * (-hhxx * hhyy * _f[i * _m + j] + hhyy * (__unew[i * _m + j - 1] + __unew[i * _m + j + 1]) + hhxx * (__unew[(i - 1) * _m + j] + __unew[(i + 1) * _m + j])) / (hhyy + hhxx);

						__unew[i * _m + j] = val;
						Error += (__uold[i * _m + j] - val) * (__uold[i * _m + j] - val);
					}
				}
			}
		}
	};

	//one-dimensional partitioning by rows? yeah, by rows
	auto FreduceRed = [&_uold, &_unew, _f, _m, hhxx, hhyy](const long i, double &Error)
	{
		double *__restrict__ __uold __attribute__((aligned(ALIGNMENT))) = _uold;
		double *__restrict__ __unew __attribute__((aligned(ALIGNMENT))) = _unew;
		double *__restrict__ __f __attribute__((aligned(ALIGNMENT))) = _f;

#pragma ivdep
		for (long j = 1; j < _m - 2; j++)
			if ((i + j) % 2 == 1)
			{
				const long double val = 0.5 * (hhxx * hhyy * _f[i * _m + j] + hhyy * (_uold[i * _m + j - 1] + _uold[i * _m + j + 1]) + hhxx * (_uold[(i - 1) * _m + j] + _uold[(i + 1) * _m + j])) / (hhyy + hhxx);

				_unew[i * _m + j] = val;
				Error += (_uold[i * _m + j] - val) * (_uold[i * _m + j] - val);
			}
	};

	auto FreduceBlack = [&_unew, _f, _m, hhxx, hhyy](const long i, double &Error)
	{
		//double * __restrict__ _uold = uold;
		//double * __restrict__ _u    = u;
		//double * __restrict__ _f    = f;

		//double * __restrict__ __uold __attribute__((aligned(ALIGNMENT)))  = _uold;
		double *__restrict__ __unew __attribute__((aligned(ALIGNMENT))) = _unew;
		double *__restrict__ __f __attribute__((aligned(ALIGNMENT))) = _f;

		//const long _n   = n-1;
		//const long mj   = m*j;
		//const long mjp1 = mj+m;
		//const long mjm1 = mj-m;

#pragma ivdep
		for (long j = 1; j < _m - 2; j++)
			if ((i + j) % 2 == 0)
			{
				const long double val = 0.5 * (hhxx * hhyy * _f[i * _m + j] + hhyy * (_unew[i * _m + j - 1] + _unew[i * _m + j + 1]) + hhxx * (_unew[(i - 1) * _m + j] + _unew[(i + 1) * _m + j])) / (hhxx + hhxx);

				Error += (_unew[i * _m + j] - val) * (_unew[i * _m + j] - val);
				_unew[i * _m + j] = val;
			}
	};

	while (getActualNumberOfIterations() < getMaxNumberOfIterations() && Error > getTolerance())
	{

		Error = 0.0;

		pfr.parallel_for(0, _n - 1, 1, getChunk(), Fcopy, getNumberOfWorkers());

		pfr.parallel_reduce(Error, 0.0, 1, _n - 2, 1, getChunk(), FreduceRed, Fsum, getNumberOfWorkers());
		//interessante osservare come l'errore vada resettato a zero
		pfr.parallel_reduce(Error, 0.0, 1, _n - 2, 1, getChunk(), FreduceBlack, Fsum, getNumberOfWorkers());

		//pfr.parallel_reduce(Error, 0.0, 1, _n-1, CI, getChunk(), FreduceRedBlocked, Fsum, getNumberOfWorkers());
		//pfr.parallel_reduce(Error, Error, 1,_n-1, CI, getChunk(), FreduceBlackBlocked, Fsum, getNumberOfWorkers());

		Error = sqrt(Error) / sqrt(_n * _m);

		//pfr.parallel_for(0, _n-1, 1, getChunk(), Fcopy, getNumberOfWorkers());

		incrActualNumberOfIterations();

		printf("%i %e\n", getActualNumberOfIterations(), Error);
	}

	//std::swap(_unew, _uold);
	grid->setError(Error);
	free(_uold);
}
