#include "_2DParallelJacobi.hpp"
#include <functional>
#include "opt.hpp"

using namespace std;
using namespace ff;


_2DParallelJacobi::_2DParallelJacobi(double tol, 
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

inline int _2DParallelJacobi::getNumberOfWorkers(){ return numberOfWorkers; }

inline int _2DParallelJacobi::getChunk(){ return chunk; }


inline void _2DParallelJacobi::operator()(_2DPoissonEquation * eq){

	_2DGrid * grid = eq->getGrid();

	long int _n = grid->getXlen(),
			 _m = grid->getYlen();

	double * __restrict__ _unew = grid->getU(),
		   * __restrict__ _uold = (double *)malloc(sizeof(double) * _n *_m);

	double * _f = eq->getF();

	double hhxx = grid->getXspacing() * grid->getXspacing(),
		   hhyy = grid->getYspacing() * grid->getYspacing();

	double Error = 10 * getTolerance();

	long int CI = ci,
	     	 CJ = cj;

	ParallelForReduce<double> pfr(getNumberOfWorkers(), (getNumberOfWorkers() < ff_numCores()));

	auto Fsum = [](double& v, const double elem) { v += elem; };

	auto Fcopy = [&_uold,&_unew,_m](const long i) {
    	//double * __restrict__ __attribute__((aligned(ALIGNMENT))) __uold = _uold;
    	//double * __restrict__ __attribute__((aligned(ALIGNMENT))) __unew = _unew;
		double * __restrict__ __uold = (double *)__builtin_assume_aligned(_uold, ALIGNMENT);
    	double * __restrict__ __unew = (double *)__builtin_assume_aligned(_unew, ALIGNMENT);


    	#pragma ivdep
    	for (int j=0; j < _m; j++)
      		__uold[_m*i + j] = __unew[_m*i + j];
    };

    //one-dimensional partitioning by rows? yeah, by rows
    /*auto Freduce = [&_uold,&_unew,_f,_m,hhxx,hhyy](const long i, double& Error) {
    	//double * __restrict__ __uold = (double *)__builtin_assume_aligned(_uold, ALIGNMENT);
    	//double * __restrict__ __unew = (double *)__builtin_assume_aligned(_unew, ALIGNMENT);
    	double * __restrict__ __f    = (double *)__builtin_assume_aligned(_f, ALIGNMENT);

    	//double * __restrict__ __uold __attribute__((aligned(ALIGNMENT)))  = _uold;
    	//double * __restrict__ __unew __attribute__((aligned(ALIGNMENT)))  = _unew;
    	//double * __restrict__ __f    __attribute__((aligned(ALIGNMENT)))  = _f;

    	
    	double cxcy[VEC_SIZE] __attribute__((aligned(ALIGNMENT)));
    	double cx[VEC_SIZE]   __attribute__((aligned(ALIGNMENT)));
    	double cy[VEC_SIZE]   __attribute__((aligned(ALIGNMENT)));

    	#pragma GCC ivdep
    	for(int i = 0; i < VEC_SIZE; i++){
    		cxcy[i] = hhxx*hhyy;
    		cx[i] = hhxx;
    		cx[i] = hhyy;
    	}

    	const int align_offset = (ALIGNMENT - ((uintptr_t)_uold) % ALIGNMENT)/ sizeof(double);

    	for (int j=1; j< align_offset; j++){
    		_unew[i*_m + j] = 0.5 * ( hhxx*hhyy*_f[i*_m + j] - hhyy*(_uold[i*_m + j-1] + _uold[i*_m + j+1]) -
								  hhxx*(_uold[(i-1)*_m + j] + _uold[(i+1)*_m + j]) )/(hhyy + hhxx);
      		

      		//__unew[i*_m + j] = val;
			Error += (_uold[i*_m + j] - _unew[i*_m + j]) * (_uold[i*_m + j] - _unew[i*_m + j]);
    	}

    	double * __restrict__ __uold = (double *)__builtin_assume_aligned(_uold+align_offset, ALIGNMENT);
    	double * __restrict__ __unew = (double *)__builtin_assume_aligned(_unew+align_offset, ALIGNMENT);

    	const int vec_loop_trip = (_m - align_offset) / VEC_SIZE;

    	#pragma GCC ivdep
    	for (int j=1; j < vec_loop_trip-1; j++){
    		for(int k = 0; k < VEC_SIZE; k++)
    			__unew[i*_m + j] = 0.5 * ( cxcy[k]*_f[i*_m + j] - cy[k]*(__uold[i*_m + j-1] + __uold[i*_m + j+1]) -
							  	cx[k]*(__uold[(i-1)*_m + j] + __uold[(i+1)*_m + j]) )/(hhyy + hhxx);

    		Error += (__uold[i*_m + j] - __unew[i*_m + j]) * (__uold[i*_m + j] - __unew[i*_m + j]);

    	}

    	for (int j=align_offset+VEC_SIZE*vec_loop_trip; j<_m-1; j++){
    		__unew[i*_m + j] = 0.5 * ( hhxx*hhyy*_f[i*_m + j] - hhyy*(__uold[i*_m + j-1] + __uold[i*_m + j+1]) -
						  hhxx*(__uold[(i-1)*_m + j] + __uold[(i+1)*_m + j]) )/(hhyy + hhxx);

    		Error += (__uold[i*_m + j] - __unew[i*_m + j]) * (__uold[i*_m + j] - __unew[i*_m + j]);
    	}	

    	//#pragma GCC ivdep
    	//for (long j=1; j < _m-1; j++){
        //	__unew[i*_m + j] = 0.5 * ( hhxx*hhyy*_f[i*_m + j] - hhyy*(__uold[i*_m + j-1] + __uold[i*_m + j+1]) -
		//						  hhxx*(__uold[(i-1)*_m + j] + __uold[(i+1)*_m + j]) )/(hhyy + hhxx);
      		

      		//__unew[i*_m + j] = val;
		//	Error += (__uold[i*_m + j] - __unew[i*_m + j]) * (__uold[i*_m + j] - __unew[i*_m + j]);
		//}
			
    };*/

    	auto Freduce = [&_uold,&_unew,_f,_m, _n,hhxx,hhyy](const long i, double& Error){

    		double * __restrict__ __uold = (double *)__builtin_assume_aligned(_uold, ALIGNMENT);
	    	double * __restrict__ __unew = (double *)__builtin_assume_aligned(_unew, ALIGNMENT);
	    	double * __restrict__ __f    = (double *)__builtin_assume_aligned(_f, ALIGNMENT);

	    	#pragma ivdep
	    	for(int j = 1; j < _m-1; j++){
	    		const long double val = 0.5 * ( hhxx*hhyy*_f[i*_m + j] + hhyy*(__uold[i*_m + j-1] + __uold[i*_m + j+1]) +
									  hhxx*(__uold[(i-1)*_m + j] + __uold[(i+1)*_m + j]) )/(hhyy + hhxx);
	      		//const long double val = 0.25 * ( -hhxx*_f[i*_m + j]  - __uold[i*_m + j-1] -  __uold[i*_m + j+1] +
				//					  		   - __uold[(i-1)*_m + j] - __uold[(i+1)*_m + j] );
	      		__unew[i*_m + j] = val;
				Error += (__uold[i*_m + j] - __unew[i*_m + j]) * (__uold[i*_m + j] - __unew[i*_m + j]);
	    	}


    	};

        auto FreduceBlocked = [&_uold,&_unew,_f,_m, _n,hhxx,hhyy, CI, CJ](const long ii, double& Error) {

	    	double * __restrict__ __uold = (double *)__builtin_assume_aligned(_uold, ALIGNMENT);
	    	double * __restrict__ __unew = (double *)__builtin_assume_aligned(_unew, ALIGNMENT);
	    	double * __restrict__ __f    = (double *)__builtin_assume_aligned(_f, ALIGNMENT);

	    	#pragma ivdep
	    	for(long jj = 1; jj < _m; jj+=CJ){
	        	for(long i = ii; i < min(ii + CI, _n-1); i++){
	          		for (long j = jj; j < min(jj + CJ, _m-1); j++){

	        			const long double val = 0.5 * ( hhxx*hhyy*_f[i*_m + j] + hhyy*(__uold[i*_m + j-1] + __uold[i*_m + j+1]) +
									  hhxx*(__uold[(i-1)*_m + j] + __uold[(i+1)*_m + j]) )/(hhyy + hhxx);
	      		
	      				__unew[i*_m + j] = val;
						Error += (__uold[i*_m + j] - __unew[i*_m + j]) * (__uold[i*_m + j] - __unew[i*_m + j]);

			        }
				} 
			}	  
		};


	while(getActualNumberOfIterations() < getMaxNumberOfIterations() && Error > getTolerance()){

		Error = 0.0;

		pfr.parallel_for(0,_n-1,1, getChunk(), Fcopy, getNumberOfWorkers());

		pfr.parallel_reduce(Error, 0.0, 1, _n-2, 1, getChunk(), Freduce, Fsum, getNumberOfWorkers());

		// pfr.parallel_reduce(Error, 0.0, 1, _n-1, CI, getChunk(), FreduceBlocked, Fsum, getNumberOfWorkers());
		
		Error = sqrt(Error)/sqrt(_n*_m);

		// pfr.parallel_for(0,_n-1,1, getChunk(), Fcopy, getNumberOfWorkers());

		incrActualNumberOfIterations();

		// printf("%i %e\n", getActualNumberOfIterations(), Error);

	}

	//std::swap(_unew, _uold);
	grid->setError(Error);
	free(_uold);

}
