#include "_2DGrid.hpp"
#include "_2DDirichlet.hpp"
#include "_2DSequentialJacobi.hpp"
#include "_2DSequentialGaussSeidel.hpp"
#include "_2DSequentialGaussSeidelRedBlack.hpp"
#include "_2DParallelJacobi.hpp"
#include "_2DParallelGaussSeidelRedBlack.hpp"
#include <cmath>
#include <math.h>
#include <stdlib.h>
#include <ff/parallel_for.hpp>

using namespace ff;

int main(int argc, char *argv[])
{

	int n = atoi(argv[1]),
		m = atoi(argv[2]);

	double epsilon = atof(argv[3]);

	int nw = atoi(argv[4]),
		chunk = atoi(argv[5]);

	int ci = atoi(argv[6]),
		cj = atoi(argv[7]),
		it = atoi(argv[8]);

	/*
	*	nabla^2 u = - f
	*  u = sin(pi*x)*sin(pi*y)
	   f = 2*pi*pi*sin(pi*x)*sin(pi*y)
	*/
	//auto g = [](double x, double y)->double { return exp(x*y); };
	//auto e = [](double x, double y)->double { return exp(x*y); };
	//auto f = [](double x, double y)->double { return (x*x + y*y)*exp(x*y); };
	auto g = [](double x, double y) -> double
	{ return 0; };
	auto f = [](double x, double y) -> double
	{ return 2 * M_PI * M_PI * sin(M_PI * x) * sin(M_PI * y); };
	auto e = [](double x, double y) -> double
	{ return sin(M_PI * x) * sin(M_PI * y); };

	//auto g = [](double x, double y)->double { return exp(x)*sin(y) + 0.25*(x*x + y*y); };
	//auto e = [](double x, double y)->double { return exp(x)*sin(y) + 0.25*(x*x + y*y); };
	//auto f = [](double x, double y)->double { return 1; };

#ifdef SEQ_GAUSS

	//printf("%-10s%-10s%-10s%-15s%-10s%-10s%-10s%-10s%-20s%-20s%-20s\n", "n", "m", "#workers", "time (msec)", "epsilon", "chunk", "c_i", "c_j", "residual", "e", "#iters");

	_2DPoissonEquation *LaPlaceEquation1 = new _2DPoissonEquation(new _2DGrid(0., 1., 0., 1., n, m),
																  new _2DDirichlet(g), f, e);
	_2DSequentialGaussSeidel Gauss(epsilon, it);

	ffTime(START_TIME);
	Gauss(LaPlaceEquation1);
	ffTime(STOP_TIME);

	printf("\n");
	printf("error = %f\n", LaPlaceEquation1->getSolutionError());
	printf("wrt exact solution error = %f\n", LaPlaceEquation1->checkExactSolution(g));
	printf("no. iterations = %d\n", Gauss.getActualNumberOfIterations());
	printf("time (ms): %g\n", ffTime(GET_TIME));
	printf("\n");

#endif

#ifdef SEQ_RBGAUSS

	_2DPoissonEquation *LaPlaceEquation2 = new _2DPoissonEquation(new _2DGrid(0., 1., 0., 1., n, m),
																  new _2DDirichlet(g),
																  f, e);

	_2DSequentialGaussSeidelRedBlack GaussRB(epsilon, it);

	ffTime(START_TIME);
	GaussRB(LaPlaceEquation2);
	ffTime(STOP_TIME);

	printf("\n");
	printf("error = %f\n", LaPlaceEquation2->getSolutionError());
	printf("wrt exact solution error = %f\n", LaPlaceEquation2->checkExactSolution(g));
	printf("no. iterations = %d\n", GaussRB.getActualNumberOfIterations());
	printf("time (ms): %g\n", ffTime(GET_TIME));
	printf("\n");

#endif

#ifdef SEQ_JAC

	_2DPoissonEquation *LaPlaceEquation3 = new _2DPoissonEquation(new _2DGrid(0., 1., 0., 1., n, m),
																  new _2DDirichlet(g),
																  f, e);

	_2DSequentialJacobi Jacobi(epsilon, it);

	ffTime(START_TIME);
	Jacobi(LaPlaceEquation3);
	ffTime(STOP_TIME);

	printf("\n");
	printf("error = %f\n", LaPlaceEquation3->getSolutionError());
	printf("wrt exact solution error = %f\n", LaPlaceEquation3->checkExactSolution(g));
	printf("no. iterations = %d\n", Jacobi.getActualNumberOfIterations());
	printf("time (ms): %g\n", ffTime(GET_TIME));
	printf("\n");
	//LaPlaceEquation3->printU();

#endif

	//*********************** PARALLEL SOLVING *********************

#ifdef PAR_JAC
	_2DPoissonEquation *LaPlaceEquation4 = new _2DPoissonEquation(new _2DGrid(0., 1., 0., 1., n, m),
																  new _2DDirichlet(g),
																  f, e);
	_2DParallelJacobi parJacobi(epsilon, it, nw, chunk, ci, cj);

	ffTime(START_TIME);
	parJacobi(LaPlaceEquation4);
	ffTime(STOP_TIME);

	printf("\n");
	printf("error = %f\n", LaPlaceEquation4->getSolutionError());
	printf("wrt exact solution error = %f\n", LaPlaceEquation4->checkExactSolution(g));
	printf("no. iterations = %d\n", parJacobi.getActualNumberOfIterations());
	printf("time (ms): %g\n", ffTime(GET_TIME));
	printf("\n");

#endif

#ifdef PAR_GAUSS
	_2DPoissonEquation *LaPlaceEquation5 = new _2DPoissonEquation(new _2DGrid(0., 1., 0., 1., n, m),
																  new _2DDirichlet(g),
																  f, e);
	_2DParallelGaussSeidelRedBlack parGaussSeidel(epsilon, it, nw, chunk, ci, cj);

	ffTime(START_TIME);
	parGaussSeidel(LaPlaceEquation5);
	ffTime(STOP_TIME);

	printf("\n");
	printf("error = %f\n", LaPlaceEquation5->getSolutionError());
	printf("wrt exact solution error = %f\n", LaPlaceEquation5->checkExactSolution(g));
	printf("no. iterations = %d\n", parGaussSeidel.getActualNumberOfIterations());
	printf("time (ms): %g\n", ffTime(GET_TIME));

#endif
}