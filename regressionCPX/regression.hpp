#include <ilcplex/ilocplex.h>
#include <string>
#include <iostream>
#include <fstream>
#include <queue>
#include <map>
#include "knot.hpp"

using namespace std;

class regression
{
/****************************** WITHOUT CPLEX ********************************/
	private:
	void pavKnots( vector< Knot > &knots, int *g, int *f, int *b, double *xd, bool inform );
	// regression by an isotonic series of "knots" (needed by lpav)
	double lpmrPartition( int N, int k, double dMin, double dMax, bool considerSteepness, double *td, double *xd, int *T );
	// optimum partition for piecewise monotonic regression with (optional) steepness bounds

	public:
	double pav( int N, int sign, double *xd, double *xr );
	// regression by an isotonic function
	double lpav( int N, int sign, double L, double *td, double *xd, double *xr, bool inform );
	// regression by an isotonic function with steepness bound
	double lpmr( int N, int k, int mode, double dMin, double dMax, bool considerSteepness, double *td, double *xd, double *xr );
	// piecewise monotonic regression with (optional) steepness bounds

/******************************* WITH CPLEX **********************************/
	double isoReg( int N, int sign, double *xd, double *xr );
	// regression by an isotonic function
	double slopeReg( int N, double dMin, double dMax, double *td, double *xd, double *xr );
	// regression by a function that has limited slope (generalization of isoReg)
	double minMaxReg( int N, int kMinA, int kMinB, int kMaxA, int kMaxB, double dMin, double dMax, double T, double *td, double *xd, double *xr );
	// regression by a (periodic) function that has limited slope and exactly one
	// local minimum and exactly one local maximum (per period)
	double lpmrIQP( int N, int nMin, int nMax, int sign, double dMin, double dMax, double T, double *td, double *xd, double *xr );
	// piecewise monotonic regression as IQP
};
