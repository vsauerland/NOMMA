#include "regression.hpp"

using namespace std;


int main( int argc, char **argv )
{
	string instanceName;
//	instanceName = "sineLikeNoise_0.2_200.dat";
	instanceName = "../regression/SMHI_BY5_phyt_seas_2columns.dat"; // VSA: nur kurz
	// instanceName = argv[ 1 ];
	double *td;
	double *xd;
	double *xr;
	double dMin, dMax;
	int n;

	// READ DATA 
	string line;
	vector<string> lines;
	stringstream ss;
	ifstream file( instanceName.c_str() );
	while ( getline( file, line ) )
	{
		lines.push_back( line );
	}
	file.close(); 
	n = lines.size();
	printf( "the data file has %i lines\n", n );
	td = ( double * )calloc( n, sizeof( double ) );
	xd = ( double * )calloc( n, sizeof( double ) );
	xr = ( double * )calloc( n, sizeof( double ) );
	for ( int i = 0; i < n; i++ )
	{
		ss.str( "" );
		ss.clear();
		ss << lines[ i ];
		ss >> td[ i ];
		ss >> xd[ i ];
	}

	// DETERMINE AND DISPLAY DATA STEEPNESS LIMITS
	dMin = ( xd[ 1 ] - xd[ 0 ] ) / ( td[ 1 ] - td[ 0 ] );
	dMax = dMin;
	for ( int i = 0; i < n; i++ )
	{
		double d = ( xd[ i + 1 ] - xd[ i ] ) / ( td[ i + 1 ] - td[ i ] );
		if ( d < dMin )
		{
			dMin = d;
		}
		if ( d > dMax )
		{
			dMax = d;
		}
	}
	printf( "data slope is between %f and %f\n\n", dMin, dMax );

	// CALCULATE REGRESSIONS
	double r;
	int n1 = n / 4;
//	double L = 1.44; // Lipschitz constant / steepness bound
	double L = 0.14; // VSA: nur kurz Lipschitz constant / steepness bound
	// sine model derivation limit
	// according to f^(j) = b * d^j * g(dx-c), g in {-sin,sin,-cos,cos},
	// assuming |b| <= 1.2 and |d| <= 1.2, i.e. |b*d^j| <= 1.2^(j+1):

	printf( "Fitting synthetic measurements y = f(x) = sin(x) + 0.3*sin(2*x) + normal(0;0.45)\n" );
	// (sigma of normal distibution is 0.2*range(f(x)) = 0.2*2.273 ~ 0.45

	r = pav( n1, 1, td, xd, xr );
	printf( "PAV (mon. incr.) for 1st quarter: SSE = %f, RMSE = %f\n", r, sqrt( r / n1 ) );
	r = isoReg( n1, 1, td, xd, xr );
	printf( "isoReg (mon. incr.) for 1st quarter: SSE = %f, RMSE = %f\n\n", r, sqrt( r / n1 ) );

	r = pav( n1, -1, td + n1, xd + n1, xr );
	printf( "PAV (mon. decr.) for 2nd quarter: SSE = %f, RMSE = %f\n\n", r, sqrt( r / n1 ) );

	r = lpav( n1, 1, L, td, xd, xr, false );
	printf( "LPAV (mon. incr., L=%f) for 1st quarter: SSE = %f, RMSE = %f\n", L, r, sqrt( r / n1 ) );

	r = slopeReg( 0, n1 - 1, 0, L, td, xd, xr );
	printf( "slopeReg (mon. incr., L=%f) for 1st quarter: SSE = %f, RMSE = %f\n\n", L, r, sqrt( r / n1 ) );

	r = lpav( n1, -1, L, td + n1, xd + n1, xr, false );
	printf( "LPAV (mon. decr., L=%f) for 2nd quarter: SSE = %f, RMSE = %f\n\n", L, r, sqrt( r / n1 ) );

//	double T = 6.2831853;
	double T = 365.0; // VSA: nur kurz
	int k = 3; // number of pieces for piece-wise monotonic regression
//	int mode = 0; // 0: optimal, 1: begin increasing, -1: begin decreasing
	int mode = 1; // VSA: nur kurz
	bool considerL = false; // do not apply steepness bounds
	r = lpmr( n, k, mode, -L, L, considerL, td, xd, xr );
	printf( "PMR with %i \"extremes\": SSE = %f, RMSE = %f\n\n", k - 1, r, sqrt( r / n ) );

	considerL = true; // apply steepness bounds
	r = lpmr( n, k, mode, -L, L, considerL, td, xd, xr );
	printf( "LPMR with %i \"extremes\": SSE = %f, RMSE = %f\n\n", k - 1, r, sqrt( r / n ) );
/* VSA: nur kurz
	// WRITE REGRESSION TIME SERIES TO FILE
	ofstream f2( "regressionLPMR.txt" );
	for ( int i = 0; i < n; i++ )
	{
		f2 << td[ i ] << " " << xr[ i ] << "\n";
	}
	f2.close();*/
	printf( "LPMR regression output written to \"regressionLPMR.txt\" \n" );

	k = 5;
	r = lpmr( n, k, mode, -L, L, considerL, td, xd, xr );
	printf( "LPMR with %i \"extremes\": SSE = %f, RMSE = %f\n\n", k - 1, r, sqrt( r / n ) );

	r = lpmrIQP( n, 1, 1, 0, -L, L, 0, td, xd, xr );
	printf( "LPMR_IQP with 2 extremes: SSE = %f, RMSE = %f\n\n", r, sqrt( r / n ) );
	
	r = lpmrIQP( n, 1, 1, 0, -L, L, T, td, xd, xr );
	printf( "LPMR_IQP with 2 extremes: SSE = %f, RMSE = %f\n\n", r, sqrt( r / n ) );

	r = lpmrIQP( n, 2, 2, 0, -L, L, T, td, xd, xr );
	printf( "LPMR_IQP with 4 extremes: SSE = %f, RMSE = %f\n\n", r, sqrt( r / n ) );

	// VSA: nur kurz
	r = lpmrIQP( n, 2, 2, -1, -L, L, T, td, xd, xr );
	printf( "LPMR_IQP with 4 extremes: SSE = %f, RMSE = %f\n\n", r, sqrt( r / n ) );
	ofstream f2( "regressionLPMR.txt" );
	for ( int i = 0; i < n; i++ )
	{
		f2 << td[ i ] << " " << xr[ i ] << "\n";
	}
	f2.close();

	free( td );
	free( xd );
	free( xr );
}
