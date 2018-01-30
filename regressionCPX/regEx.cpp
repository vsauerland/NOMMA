#include "regression.hpp"

using namespace std;


int main( int argc, char **argv )
{
	string instanceName;
	instanceName = "sineLikeNoise_0.45_200.dat";
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
	ifstream inFile( instanceName.c_str() );
	while ( getline( inFile, line ) )
	{
		lines.push_back( line );
	}
	inFile.close(); 
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
	int n4 = n / 4;
	double L = 1.44; // Lipschitz constant / steepness bound
	// sine model derivation limit
	// according to f^(j) = b * d^j * g(dx-c), g in {-sin,sin,-cos,cos},
	// assuming |b| <= 1.2 and |d| <= 1.2, i.e. |b*d^j| <= 1.2^(j+1):

	printf( "Fitting synthetic measurements y = f(x) = sin(x) + 0.3*sin(2*x) + normal(0;0.45)\n\n" );
	// (sigma of normal distibution is 0.2*range(f(x)) = 0.2*2.273 ~ 0.45

	printf( "-------------------------------------------------------------------------------\n\n" );

	regression reg;
	ofstream ofFile;

	r = reg.pav( n4, 1, xd, xr );
	printf( "PAV (mon. incr.) for 1st quarter: SSE = %f, RMSE = %f\n", r, sqrt( r / n4 ) );
	ofFile.open( "regression_PAV_Q1.txt" );
	for ( int i = 0; i < n; i++ ) ofFile << td[ i ] << " " << xr[ i ] << "\n";
	ofFile.close();

	r = reg.isoReg( n4, 1, xd, xr );
	printf( "isoReg (mon. incr.) for 1st quarter: SSE = %f, RMSE = %f\n\n", r, sqrt( r / n4 ) );

	r = reg.pav( n4, -1, xd + n4, xr );
	printf( "PAV (mon. decr.) for 2nd quarter: SSE = %f, RMSE = %f\n\n", r, sqrt( r / n4 ) );
	ofFile.open( "regression_PAV_Q2.txt" );
	for ( int i = 0; i < n; i++ ) ofFile << td[ i ] << " " << xr[ i ] << "\n";
	ofFile.close();

	r = reg.lpav( n4, 1, L, td, xd, xr, false );
	printf( "LPAV (mon. incr., L=%f) for 1st quarter: SSE = %f, RMSE = %f\n", L, r, sqrt( r / n4 ) );
	ofFile.open( "regression_LPAV_Q1.txt" );
	for ( int i = 0; i < n; i++ ) ofFile << td[ i ] << " " << xr[ i ] << "\n";
	ofFile.close();

	r = reg.slopeReg( n4, 0, L, td, xd, xr );
	printf( "slopeReg (mon. incr., L=%f) for 1st quarter: SSE = %f, RMSE = %f\n\n", L, r, sqrt( r / n4 ) );

	r = reg.lpav( n4, -1, L, td + n4, xd + n4, xr, false );
	printf( "LPAV (mon. decr., L=%f) for 2nd quarter: SSE = %f, RMSE = %f\n\n", L, r, sqrt( r / n4 ) );
	ofFile.open( "regression_LPAV_Q2.txt" );
	for ( int i = 0; i < n; i++ ) ofFile << td[ i ] << " " << xr[ i ] << "\n";
	ofFile.close();

	printf( "-------------------------------------------------------------------------------\n\n" );

	int k = 3; // number of pieces for piece-wise monotonic regression
	int mode = 0; // 0: optimal, 1: begin increasing, -1: begin decreasing
	bool considerL = false; // do not apply steepness bounds

	r = reg.slopeReg( n, -L, L, td, xd, xr );
	printf( "slopeReg (steepness in [-L,L]): SSE = %f, RMSE = %f\n\n", r, sqrt( r / n ) );
	ofFile.open( "regression_Bounded_Steepness.txt" );
	for ( int i = 0; i < n; i++ ) ofFile << td[ i ] << " " << xr[ i ] << "\n";
	ofFile.close();

	r = reg.lpmr( n, k, mode, -L, L, considerL, td, xd, xr );
	printf( "PMR with %i \"extremes\": SSE = %f, RMSE = %f\n\n", k - 1, r, sqrt( r / n ) );
	ofFile.open( "regression_PMR.txt" );
	for ( int i = 0; i < n; i++ ) ofFile << td[ i ] << " " << xr[ i ] << "\n";
	ofFile.close();

	considerL = true; // apply steepness bounds, now
	r = reg.lpmr( n, k, mode, -L, L, considerL, td, xd, xr );
	printf( "LPMR with %i \"extremes\": SSE = %f, RMSE = %f\n\n", k - 1, r, sqrt( r / n ) );
	ofFile.open( "regression_LPMR_2EX.txt" );
	for ( int i = 0; i < n; i++ ) ofFile << td[ i ] << " " << xr[ i ] << "\n";
	ofFile.close();

	printf( "-------------------------------------------------------------------------------\n\n" );

	k = 5; // new number of monotonic pieces
	r = reg.lpmr( n, k, mode, -L, L, considerL, td, xd, xr );
	printf( "LPMR with %i \"extremes\": SSE = %f, RMSE = %f\n\n", k - 1, r, sqrt( r / n ) );
	ofFile.open( "regression_LPMR_4EX.txt" );
	for ( int i = 0; i < n; i++ ) ofFile << td[ i ] << " " << xr[ i ] << "\n";
	ofFile.close();

	double T = 6.2831853;
	r = reg.lpmrIQP( n, 1, 1, 0, -L, L, T, td, xd, xr );
	printf( "LPMR_IQP with 2 extremes: SSE = %f, RMSE = %f\n\n", r, sqrt( r / n ) );
	ofFile.open( "regression_LPMR_IQP_2EX.txt" );
	for ( int i = 0; i < n; i++ ) ofFile << td[ i ] << " " << xr[ i ] << "\n";
	ofFile.close();

	r = reg.lpmrIQP( n, 2, 2, 0, -L, L, T, td, xd, xr );
	printf( "LPMR_IQP with 4 extremes: SSE = %f, RMSE = %f\n\n", r, sqrt( r / n ) );
	ofFile.open( "regression_LPMR_IQP_4EX.txt" );
	for ( int i = 0; i < n; i++ ) ofFile << td[ i ] << " " << xr[ i ] << "\n";
	ofFile.close();

	printf( "-------------------------------------------------------------------------------\n\n" );

	printf( "Note that LPMR solutions might violate steepness bounds in turning points,\n" );
	printf( "while LPMR_IQP solutions satisfy steepness bounds everywhere.\n\n" );

	free( td );
	free( xd );
	free( xr );
}
