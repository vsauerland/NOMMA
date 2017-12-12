#include "regression.hpp" 

double regression::pav( int N, int sign, double *td, double *xd, double *xr )
// regression by an isotonic function,
// using "pool adjacent violators (PAV) algorithm" [BBBB72]
//
// INPUT
// N: length of the time series that is to be fit
// sign: indicates if fit must be monotonically increasing (sign=1)
//       or monotonically decreasing (sign=-1)
// td: time series times
// xd: time series values
// 
// OUTPUT
// xr: resulting data-fit time series
// returns SSE between xd and xr
//
// [BBBB72] Barlow, R. E., Bartholomew, D. J., Bremmer, J. M., Brunk, H. D.
// Statistical Inference under order restrictions. The theory and
// applications of isotonic regression. John Wiley & Sons, 1972
{
	double *w, *xrw;
	w = ( double* )calloc( N, sizeof( double ) );
	xrw = ( double* )calloc( N, sizeof( double ) );
	// if monotonic decresing regression is requested,
	// apply monotonic increasing regression to negated data
	// (and negate the result in the end)
	if ( sign != -1 ) sign = 1;
	for ( int i = 0; i < N; i++ )
	{
		xr[ i ] = sign * xd[ i ]; // initial solution = data
		w[ i ] = 1; // initial block lengths = 1
		xrw[ i ] = xr[ i ]; // initial block means = data
	}
	// initialize counters:
	int n = N; // number of blocks with locally constant fit
	int m = 0; // number of steps that violate monotony property
	for ( int i = 0; i < n - 1; i++ ) if ( xrw[ i ] > xrw[ i + 1] ) m++;

	while ( m > 0 )
	{
		int i = 0;
		int j = 0;
		// merge "consecutive downstairs blocks" with mean values
		while ( j < n )
		{
			int k = 0;
			while ( j + k + 1 < n &&  xrw[ j + k ] > xrw[ j + k + 1 ] ) k++;
			double wsum = 0;
			double xrsum = 0;
			for ( int l = j; l <= j + k; l++ )
			{
				wsum += w[ l ];
				xrsum += xr[ l ];
			}
			w[ i ] = wsum;
			xr[ i ] = xrsum;
			i++;
			j = j + k + 1;
		}
		n = i;
		for ( int i = 0; i < n; i++ ) xrw[ i ] = xr[ i ] / w[ i ];
		m = 0;
		for ( int i = 0; i < n - 1; i++ ) if ( xrw[ i ] > xrw[ i + 1 ] ) m++;
	}
	// derive solution
	int k = 0;
	for ( int i = 0; i < n; i++ )
	{
		for ( int j = 0; j < w[ i ]; j++ )
		{
			xr[ k ] = sign * xrw[ i ];
			k++;
		}
	}
	free( w );
	free( xrw );
	double r = 0;
	for ( int i = 0; i < N; i++ )
	{
		r = r + ( xr[ i ] - xd[ i ] ) * ( xr[ i ] - xd[ i ] );
	}
    return( r );
}

void regression::pavKnots( vector< Knot > &knots, int *g, int *f, int *b, double *xd, bool inform )
// the "pool adjacent violators (PAV) algorithm" adaption to knots [YW09]
//
// ATTENTION: This implementation as base for the LPAV algorithm
// DOES NOT WORK if initial consecutive knots have equal knot solutions
//
// INPUT/OUTPUT
// knots: a vector of knots that is manipulated to be monotonically
//        increasing (considering knot heights)
//
// OUTPUT
// g:  number of blocks in the final vector of knots
// f:  for each start index of a block the start index of the following block
//     (starting with 0 it will finally contain N - g irrelevant values)
// b:  for each end index of a block the end index of the following block
//     (starting with N it will finally contain N - g irrelevant values)
// xd: given data (e.g. time series of measurements)
///
// [YW09] Yeganova, L., Wilbur, W. J.
// Isotonic Regression under Lipschitz Constraint.
// J. Optim. Theory Allp. 141:429-443. doi:10.1007/s10957-008-9477-0
{
	int N = knots.size();
	if ( inform ) printf( "PAV_Knots: N = %i\n", N );
	double *v, *h;
	v = ( double* )calloc( N, sizeof( double ) );
	h = ( double* )calloc( N, sizeof( double ) );
	// Step 1 of modified PAV algorithm in [YW09]:
	for ( int i = 0; i < N; i++ )
	{
		f[ i ] = i + 1;
		b[ i ] = i - 1;
		v[ i ] = knots[ i ].getX(); // the "block solutions"
		h[ i ] = knots[ i ].getHeight(); // the "block heights"
	}
	f[ N ] = N + 1;
	b[ N ] = N - 1;
	// Step 2 of adjusted PAV algorithm [YW09]:
	queue< int > q1;
	queue< int > q2;
	for ( int i = 0; i < N; i++ ) q1.push( i );
	int m, flag; // Step 3 of adjusted PAV algorithm [YW09]
	while ( !q1.empty() ) // Step 4 of modified PAV algorithm [YW09]
	{
		m = 0; // Step 5 of adjusted PAV algorithm [YW09]
		while ( !q1.empty() ) // Step 6 of adjusted PAV algorithm [YW09]
		{
			// Step 6A of adjusted PAV algorithm [YW09]:
			int k = q1.front();
			q1.pop();
			if ( inform ) printf( "PAV_Knots: got k = %i from queue\n", k );
			if ( m <= k ) // Step 6B of adjusted PAV algorithm [YW09]
			{
				flag = 0; // Step 6C of adjusted PAV algorithm [YW09]
				while ( f[ k ] != N && v[ k ] + h[ k ] > v[ f[ k ] ] ) // Step 6D
				{
					if ( inform ) printf( "PAV_Knots: pooling:\n" );
					// Step 6D.I, pool blocks beginning at k anf f[ k ]:
					int u = f[ f[ k ] ]; // Step 6D.II
					if ( inform ) printf( "PAV_Knots: first block starts with knot %i, second block ends with knot %i\n", k, u - 1 );
					int aB = knots[ k ].getA();
					int bB = knots[ u - 1 ].getB();
					int nH = bB - aB;
					if ( inform ) printf( "PAV_Knots: the corresponding time series indices are %i and %i, so we will obtain %i h-values\n", aB, bB, nH );
					double *hB;
					hB = ( double* )calloc( nH, sizeof( double ) ); 
					double hSum = 0;
					int ii = 0;
					for ( int l = 0; l < u - k; l++ )
					{
						int al = knots[ k + l ].getA();
						int bl = knots[ k + l ].getB();
						if ( inform ) printf( "PAV_Knots: knot %i has solution %f, starts at index %i and ends at index %i\n", k + l, knots[ k + l ].getX(), al, bl );
						for ( int i = al; i < bl; i++ )
						{
							hB[ ii ] = knots[ k + l ].getH()[ i - al ];
							if ( inform ) printf( "PAV_Knots: the %i-th h-value is, thus, %f\n", ii, hB[ ii  ] );
							ii++;
						}
						hSum = hSum + knots[ k + l ].getHeight();
						if ( l < u - k - 1 )
						{
							hB[ ii ] = 0;
							if ( inform ) printf( "PAV_Knots: the %i-th h-value is, thus, 0\n", ii );
							ii++;
						}
					}
					if ( inform ) printf( "PAV_Knots: update the block solution v[ %i ]\n", k );
					v[ k ] = 0;
					if ( inform ) printf( "PAV_Knots: adding data values" );
					for ( int i = aB; i <= bB; i++ )
					{
						if ( inform ) printf( " %f", xd[ i ] );
						v[ k ] = v[ k ] + xd[ i ];
					}
					if ( inform ) printf( "\n" );
					if ( inform ) printf( "PAV_Knots: initialized v[ %i ] with data sum %f\n", k, v[ k ] );
					if ( inform ) printf( "PAV_Knots: now, subtracting multiples of hB values" );
					for ( int i = aB + 1; i <= bB; i++ )
					{
						if ( inform ) printf( "  %i * %f,", bB - i + 1, hB[ i - aB - 1 ] );
						v[ k ] = v[ k ] - ( bB - i + 1 ) * hB[ i - aB - 1 ];
					}
					if ( inform ) printf( "\n" );
					if ( inform ) printf( "PAV_Knots: after reducing by multiples of heights v[ %i ] is %f\n", k, v[ k ] );
					v[ k ] = v[ k ] / ( bB - aB + 1 );
					if ( inform ) printf( "PAV_Knots: new block solution v[ %i ] is %f\n", k, v[ k ] );
					if ( inform ) printf( "PAV_Knots: will free hB\n" );
					free( hB );
					if ( inform ) printf( "PAV_Knots: freed hB\n" );
					h[ k ] = h[ k ] + h[ f[ k ] ]; // new block height
					// update f and b
					f[ k ] = u; // Step 6D.III of adjusted PAV algorithm [YW09]
					if ( u < N + 1 ) b[ u ] = k; // Step 6D.IV
					m = u; // Step 6D.V
					flag = 1; // Step 6D.VI
				}
			}
			if ( flag == 1 && b[ k ] >= 0 ) q2.push( b[ k ] ); // Step 6E
		}
		q1.swap( q2 ); // Step 7, interchange q1 and q2
	}
	free( h );
	// derive number of blocks and associated index pointers;
	*g = 0;
	m = 0;
	while ( m < N )
	{
		double he = 0; // partial sum of heights of "knots" in block g
		if ( inform ) printf( "PAV_Knots: %i-th block solution is %f\n", *g, v[ m ] );
		for ( int i = m; i < f[ m ]; i++ )
		{
			knots[ i ].setX( v[ m ] + he );
			he = he + knots[ i ].getHeight();
			if ( inform ) printf( "PAV_Knots: corresponding solution of knot %i is %f\n", i, knots[ i ].getX() );
		}
		*g = *g + 1;
		m = f[ m ];
	}
	free( v );
}

double regression::lpav( int N, int sign, double L, double *td, double *xd, double *xr, bool inform )
// regression by an isotonic function with limited steepness
// using "Lipschitz pool adjacent violators (LPAV) algorithm" [YW09]
//
// INPUT
// N: length of the time series that is to be fit
// sign: indicates if fit must be monotonically increasing (sign=1)
//       or monotonically decreasing (sign=-1)
// L: maximum allowed steepness (absolute value of the first derivative)
// td: time series times
// xd: time series values
// 
// OUTPUT
// xr: resulting data-fit time series
// returns SSE between xd and xr the
//
// USES pavKnots
//
// [YW09] Yeganova, L., Wilbur, W. J.
// Isotonic Regression under Lipschitz Constraint.
// J. Optim. Theory Allp. 141:429-443. doi:10.1007/s10957-008-9477-0
{
	int g, *F, *B;
	double *x;
	vector< Knot > knots;
	x = ( double* )calloc( N, sizeof( double ) );
	F = ( int* )calloc( N + 1, sizeof( int ) );
	B = ( int* )calloc( N + 1, sizeof( int ) );

	// initialize
	for ( int i = 0; i < N; i++ )
	{
		x[ i ] = xd[ i ];
		if ( sign == -1 ) x[ i ] = -x[ i ];
		Knot knot( i, i, x );
		knots.push_back( knot );
	}
	vector< Knot > knotsCopy( knots );
	pavKnots( knots, &g, F, B, x, false ); // Step 1 of LPAV algorithm [YW09]
	// Step 2 of LPAV algorithm [YW09] skipped
	// (we shrink the vector "knots" instead of using K)
	// Step 3 of LPAV algorithm [YW09]:
	queue< int > q;
	int j = 0;
	if ( inform ) printf( "LPAV: enqueuing block starts:" );
	for ( int i = 1; i < g; i++ )
	{
		j = F[ j ];
		if ( inform ) printf( " %i", j );
		q.push( j );
	}
	if ( inform ) printf( "\n" );
	int flag = 1;
	while ( flag == 1 ) // Step 4 of LPAV algorithm [YW09]
	{
		if ( inform ) printf( "LPAV: searching for violated Lipschitz constraints\n" );
		flag = 0; // Step 5 of LPAV algorithm [YW09]
		int nTied = 0; // not working with K we need to count tied knots
		while ( !q.empty() ) // Step 6 of LPAV algorithm [YW09]
		{
			// Step 6A of LPAV algorithm [YW09]:
			j = q.front() - nTied;
			q.pop();
			if ( inform ) printf( "LPAV: step 6A, j = %i - %i = %i\n", j + nTied, nTied, j );
			// Step 6B of LPAV algorithm [YW09]:
			Knot knotA = knots[ j - 1 ];
			Knot knotB = knots[ j ];
			double gap = knotB.getX() - knotA.getX() - knotA.getHeight();
			if ( knotA.getB() + 1 != knotB.getA() )
			{
				if ( inform ) printf( "LPAV: last index of knotA plus 1 should equal first index of knotB\n" );
				if ( inform ) printf( "LPAV: but the indices are %i and %i\n", knotA.getB(), knotB.getA() );
				exit( 1 );
			}
			double maxGap = L * ( td[ knotB.getA() ] - td[ knotA.getB() ] );
			if ( gap > maxGap )
			{
				if ( inform ) printf( "LPAV: tying knots at current junction [%i,%i]\n", j - 1, j );
				if ( inform ) printf( "LPAV: (since gap %f exceeds maxGap %f)\n", gap, maxGap );
				if ( inform )
				{
					printf( "LPAV: first knot:\n" );
					knotA.display();
					printf( "LPAV: second knot:\n" );
					knotB.display();
				}
				// Step 6B.I of LPAV algorithm [YW09] (tie knots):
				knotA.tie( knotB, maxGap, td, x );
				if ( inform )
				{
					printf( "LPAV: new knot:\n" );
					knotA.display();
				}
				// Step 6B.II of LPAV algorithm [YW09]
				// (shrink vector "knots" instead of using K):
				knots[ j - 1 ] = knotA;
				knots.erase( knots.begin() + j );
				nTied++;
				flag = 1; // Step 6B.III of LPAV algorithm [YW09]
			}
		} // while ( !q.empty() )
		if ( flag == 1 )
		{
// [YW09] mean to only keep the latest tied knots at their current solutions
// while all other knots shift back to the solutions they had
// before the last application of PAV:
// ----------------------------------------------------------------------------
			j = 0;
			for ( int i = 0; i < knots.size(); i++ )
			{
				while ( j < knotsCopy.size() && knotsCopy[ j ].getB() < knots[ i ].getB() ) j++;
				if ( knotsCopy[ j ].getB() != knots[ i ].getB() )
				{
					if ( inform ) printf( "LPAV: error while checking up untied knots\n" );
					exit( 1 );
				}
				if ( knotsCopy[ j ].getA() == knots[ i ].getA() )
				{
					knots[ i ].setX( knotsCopy[ j ].getX() );
				}
			}
// ----------------------------------------------------------------------------
			knotsCopy = knots;
			if ( inform )
			{
				printf( "LPAV: array of knots before next PAV_Knots:\n" );
				for ( int i = 0; i < knots.size(); i++ )
				{
					printf( "%i) ", i );
					knots[ i ].display();
				}
			}
			pavKnots( knots, &g, F, B, x, false ); // Step 7 of LPAV [YW09]
		}
		else if ( inform ) printf( "LPAV: no Lipschitz constraint violating knot, anymore\n" ); 
		// Step 8 of LPAV algorithm [YW09]:
		j = 0;
		if ( inform ) printf( "LPAV: enqueuing block starts:" );
		for ( int i = 1; i < g; i++ )
		{
			j = F[ j ];
			if( inform ) printf( " %i", j );
			q.push( j );
		}
		if( inform ) printf( "\n" );
	} // while ( flag == 1 ) (while fixed some violated Lipschitz constraint)

	// obtain regression time series from vector of knots:
	if ( inform ) printf( "LPAV: xr =" );
	for ( int i = 0; i < knots.size(); i++ )
	{
		Knot knot = knots[ i ];
		xr[ knot.getA() ] = knot.getX();
		if ( inform ) printf( " %f", xr[ knot.getA() ] );
		double sumH = 0;			for ( j = 1; j <= knot.getB() - knot.getA(); j++ )
		for ( j = 1; j <= knot.getB() - knot.getA(); j++ )
		{
			sumH = sumH + knot.getH()[ j - 1 ];
			xr[ knot.getA() + j ] = knot.getX() + sumH;
			if ( inform ) printf( " %f", xr[ knot.getA() + j ] );
		}
	}
	if ( inform ) printf( "\n" );
	if ( sign == -1 ) for ( int i = 0; i < N; i++ ) xr[ i ] = -xr[ i ];
	free( x );
	free( F );
	free( B );
	double r = 0;
	for ( int i = 0; i < N; i++ ) r = r + ( xr[ i ] - xd[ i ] ) * ( xr[ i ] - xd[ i ] );
	if ( inform ) printf( "LPAV: SSE = %f\n", r );
	return( r );
}

double regression::lpmrPartition( int N, int k, double dMin, double dMax, bool considerSteepness, double *td, double *xd, int *T )
// (nearly) optimum partition for piecewise monotonic regression
// with bounded steepness, using dynamic programming algorithm [DP91]
//
// note, that a solution w.r.t. the obtained partition might violate
// the slope limits in the "turning points" of the partition but provides a
// lower bound on the misfit of an actual solution to the problem
//
// INPUT
// N:    time series length
// k:    number of local extremes
// dMin: lower steepness bound (<0)
// dMax: upper steepness bound (>0)
// condiderSteepness: indicates if steepness bounds dMin, dMax are considered
//       true => (xr[i+1]-xr[i] )/(td[i+1]-td[i]) in [dMin,dMax]
// td:   time series times
// xd:   time series values
//
// OUTPUT
// T: k + 1 indices of the optimum partition of xd into monotonic segments
// returns SSE between xd and the resulting regression time series
//
// USES pav or lpav
//
// [DP91] I. C. Demetriou. M. J. D. Powell.
// Least Squares Smoothing of Univariate Data to achieve Piecewise Monotonicity.
// IMA Journal of Numerical Analysis 11:411-432 (1991)
// DOI 10.1093/imanum/11.3.411
{
	assert( k > 1 );
	double xnorm = 0;
	for ( int i = 0; i < N; i++ ) xnorm += xd[ i ] * xd[ i ];
	// STEP 0 of Algorithm 2 in [DP91]
	vector<int> U; // indices of local maxima in xd
	vector<int> L; // indices of local minima in xd
	int a = 0;
	while ( xd[ a + 1 ] == xd[ a ] && a < N - 1 ) a++;
	int b = N;
	while ( xd[ b - 1 ] == xd[ b ] && b > 0 ) b--;
	L.push_back( 0 );
	if ( xd[ a ] > xd[ a + 1 ] ) U.push_back( 0 );
	for ( int i = a + 1; i < b; i++ )
	{
	    int q = i;
    	while ( xd[ q + 1 ] == xd[ q ] && q < N - 2 ) q++;
    	if ( xd[ i - 1 ] > xd[ i ] && xd[ q ] < xd[ q + 1 ] ) L.push_back( i );
    	if ( xd[ i - 1 ] < xd[ i ] && xd[ q ] > xd[ q + 1 ] ) U.push_back( i );
	}
	if ( ( k % 2 ) == 0 || xd[ b - 1 ] > xd[ b ] ) L.push_back( N - 1 );
	if ( ( k % 2 ) == 1 || xd[ b - 1 ] < xd[ b ] ) U.push_back( N - 1 );
	// STEP 1 of Algorithm 2 in [DP91]
	// calculate best monotonically incr. fit to required data segments
	double *xr;
	xr = ( double* )calloc( N, sizeof( double ) );
	map< pair<int,int>, double > alpha;
    for ( int i = 0; i < L.size(); i++ )
	{
        int p = L[ i ];
        for ( int j = 0; j < U.size(); j++ )
		{
            int q = U[ j ];
            if ( p <= q )
			{
				double value;
				if ( !considerSteepness ) // original pmr algorithm [DP91] case
				{
					value = pav( q + 1 - p, 1, td + p, xd + p, xr + p );
				}
				else if ( p == 0 )
				{
					value = lpav( q + 1 - p, 1, dMax, td + p, xd + p, xr + p, false );
				}
				else
				{
					value = lpav( q - p, 1, dMax, td + p + 1, xd + p + 1, xr + p + 1, false );
				}
                alpha.insert( make_pair( make_pair( p, q ), value ) );
			}
        }
    }
    // calculate best monotonically decr. fit to required data segments
	map< pair<int,int>, double > beta;
    for ( int i = 0; i < U.size(); i++ )
	{
        int p = U[ i ];
        for ( int j = 0; j < L.size(); j++ )
		{
            int q = L[ j ];
            if ( p <= q )
			{
				double value;
				if ( !considerSteepness ) // original pmr algorithm [DP91] case
				{
					value = pav( q + 1 - p, -1, td + p, xd + p, xr + p );
				}
				else if ( p == 0 )
				{
					value = lpav( q + 1 - p, -1, -dMin, td + p, xd + p, xr + p, false );
				}
				else
				{
					value = lpav( q - p, -1, -dMin, td + p + 1, xd + p + 1, xr + p + 1, false );
				}
                beta.insert( make_pair( make_pair( p, q ), value ) );
			}
        }
    }
	free( xr );
	// STEP 2 of Algorithm 2 in [DP91]
	map< pair<int,int>, double > nu;
	map< pair<int,int>, int > tau;
	for ( int i = 0; i < U.size(); i++ )
	{
		double value = alpha.at( make_pair( 0, U[ i ] ) );
		nu.insert( make_pair( make_pair( U[ i ], 1 ), value ) );
	}
	// STEP 3 and STEP 4 of Algorithm 2 in [DP91]
	for ( int m = 2; m < k; m++ )
	{
		if ( m % 2 == 0 )
		{
			for ( int i = 0; i < L.size(); i++ )
			{
				int l = L[ i ];
				// formula(4.10):
				nu.insert( make_pair( make_pair( l, m ), xnorm ) );
				tau.insert( make_pair( make_pair( l, m ), l ) );
				for ( int j = 0; j < U.size(); j++ )
				{
					int t = U[ j ];
					if ( t <= l )
					{
						double value = nu.at( make_pair( t, m - 1 ) );
						value = value + beta.at( make_pair( t, l ) );
						if ( value < nu.at( make_pair( l, m ) ) )
						{
							nu[ make_pair( l, m ) ] = value;
							tau[ make_pair( l, m ) ] = t;
						}
					}
				}
			}
		}
		else
		{
			for ( int i = 0; i < U.size(); i++ )
			{
				int l = U[ i ];
				// formula(4.11):
				nu.insert( make_pair( make_pair( l, m ), xnorm ) );
				tau.insert( make_pair( make_pair( l, m ), l ) );
				for ( int j = 0; j < L.size(); j++ )
				{
					int t = L[ j ];
					if ( t <= l )
					{
						double value = nu.at( make_pair( t, m - 1 ) );
						value = value + alpha.at( make_pair( t, l ) );
						if( value < nu.at( make_pair( l, m ) ) )
						{
							nu[ make_pair( l, m ) ] = value;
							tau[ make_pair( l, m ) ] = t;
						}
					}
				}
			}
		}
	}
	// STEP 5 of Algorithm 2 in [DP91]
	nu.insert( make_pair( make_pair( N - 1, k ), xnorm ) );
	tau.insert( make_pair( make_pair( N - 1, k ), N - 1 ) );
	if ( k % 2 == 0 )
	{
		for ( int i = 0; i < U.size(); i++ )
		{
			int t = U[ i ];
			double value = nu.at( make_pair( t, k - 1 ) );
			value = value + beta.at( make_pair( t, N - 1 ) );
			if( value < nu.at( make_pair( N - 1, k ) ) )
			{
				nu[ make_pair( N - 1, k ) ] = value;
				tau[ make_pair( N - 1, k ) ] = t;
			}
		}
	}
	else
	{
		for ( int i = 0; i < L.size(); i++ )
		{
			int t = L[ i ];
			double value = nu.at( make_pair( t, k - 1 ) );
			value = value + alpha.at( make_pair( t, N - 1 ) );
			if( value < nu.at( make_pair( N - 1, k ) ) )
			{
				nu[ make_pair( N - 1, k ) ] = value;
				tau[ make_pair( N - 1, k ) ] = t;
			}
		}
	}
	// STEP 6 of Algorithm 2 in [DP91]
	T[ 0 ] = 0;
	T[ k ] = N - 1;
	double r = 0;
	for ( int i = k - 1; i > 0; i-- )
	{
		T[ i ] = tau.at( make_pair( T[ i + 1 ], i + 1 ) );
		if ( i % 2 == 0 ) r = r + alpha.at( make_pair( T[ i ], T[ i + 1 ] ) );
		else r = r + beta.at( make_pair( T[ i ], T[ i + 1 ] ) );
	}
	r = r + alpha.at( make_pair( 0, T[ 1 ] ) );
	return( r );
}

double regression::lpmr( int N, int k, int mode, double dMin, double dMax, bool considerSteepness, double *td, double *xd, double *xr )
// (nearly) optimum piecewise monotonic regression with bounded steepness
// using dynamic programming algorithm (lpmrPartition) [DP91]
//
// note, that an obtained solution might violate the slope limits
// in the "turning points" of the partition but provides a lower bound
// on the misfit of an actual solution to the problem
//
// INPUT
// N:    time series length
// k:    number of local extremes
// mode: indicates how to choose the first segment of regression time series
//        1: first segment is monotonically non-decreasing
//       -1: first segment is monotonically non-increasing
//        0: first segment is chosen w.r.t minimizing regression error
// dMin: lower steepness bound (<0)
// dMax: upper steepness bound (>0)
// condiderSteepness: indicates if steepness bounds dMin, dMax are considered
//       true => (xr[i+1]-xr[i] )/(td[i+1]-td[i]) in [dMin,dMax]
// td:   time series times
// xd:   time series values
//
// OUTPUT
// xr:   regression time series
// returns SSE between xd and xr
//
// USES lpmrPartition and (pav or lpav)
//
// [DP91] I. C. Demetriou. M. J. D. Powell.
// Least Squares Smoothing of Univariate Data to achieve Piecewise Monotonicity.
// IMA Journal of Numerical Analysis 11:411-432 (1991)
// DOI 10.1093/imanum/11.3.411
{
	double r, r1, r2;
	int p, q, *T, *T2;
	int sign = 1;
	T = ( int* )calloc( k + 1, sizeof( int ) );
	T2 = ( int* )calloc( k + 1, sizeof( int ) );
	if ( mode == 0 || mode == 1 )
	{
		// best piecewise monotonic regression
		// that starts with non-decreasing segment
		r1 = lpmrPartition( N, k, dMin, dMax, considerSteepness, td, xd, T );
	}
	if ( mode == 0 || mode == -1 )
	{
		// best piecewise monotonic regression
		// that starts with non-increasing segment
		for ( int i = 0; i < N; i++ ) xr[ i ] = -xd[ i ];
		r2 = lpmrPartition( N, k, dMin, dMax, considerSteepness, td, xr, T2 );
	}
	// choose partition of case 2, if better or desired:
	if ( mode == -1 || ( mode == 0 && r2 < r1 ) )
	{
		sign = -1;
		for ( int i = 0; i < k + 1; i++ ) T[ i ] = T2[ i ];
	}
	// obtain resulting time series and corresponding misfit:
	r = 0;
	p = T[ 0 ];
	q = T[ 1 ];
	if ( sign == 1 )
	{
		if ( !considerSteepness )
		{
			r = r + pav( q + 1 - p, 1, td + p, xd + p, xr + p );
		}
		else // ( considerSteepness )
		{
			r = r + lpav( q + 1 - p, 1, dMax, td + p, xd + p, xr + p, false );
		}
	}
	else // ( sign == -1 )
	{
		if ( !considerSteepness )
		{
			r = r + pav( q + 1 - p, -1, td + p, xd + p, xr + p );
		}
		else // ( considerSteepness )
		{
			r = r + lpav( q + 1 - p, -1, -dMin, td + p, xd + p, xr + p, false );
		}
	}
	sign = -sign;
	for ( int m = 1; m < k; m++ )
	{
		p = T[ m ];
		q = T[ m + 1 ];
		if ( sign == 1 )
		{
			if ( !considerSteepness )
			{
				r = r + pav( q - p, 1, td + p + 1, xd + p + 1, xr + p + 1 );
			}
			else // ( considerSteepness )
			{
				r = r + lpav( q - p, 1, dMax, td + p + 1, xd + p + 1, xr + p + 1, false );
			}
		}
		else // ( sign == -1 )
		{
			if ( !considerSteepness )
			{
				r = r + pav( q - p, -1, td + p + 1, xd + p + 1, xr + p + 1 );
			}
			else // ( considerSteepness )
			{
				r = r + lpav( q - p, -1, -dMin, td + p + 1, xd + p + 1, xr + p + 1, false );
			}
		}
		sign = -sign;
	}
	free( T );
	free( T2 );
	return( r );
}
