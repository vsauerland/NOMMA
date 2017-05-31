#include "knot.hpp"
#include <stdlib.h>
#include <stdio.h>

// knot is an auxiliary construct introduced in [YW09] for their LPAV algorithm
// which solves the "monotonic regression problem s.t. Lipschitz constraints"
//
// [YW09] Yeganova, L., Wilbur, W. J.
// Isotonic Regression under Lipschitz Constraint.
// J. Optim. Theory Allp. 141:429-443. doi:10.1007/s10957-008-9477-0

Knot::Knot( int aF, int bF, double *xd )
{
	a = aF;
	b = bF;
	height = 0;
	for ( int i = a; i < b; i++ )
	{
		double he = xd[ i + 1 ] - xd[ i ];
		h.push_back( he );
		height = height + he;
	}
	// calculate knot solution according to (7) in [YW09] (all weights = 1)
	X = 0;
	for ( int i = a; i <= b; i++ ) X = X + xd[ i ];
	for ( int i = a + 1; i <= b; i++ ) X = X - ( b - i + 1 ) * h[ i ];
	X = X / ( b - a + 1 );
}

int Knot::getA()
{
	return( a );
}

int Knot::getB()
{
	return( b );
}

double Knot::getX()
{
	return( X );
}

void Knot::setX( double xF )
{
	X = xF;
}

vector< double > Knot::getH()
{
	return( h );
}

double Knot::getHeight()
{
	return( height );
}

void Knot::display()
{
	printf( "Knot attributes:\n" );
	printf( "a = %i\n", a );
	printf( "b = %i\n", b );
	printf( "h:\n" );
	for ( int i = 0; i < b - a; i++ ) printf( "%f\n", h[ i ] );
	printf( "knot height = %f\n", height );
	printf( "values obtained from solution X and heights h:\n" );
	double sumH = 0;
	printf( "%f\n", X );
	for ( int i = 0; i < b - a; i++ )
	{
		sumH = sumH + h[ i ];
		printf( "%f\n", X + sumH );
	}
	printf( "\n" );
}

void Knot::tie( Knot &K, double jH, double *td, double *xd )
// tying the knot to another knot K that must be consecutive ( b + 1 == K.getA() )
// K   the consecutive knot to be tied to
// jH  additional h-value at junction (to close the gap)
// td  times of the time series
//     (to be fit by "isotonic regression under Lipschitz constraint")
// xd  values of the time series
{
	if ( b + 1 != K.getA() )
	{
		printf( "Knot::tie: can only tie consecutive knots!s\n" );
		exit( 1 );
	}
	h.push_back( jH );
	height = height + jH;
	double he;
	for ( int i = K.getA(); i < K.getB(); i++ )
	{
		he = K.getH()[ i - K.getA() ];
		h.push_back( he );
		height = height + he;
	}
	b = K.getB();
	// calculate knot solution according to (7) in [YW09] (all weights = 1)
	X = 0;
	for ( int i = a; i <= b; i++ ) X = X + xd[ i ];
	for ( int i = a + 1; i <= b; i++ ) X = X - ( b + 1 - i ) * h[ i - a - 1 ];
	X = X / ( b + 1 - a );
}
