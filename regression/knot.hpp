#include <vector>

using namespace std;

class Knot
// Knot is an auxiliary construct introduced in [YW09] for their LPAV algorithm
// which solves the "monotonic regression problem with Lipschitz constraints"
//
// [YW09] Yeganova, L., Wilbur, W. J.
// Isotonic Regression under Lipschitz Constraint.
// J. Optim. Theory Allp. 141:429-443. doi:10.1007/s10957-008-9477-0
{
	private:
	int a;				 // first index covered by the knot
	int b;				 // last index covered by the knot
	double X;			 // knot solution
	vector< double > h;  // (b-a) single positive heights
	double height;       // total height of the knot

	public:
	Knot( int aF, int bF, double *xd );
	int getA();
	int getB();
	double getX();
	void setX( double xF );
	vector< double > getH();
	double getHeight();
	void display();
	void tie( Knot &K, double jH, double *td, double *xd ); // tying the knot with another knot K that must be consecutive ( b + 1 == K.getA() ) 
};
