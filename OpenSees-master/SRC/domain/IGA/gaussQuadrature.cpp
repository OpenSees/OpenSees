/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 1999, The Regents of the University of California    **
** All Rights Reserved.                                               **
**                                                                    **
** Commercial use of this program without express permission of the   **
** University of California, Berkeley, is strictly prohibited.  See   **
** file 'COPYRIGHT'  in main directory for information on usage and   **
** redistribution,  and for a DISCLAIMER OF ALL WARRANTIES.           **
**                                                                    **
** Developed by:                                                      **
**   Frank McKenna (fmckenna@ce.berkeley.edu)                         **
**   Gregory L. Fenves (fenves@ce.berkeley.edu)                       **
**   Filip C. Filippou (filippou@ce.berkeley.edu)                     **
**                                                                    **
** ****************************************************************** */

// The function quadrature returns a n x 1 column vector W of quadrature
// weights and a n x dim matrix of quadrature points, where n is the
// number of quadrature points.  The function is called as follows:
//
// [W,Q]=quadrature( nint, type, dim )
//
// nint is the quadrature order, type is the type of quadrature
// (i.e. gaussian, triangular, etc.. ) and dim is the number of spacial
// dimentions of the problem.  The default for type is GAUSS and the
// default for dim is unity.
//
// wrQ=quadrature(nint,'TRIANGULAR',2);itten by Jack Chessa
//            j-chessa@northwestern.edu
// Department of Mechanical Engineering
// Northwestern University

// Adapted to OpenSees by Felipe Elgueta and José A. Abell (UANDES, Chile) www.joseabell.com
// Only using gaussian quadrature, spatial dimension 2

#include <IGASurfacePatch.h>
#include <elementAPI.h>
#include <ID.h>
#include <Vector.h>
#include <Matrix.h>
#include <Subdomain.h>
#include <OPS_Globals.h>

#include <math.h>

void gaussQuad(int order, Vector& pt, Vector& wt)
{



	int n = order;
	int x1 = 1;
	int x2 = -1;

	double eps = 1e-15;
	double m = (n + 1) / 2;
	double xm = (1 / 2) * (x2 + x1);
	double xl = (1 / 2) * (x2 - x1);

	// double pi = 3.1415926535897932846;
	double pi = acos(-1);

	double z1 = 0;
	double z;
	double p1;
	double p2;
	double p3;
	double pp;

	for (int i = 0; i < m; ++i)
	{
		z = cos(pi * (i + 1 - 0.25) / (n + 0.5));
		p1 = 1.0;
		p2 = 0.0;

		while (abs(z - z1) > eps)
		{
			p1 = 1.0;
			p2 = 0.0;
			for (int j = 0; j < n; ++j)
			{
				p3 = p2;
				p2 = p1;
				p1 = ((2 * (j + 1) - 1) * z * p2 - ((j + 1) - 1) * p3) / (j + 1);
			}
			pp = n * (z * p1 - p2) / (z * z - 1);
			z1 = z;
			z = z1 - p1 / pp;
		}
		pt(i) = xm - x1 * z;
		pt(n + 1 - (i + 1 + 1)) = xm + x1 * z;
		wt(i) = 2 * x1 / ((1 - z * z) * pp * pp);
		wt(n + 1 - (i + 1 + 1)) = wt(i);
	}



	if (order == 1) {
		pt(0) = 0.000000000000000;
		wt(0) = 2.000000000000000;

	}
	else if (order == 2) {
		pt(0) = 0.577350269189626;
		pt(1) = -0.577350269189626;

		wt(0) = 1.000000000000000;
		wt(1) = 1.000000000000000;

	}
	else if (order == 3) {
		pt(0) = 0.774596669241483;
		pt(1) = -0.774596669241483;
		pt(2) = 0.000000000000000;

		wt(0) = 0.555555555555556;
		wt(1) = 0.555555555555556;
		wt(2) = 0.888888888888889;

	}
	else if (order == 4) {
		pt(0) = 0.861134311594053;
		pt(1) = -0.861134311594053;
		pt(2) = 0.339981043584856;
		pt(3) = -0.339981043584856;

		wt(0) = 0.347854845137454;
		wt(1) = 0.347854845137454;
		wt(2) = 0.652145154862546;
		wt(3) = 0.652145154862546;

	}
	else if (order == 5) {
		pt(0) = 0.906179845938664;
		pt(1) = -0.906179845938664;
		pt(2) = 0.538469310105683;
		pt(3) = -0.538469310105683;
		pt(4) = 0.000000000000000;

		wt(0) = 0.236926885056189;
		wt(1) = 0.236926885056189;
		wt(2) = 0.478628670499366;
		wt(3) = 0.478628670499366;
		wt(4) = 0.568888888888889;

	}
	else if (order == 6) {
		pt(0) = 0.932469514203152;
		pt(1) = -0.932469514203152;
		pt(2) = 0.661209386466265;
		pt(3) = -0.661209386466265;
		pt(4) = 0.238619186003152;
		pt(5) = -0.238619186003152;

		wt(0) = 0.171324492379170;
		wt(1) = 0.171324492379170;
		wt(2) = 0.360761573048139;
		wt(3) = 0.360761573048139;
		wt(4) = 0.467913934572691;
		wt(5) = 0.467913934572691;

	}
	else if (order == 7) {
		pt(0) = 0.949107912342759;
		pt(1) = -0.949107912342759;
		pt(2) = 0.741531185599394;
		pt(3) = -0.741531185599394;
		pt(4) = 0.405845151377397;
		pt(5) = -0.405845151377397;
		pt(6) = 0.000000000000000;

		wt(0) = 0.129484966168870;
		wt(1) = 0.129484966168870;
		wt(2) = 0.279705391489277;
		wt(3) = 0.279705391489277;
		wt(4) = 0.381830050505119;
		wt(5) = 0.381830050505119;
		wt(6) = 0.417959183673469;

	}
	else if (order == 8) {
		pt(0) = 0.960289856497536;
		pt(1) = -0.960289856497536;
		pt(2) = 0.796666477413627;
		pt(3) = -0.796666477413627;
		pt(4) = 0.525532409916329;
		pt(5) = -0.525532409916329;
		pt(6) = 0.183434642495650;
		pt(7) = -0.183434642495650;

		wt(0) = 0.101228536290376;
		wt(1) = 0.101228536290376;
		wt(2) = 0.222381034453374;
		wt(3) = 0.222381034453374;
		wt(4) = 0.313706645877887;
		wt(5) = 0.313706645877887;
		wt(6) = 0.362683783378362;
		wt(7) = 0.362683783378362;
	}
	else if (order == 9) {
		pt(0) = -0.968160239507626;
		pt(1) = -0.836031107326636;
		pt(2) = -0.613371432700590;
		pt(3) = -0.324253423403809;
		pt(4) = 0.000000000000000;
		pt(5) = 0.324253423403809;
		pt(6) = 0.613371432700590;
		pt(7) = 0.836031107326636;
		pt(8) = 0.968160239507626;

		wt(0) = 0.081274388361574;
		wt(1) = 0.180648160694857;
		wt(2) = 0.260610696402935;
		wt(3) = 0.312347077040003;
		wt(4) = 0.330239355001260;
		wt(5) = 0.312347077040003;
		wt(6) = 0.261610696402935;
		wt(7) = 0.180648160694857;
		wt(8) = 0.081274388361574;
	}
	else if (order == 10) {
		pt(0) = -0.973906528517172;
		pt(1) = -0.865063366688985;
		pt(2) = -0.679409568299024;
		pt(3) = -0.433395394129247;
		pt(4) = -0.148874338981631;
		pt(5) = 0.148874338981631;
		pt(6) = 0.433395394129247;
		pt(7) = 0.679409568299024;
		pt(8) = 0.865063366688985;
		pt(9) = 0.973906528517172;

		wt(0) = 0.066671344308688;
		wt(1) = 0.149451349150581;
		wt(2) = 0.219086362515982;
		wt(3) = 0.269266719309996;
		wt(4) = 0.295524224714753;
		wt(5) = 0.295524224714753;
		wt(6) = 0.269266719309996;
		wt(7) = 0.219086362515982;
		wt(8) = 0.149451349150581;
		wt(9) = 0.066671344308688;
	}
	else if (order == 12) {
		pt(0) = -0.981560634246719;
		pt(1) = -0.904117256370475;
		pt(2) = -0.769902674194305;
		pt(3) = -0.587317954286617;
		pt(4) = -0.367831498998180;
		pt(5) = -0.125233408511469;
		pt(6) = 0.125233408511469;
		pt(7) = 0.367831498998180;
		pt(8) = 0.587317954286617;
		pt(9) = 0.769902674194305;
		pt(10) = 0.904117256370475;
		pt(11) = 0.981560634246719;

		wt(0) = 0.047175336386512;
		wt(1) = 0.106939325995318;
		wt(2) = 0.160078328543346;
		wt(3) = 0.203167426723066;
		wt(4) = 0.233492536538355;
		wt(5) = 0.249147045813403;
		wt(6) = 0.249147045813403;
		wt(7) = 0.233492536538355;
		wt(8) = 0.203167426723066;
		wt(9) = 0.160078328543346;
		wt(10) = 0.106939325995318;
		wt(11) = 0.047175336386512;
	}


}

void gaussQuad2dNurbs(int orderU, int orderV, Matrix* quadPoint, Vector* quadWeight)
{
	Vector ptU(orderU);
	Vector ptV(orderV);
	Vector wtU(orderU);
	Vector wtV(orderV);

	gaussQuad(orderU, ptU, wtU);
	gaussQuad(orderV, ptV, wtV);

	int n = 0;
	for (int i = 0; i < orderU; ++i)
	{
		for (int j = 0; j < orderV; ++j)
		{
			(*quadPoint)(n, 0) = ptU(i); 
			(*quadPoint)(n, 1) = ptV(j); 
			(*quadWeight)(n) = wtU(i) * wtV(j); 
			n++;
		}
	}
}

