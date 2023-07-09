#include "QuadraticCyclic.h"
#include <math.h>

Matrix QuadraticCyclic::X(3,3);
Vector QuadraticCyclic::Y(3);
Vector QuadraticCyclic::A(3);

QuadraticCyclic::QuadraticCyclic(int tag, double wt, double fac_y)
:CyclicModel(tag, -1), weightFactor(wt), facty(fac_y),
 qx1(0.0), qy1(0.0), qx2(0.0), qy2(0.0), qx3(0.0), qy3(0.0)
{

}

QuadraticCyclic::~QuadraticCyclic()
{
  // does nothing
}

int QuadraticCyclic::createFullCycleTask()
{
int res = this->CyclicModel::createFullCycleTask();
    res+= createTask();

return res;
}

int QuadraticCyclic::createHalfCycleTask()
{
int res = this->CyclicModel::createHalfCycleTask();
    res+= createTask();
    
return res;
}


CyclicModel *QuadraticCyclic::getCopy()
{
CyclicModel *newModel = new QuadraticCyclic(getTag(), weightFactor, facty);	
	return newModel;
}


int QuadraticCyclic::createTask()
{
	if(f_bgn*f_end < 0)
	{
		double k0 = k_init;
		double delx0 = f_bgn/(resFactor*k0);  //say +ive

		double qy = facty;
		double x1 = d_bgn,                   y1 = f_bgn;

//		double x2 = d_bgn + delx0*(qy-1),    y2 = qy*f_bgn;
//		double x3 = qx*(d_bgn-delx0),        y3 = 0.0;

//		double x2 = qx*(d_bgn - delx0),      y2 = 0.0;
//		double x3 = d_end;
//		double y3 = f_end;



		double delx2 = delx0*(1 - qy);   // +ive
		double x2 = d_bgn - delx2;       // +ive
		double y2 = qy*f_bgn;            // +ive

		double x0 = d_bgn - delx0;       // +ive (at least 0)
		double y0 = 0.0;  // always - don't change
		
		double R  = sqrt( (x2 - x0)*(x2 - x0) + (y2 - y0)*(y2 - y0));
		double delx_end = d_end - x0;    // -ive
        double dely_end = f_end - y0;    // -ive
		
		double R_end = sqrt(delx_end*delx_end + dely_end*dely_end);

		double delx3 = delx_end*R/R_end; // -ive
		double y3 = f_end*R/R_end;       // -ive
		double x3 = x0 + delx3;          // +ive
		
//        opserr << "R = " << R << ", Re = " << R_end << ", R/Re = " << R/R_end << endln;
//        opserr << "x1, y1 = " << x1 << ", " << y1 << endln;
//        opserr << "x2, y2 = " << x2 << ", " << y2 << endln;
//        opserr << "x0, y0 = " << x0 << ", " << y0 << endln;
//        opserr << "x3, y3 = " << x3 << ", " << y3 << endln;
//
//        opserr << *this;
//        opserr << "\a";

        qx1 = x1; qy1 = y1;
        qx2 = x2; qy2 = y2;
        qx3 = x3; qy3 = y3;
                                       
		solveQuad(x1, y1, x2, y2, x3, y3);
	}

	return 0;
}


int QuadraticCyclic::solveQuad(double x1, double y1, double x2,
                               double y2, double x3, double y3)
{

	X(0,0) = x1*x1; X(0,1) = x1; X(0,2) = 1.0;
	X(1,0) = x2*x2; X(1,1) = x2; X(1,2) = 1.0;
	X(2,0) = x3*x3; X(2,1) = x3; X(2,2) = 1.0;

	Y(0) = y1; Y(1) = y2; Y(2) = y3;

	A = Y/X;
	a = A(0);
	b = A(1);
	c = A(2);
	
	// opserr << A;

	return 0;
}

double QuadraticCyclic::getQuadFactor(double x1, double y1, double dx)
{
// double dx = (d_curr - d_hist)/2;
double x_nxt = x1 + dx;
double y_nxt = a*x_nxt*x_nxt + b*x_nxt + c;
double x_prv = x1 - dx;
double y_prv = a*x_prv*x_prv + b*x_prv + c;

//	opserr << "x1, y1 = " << x_nxt << ", " <<  y_nxt << endln;
//	opserr << "x2, y2 = " << x_prv << ", " <<  y_prv << endln;
	return (rationalize(x_prv, y_prv, x_nxt, y_nxt));

}


double QuadraticCyclic::getTaskFactor()
{
double tfactor;
//	// redundant - only for print
//	if(d_curr >= 0 && !initYieldPos)
//		return 1.0;
//	if(d_curr  < 0 && !initYieldNeg)
//		return 1.0;
//	// end redundant



	if(yielding /*&& fabs(d_curr) >= fabs(d_end) */)
//		return resFactor; // will eventually unload
		tfactor = cycFactor_hist;
    else
    {
    	if(f_bgn*f_end < 0)
    	{
    		// if(contains(0.0, f_bgn, f_curr))
			if(contains(qy1, qy3, f_curr))  
    			tfactor = getQuadFactor(d_curr, f_curr, (d_curr - d_hist)/2);
    		else
			{
				tfactor=rationalize(d_curr, f_curr, d_end, f_end);
				tfactor = weightFactor*tfactor + (1 - weightFactor)*resFactor;
		    }
    	}
    	else // half-cycle
		{
    		tfactor = rationalize(d_bgn, f_bgn, d_end, f_end);
            tfactor = weightFactor*tfactor + (1 - weightFactor)*resFactor;
		}
	}

//	opserr << "tfactor = " << tfactor << endln;
	return tfactor;
}

void QuadraticCyclic::Print (OPS_Stream &s, int flag)
{
	this->CyclicModel::Print (s, flag);
	s << "+QuadraticCyclic\n";
	s << "   taskFactor = " << cycFactor << endln;
	s << "   a=" << a <<", b=" << b <<", c=" << c << endln;
	s << "----------------------------------------"
      << "----------------------------------------"
	  << endln;
}



