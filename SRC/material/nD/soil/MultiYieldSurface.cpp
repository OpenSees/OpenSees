// $Revision: 1.8 $
// $Date: 2003-02-14 23:01:30 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/nD/soil/MultiYieldSurface.cpp,v $
                                                                        
// Written: ZHY
// Created: August 2000
//
// MultiYieldSurface.cpp
// ---------------------
//

#include <math.h>
#include <stdlib.h>

#include <MultiYieldSurface.h>


// YieldSurface class methods
MultiYieldSurface::MultiYieldSurface():
theSize(0.0), theCenter(6), plastShearModulus(0.0)
{

}

MultiYieldSurface::MultiYieldSurface(const Vector & theCenter_init, 
                                     double theSize_init, double plas_modul):
theSize(theSize_init), theCenter(theCenter_init), plastShearModulus(plas_modul)
{

}

MultiYieldSurface::~MultiYieldSurface()
{

}

void MultiYieldSurface::setData(const Vector & theCenter_init, 
                                double theSize_init, double plas_modul)
{
  theSize = theSize_init;
  theCenter = theCenter_init;
  plastShearModulus = plas_modul;
}

void MultiYieldSurface::setCenter(const Vector & newCenter)
{
  if (newCenter.Size() != 6) {
    opserr << "FATAL:MultiYieldSurface::setCenter(Vector &): vector size not equal 6" << endln;
    exit(-1);
  }

  theCenter = newCenter;
}


/**********************************************
ostream & operator<< (ostream & os, const MultiYieldSurface & a)
{
  os << "  theSize = " << a.theSize << endln 
     << "  theCenter = " << a.theCenter << endln
     << "  plastShearModulus = " << a.plastShearModulus << endln;
  
  return os;
}


istream & operator>> (istream & is, MultiYieldSurface & a)
{
  is >> a.theSize >> a.theCenter >> a.plastShearModulus;

  return is;
}
*********************************************/

double secondOrderEqn(double A, double B, double C, int i)
{
  if(A == 0){
    opserr << "FATAL:second_order_eqn: A=0." << endln;
    if(i==0) opserr << " when finding reference point on outer surface." <<endln;
    else opserr << " when moving active surface." <<endln;
    exit(-1);   
  }
  if(C == 0) return 0;
  if(B == 0){
    if(C/A > 0){
      opserr << "FATAL:second_order_eqn: Complex roots.\n";
      exit(-1);
    } 
    return sqrt(-C/A);
  }

	double determ, val1, val2, val;
  determ = B*B - 4.*A*C; 
  if(determ < 0){
    opserr << "FATAL:second_order_eqn: Complex roots.\n";
    if(i==0) opserr << " when finding reference point on outer surface." <<endln;
    else opserr << " when moving active surface." <<endln;
    opserr << "B2=" << B*B << " 4AC=" << 4.*A*C <<endln; 
    exit(-1);
  }
  
  if (B > 0) val1 = (-B - sqrt(determ)) / (2.*A);
  else val1 = (-B + sqrt(determ)) / (2.*A);
  val2 = C / (A * val1);

  if (val1 < 0 && val2 < 0){
		if (fabs(val1) < LOW_LIMIT) val1 = 0.;
		else if (fabs(val2) < LOW_LIMIT) val2 = 0.;
	}

  if (val1 < 0 && val2 < 0){
    opserr << "FATAL:second_order_eqn: Negative roots.\n";
    if(i==0) opserr << " when finding reference point on outer surface." <<endln;
    else opserr << " when moving active surface." <<endln;
		opserr << "A=" << A << " B=" << B << " C=" << C << " det=" << determ << 
			" x1=" << val1 << " x2=" << val2 << endln;  
    exit(-1);   
  }
  
  if (val1 < 0) return  val2;
  else if (val2 < 0) return  val1;
  else{
    val = val1;
    if(val > val2)  val = val2;
		return val;
  }
}
