// $Revision: 1.1 $
// $Date: 2000-12-19 03:35:02 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/nD/soil/MultiYieldSurface.cpp,v $
                                                                        
// Written: ZHY
// Created: August 2000
//
// MultiYieldSurface.cpp
// ---------------------
//
#include <iostream.h>
#include <fstream.h>
#include <iomanip.h>
#include <math.h>
#include <stdlib.h>

#include "MultiYieldSurface.h"

const Vector zeroVector = Vector(6);


// YieldSurface class methods
MultiYieldSurface::MultiYieldSurface()
{
  theSize = plastShearModulus = 0.;
  theCenter = zeroVector;
}

MultiYieldSurface::MultiYieldSurface(const Vector & theCenter_init, 
                                     double theSize_init, double plas_modul) 
{
  theCenter = theCenter_init;
  theSize = theSize_init;
  plastShearModulus = plas_modul;
}

MultiYieldSurface::~MultiYieldSurface()
{
}

void MultiYieldSurface::setCenter(const Vector & newCenter)
{
  if (newCenter.Size() != 6) {
	cerr << "FATAL:MultiYieldSurface::setCenter(Vector &): vector size not equal 6" << endl;
	g3ErrorHandler->fatal("");
  }

  theCenter = newCenter;
}


ostream & operator<< (ostream & os, const MultiYieldSurface & a)
{
  os << "  theSize = " << a.theSize << endl 
     << "  theCenter = " << a.theCenter << endl
     << "  plastShearModulus = " << a.plastShearModulus << endl;
  
  return os;
}


istream & operator>> (istream & is, MultiYieldSurface & a)
{
  is >> a.theSize >> a.theCenter >> a.plastShearModulus;

  return is;
}


double secondOrderEqn(double A, double B, double C, int i)
{
  if(A == 0){
    cerr << "FATAL:second_order_eqn: A=0." << endl;
    if(i==0) cerr << " when finding reference point on outer surface." <<endl;
    else cerr << " when moving active surface." <<endl;
    g3ErrorHandler->fatal("");   
  }
  if(C == 0) return 0;
  if(B == 0){
    if(C/A > 0){
      cerr << "FATAL:second_order_eqn: Complex roots.\n";
      g3ErrorHandler->fatal("");
    } 
    return sqrt(-C/A);
  }

	double determ, val1, val2, val;
  determ = B*B - 4.*A*C; 
  if(determ < 0){
    cerr << "FATAL:second_order_eqn: Complex roots.\n";
    if(i==0) cerr << " when finding reference point on outer surface." <<endl;
    else cerr << " when moving active surface." <<endl;
    cerr << "B2=" << B*B << " 4AC=" << 4.*A*C <<endl; 
    g3ErrorHandler->fatal("");
  }
  
  if (B > 0) val1 = (-B - sqrt(determ)) / (2.*A);
  else val1 = (-B + sqrt(determ)) / (2.*A);
  val2 = C / (A * val1);

  if (val1 < 0 && val2 < 0){
		if (fabs(val1) < LOW_LIMIT) val1 = 0.;
		else if (fabs(val2) < LOW_LIMIT) val2 = 0.;
	}

  if (val1 < 0 && val2 < 0){
    cerr << "FATAL:second_order_eqn: Negative roots.\n";
    if(i==0) cerr << " when finding reference point on outer surface." <<endl;
    else cerr << " when moving active surface." <<endl;
		cerr << "A=" << A << " B=" << B << " C=" << C << " det=" << determ << 
			" x1=" << val1 << " x2=" << val2 << endl;  
    g3ErrorHandler->fatal("");   
  }
  
  if (val1 < 0) return  val2;
  else if (val2 < 0) return  val1;
  else{
    val = val1;
    if(val > val2)  val = val2;
		return val;
  }
}
