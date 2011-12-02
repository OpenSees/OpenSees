// $Revision: 1.1 $
// $Date: 2000-12-19 03:35:02 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/nD/soil/MultiYieldSurface.h,v $
                                                                        
// Written: ZHY
// Created: August 2000
//
// MultiYieldSurface.h
// -------------------
//

#ifndef _MultiYieldSurface_H_
#define _MultiYieldSurface_H_

#include "T2Vector.h"

#define UP_LIMIT    1.0e+30
#define LOW_LIMIT   1.0e-15
#define LOCK_VALUE  1.0e+30

// global function to find the roots of a second order equation
double secondOrderEqn(double A, double B, double C, int i);

// define yield surface in stress space
class MultiYieldSurface
{
 
public:
  //constructors
  MultiYieldSurface();
  MultiYieldSurface(const Vector & center_init, double size_init, 
                    double plas_modul); 
  ~MultiYieldSurface();
	const Vector & center() const {return theCenter; }
	double size() const {return theSize; }
	double modulus() const {return plastShearModulus; }
  void  setCenter(const Vector & newCenter);
  friend ostream & operator<< (ostream & os, const MultiYieldSurface & );  
	friend istream & operator>> (istream & is, MultiYieldSurface & );

protected:

private:
  double theSize;
  Vector theCenter;  
  double plastShearModulus;

};

#endif
