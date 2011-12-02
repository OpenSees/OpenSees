// $Revision: 1.1 $
// $Date: 2000-12-19 03:35:02 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/nD/soil/T2Vector.cpp,v $
                                                                        
// Written: ZHY
// Created: August 2000

//
// T2Vector.cpp
// ----------
//
#include <iostream.h>
#include <fstream.h>
#include <iomanip.h>
#include <math.h>
#include <stdlib.h>
#include "T2Vector.h"

#define UP_LIMIT    10e+30
#define LOW_LIMIT   10e-15


double operator && (const Vector & a, const Vector & b)
{
	if (a.Size() !=6 || b.Size() !=6) {
	cerr << "FATAL:operator && (Vector &, Vector &): vector size not equal 6" << endl;
	g3ErrorHandler->fatal("");
  }

  double result = 0.;  

  for (int i=0; i<3; i++)
    result += a[i]*b[i] + 2*a[i+3]*b[i+3];
  return result;
}


// T2Vector class methods
T2Vector::T2Vector() : theT2Vector(6), theDeviator(6), theVolume(0.)
{ 
}


T2Vector::T2Vector(const Vector & Vector_init) : theT2Vector(6), theDeviator(6)
{
	if ( Vector_init.Size() != 6) {
	  cerr << "FATAL:T2Vector::T2Vector(Vector &): vector size not equal to 6" << endl;
	  g3ErrorHandler->fatal("");
  }
  theT2Vector = Vector_init;

  theVolume = (theT2Vector[0]+theT2Vector[1]+theT2Vector[2])/3.0;
  for(int i=0; i<3; i++){
    theDeviator[i] = theT2Vector[i] - theVolume;
    theDeviator[i+3] = theT2Vector[i+3];
  }
}


T2Vector::T2Vector(const Vector & deviat_init, double volume_init)
 : theT2Vector(6), theDeviator(6), theVolume(volume_init)
{
	if (deviat_init.Size() != 6) {
	  cerr << "FATAL:T2Vector::T2Vector(Vector &, double): vector size not equal 6" << endl;
	  g3ErrorHandler->fatal("");
  }

	//make sure the deviator has truely volume=0 
	double devolum = (deviat_init[0]+deviat_init[1]+deviat_init[2])/3.;

  for(int i=0; i<3; i++){
    theDeviator[i] = deviat_init[i] - devolum;
    theDeviator[i+3] = deviat_init[i+3];
    theT2Vector[i] = theDeviator[i] + theVolume;
    theT2Vector[i+3] = theDeviator[i+3]; 
  }
}


T2Vector::~T2Vector()
{
}

double T2Vector::t2VectorLength() const
{
  return sqrt(theT2Vector && theT2Vector);
}


double T2Vector::deviatorLength() const
{
  return sqrt(theDeviator && theDeviator);
}


double T2Vector::octahedralShear(int isEngrgStain) const
{
  if (isEngrgStain) 
    return 2.* sqrt(1. / 3.) * deviatorLength();
  else
    return sqrt(1. / 3.) * deviatorLength();
}


double T2Vector::deviatorRatio(double residualPress) const
{
  if ((fabs(theVolume)+fabs(residualPress)) <= LOW_LIMIT) {
	cerr << "FATAL:T2Vector::deviatorRatio(): volume <=" << LOW_LIMIT << endl;
	g3ErrorHandler->fatal("");
  }
  return sqrt(3./2.* (theDeviator && theDeviator)) / (fabs(theVolume)+fabs(residualPress));
}


const Vector T2Vector::unitT2Vector() const
{
	if (t2VectorLength() <= LOW_LIMIT) {
	cerr << "WARNING:T2Vector::unitT2Vector(): vector length <=" << LOW_LIMIT << endl;
	return theT2Vector;
  }

  return theT2Vector/t2VectorLength();
}


const Vector T2Vector::unitDeviator() const
{
	if (deviatorLength() <= LOW_LIMIT) {
	cerr << "WARNING:T2Vector::unitDeviator(): vector length <=" << LOW_LIMIT << endl;
	return theDeviator;
  }

  return theDeviator/deviatorLength();
}


double T2Vector::angleBetweenT2Vector(const T2Vector & a) const
{
	if (t2VectorLength() <= LOW_LIMIT || a.t2VectorLength() <= LOW_LIMIT) {
	cerr << "FATAL:T2Vector::angleBetweenT2Vector(T2Vector &): vector length <=" << LOW_LIMIT << endl;
	g3ErrorHandler->fatal("");
  }

  double angle = (theT2Vector && a.theT2Vector) / (t2VectorLength() * a.t2VectorLength());
  if(angle > 1.) angle = 1.;
  if(angle < -1.) angle = -1.;

  return acos(angle);
}


double T2Vector::angleBetweenDeviator(const T2Vector & a) const
{
	if (deviatorLength() <= LOW_LIMIT || a.deviatorLength() <= LOW_LIMIT) {
	cerr << "FATAL:T2Vector::angleBetweenDeviator(T2Vector &): vector length <=" << LOW_LIMIT << endl;
	g3ErrorHandler->fatal("");
  }

  double angle = (theDeviator && a.theDeviator) / (deviatorLength() * a.deviatorLength());
  if(angle > 1.) angle = 1.;
  if(angle < -1.) angle = -1.;

  return acos(angle);
}


int T2Vector::operator == (const T2Vector & a) const
{
  for(int i=0; i<6; i++)
    if(theT2Vector[i] != a.theT2Vector[i]) return 0;

  return 1;
}


int T2Vector::sendSelf(int commitTag, Channel &theChannel)
{
	// Need to implement
	return 0;
}


int T2Vector::recvSelf(int commitTag, Channel &theChannel, 
		 FEM_ObjectBroker &theBroker)    
{
	// Need to implement
	return 0;
}


ostream & operator<< (ostream & os, const T2Vector & a)
{
  os.precision(16);
  os.setf(ios::showpoint);

  os << "theT2Vector = " << a.t2Vector() << endl;
  os << "theDeviator = " << a.deviator() << endl;
  os << "theVolume = " << a.volume() << endl;

  return os;
}


istream & operator>> (istream & is, T2Vector & a)
{
  Vector temp;

  is >> temp;
  a = T2Vector(temp);

  return is;
}

