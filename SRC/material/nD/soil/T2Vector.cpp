// $Revision: 1.11 $
// $Date: 2009-07-23 23:57:27 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/nD/soil/T2Vector.cpp,v $
                                                                        
// Written: ZHY
// Created: August 2000

//
// T2Vector.cpp
// ----------
//

#include <math.h>
#include <float.h>
#include <stdlib.h>
#include <T2Vector.h>
#include <Matrix.h>



Vector T2Vector::engrgStrain(6);

double operator && (const Vector & a, const Vector & b)
{
  if (a.Size() !=6 || b.Size() !=6) {
    opserr << "FATAL:operator && (Vector &, Vector &): vector size not equal 6" << endln;
    exit(-1);
  }

  double result = 0.;  

  for (int i=0; i<3; i++)
    result += a[i]*b[i] + 2*a[i+3]*b[i+3];
  return result;
}

// ---------------- add by guquan ------------------------------
// ---------------- c=a:b, c(k,l)=a(i,j)*b(i,j,k,l)-------------
void doubledotProduct (Vector & c, const Vector & a, const Matrix & b)
{
  if (c.Size() !=6 || a.Size() !=6 || b.noCols() !=6|| b.noRows() !=6) {
    opserr << "FATAL:operator && (Vector &, Matrix &): vector or Matrix size not equal 6" << endln;
    exit(-1);
  }
	
  c.Zero();
  for(int j=0;j<6;j++){
	  for (int i=0; i<3; i++){
		c[j] += a[i]*b(i,j) + 2*a[i+3]*b(i+3,j);
	  }
  }
  return;
}


// ---------------- add by guquan ------------------------------
// ---------------- c=a:b, c(i,j,k,l)=a(i,j,m,n)*b(m,n,k,l)-------------
void doubledotMatrixProduct (Matrix & c, const Matrix & a, const Matrix & b)
{
  if (c.noCols() !=6 ||c.noRows() !=6 || a.noCols() !=6 ||a.noRows() !=6 || b.noCols() !=6|| b.noRows() !=6) {
    opserr << "FATAL: doubledotproduct(Matrix &, Matrix &): Matrix size not equal 6" << endln;
    exit(-1);
  }
	
  c.Zero();
  for(int i=0;i<6;i++){
	for(int j=0;j<6;j++){
	  for (int l=0; l<3; l++){
		c(i,j) += a(i,l)*b(l,j) + 2*a(i,l+3)*b(l+3,j);
	  }
	}
  }
  return;
}
 
// ---------------- add by guquan ------------------------------
// ---------------- c=a*b, c(i,j,k,l)=a(i,j)*b(k,l)-------------
void tensorProduct(Matrix & c, const Vector & a, const Vector & b)
{
  if (b.Size() !=6 || a.Size() !=6 || c.noCols() !=6|| c.noRows() !=6) {
    opserr << "FATAL:operator && (Vector &, Matrix &): vector or Matrix size not equal 6" << endln;
    exit(-1);
  }
	
  c.Zero();
  for(int j=0;j<6;j++){
	  for (int i=0; i<6; i++){
		c(i,j) = a[i]*b[j];
	  }
  }
  return;
}

double delta(int i,int j){
	if (i==j) return 1.0;
	else return 0.0;
}


// T2Vector class methods
T2Vector::T2Vector() 
:theT2Vector(6), theDeviator(6), theVolume(0.0)
{
	
}


T2Vector::T2Vector(const Vector &init, int isEngrgStrain)
:theT2Vector(6), theDeviator(6), theVolume(0)
{
  if (init.Size() != 6) {
    opserr << "FATAL:T2Vector::T2Vector(Vector &): vector size not equal to 6" << endln;
    exit(-1);
  }
  theT2Vector = init;

  theVolume = (theT2Vector[0]+theT2Vector[1]+theT2Vector[2])/3.0;
  for(int i=0; i<3; i++){
    theDeviator[i] = theT2Vector[i] - theVolume;
    theDeviator[i+3] = theT2Vector[i+3];
    if (isEngrgStrain==1) {
      theDeviator[i+3] /= 2.;
      theT2Vector[i+3] /= 2.;
    }
  }
}



T2Vector::T2Vector(const Vector & deviat_init, double volume_init)
 : theT2Vector(6), theDeviator(6), theVolume(volume_init)
{
  if (deviat_init.Size() != 6) {
    opserr << "FATAL:T2Vector::T2Vector(Vector &, double): vector size not equal 6" << endln;
    exit(-1);
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


void
T2Vector::setData(const Vector &init, int isEngrgStrain)
{
  if ( init.Size() != 6) {
    opserr << "FATAL:T2Vector::T2Vector(Vector &): vector size not equal to 6" << endln;
    exit(-1);
  }
  theT2Vector = init;

  theVolume = (theT2Vector[0]+theT2Vector[1]+theT2Vector[2])/3.0;
  for(int i=0; i<3; i++){
    theDeviator[i] = theT2Vector[i] - theVolume;
    theDeviator[i+3] = theT2Vector[i+3];
    if (isEngrgStrain==1) {
      theDeviator[i+3] /= 2.;
      theT2Vector[i+3] /= 2.;
    }
  }
}

void
T2Vector::setData(const Vector & deviat, double volume)
{
  theVolume = volume;
  
  if (deviat.Size() != 6) {
    opserr << "FATAL:T2Vector::T2Vector(Vector &, double): vector size not equal 6" << endln;
    exit(-1);
  }

  //make sure the deviator has truely volume=0 
  double devolum = (deviat[0]+deviat[1]+deviat[2])/3.;

  for(int i=0; i<3; i++){
    theDeviator[i] = deviat[i] - devolum;
    theDeviator[i+3] = deviat[i+3];
    theT2Vector[i] = theDeviator[i] + theVolume;
    theT2Vector[i+3] = theDeviator[i+3]; 
  }
}

const Vector & 
T2Vector::t2Vector(int isEngrgStrain) const
{
  if (isEngrgStrain==0) return theT2Vector;

  engrgStrain = theT2Vector;
  for(int i=0; i<3; i++){
    engrgStrain[i+3] *= 2.;
  }
  return engrgStrain;
}


const Vector & T2Vector::deviator(int isEngrgStrain) const
{
  if (isEngrgStrain==0) return theDeviator;

  engrgStrain = theDeviator;
  for(int i=0; i<3; i++){
    engrgStrain[i+3] *= 2.;
  }
  return engrgStrain;
}


double 
T2Vector::t2VectorLength() const
{
  return sqrt(theT2Vector && theT2Vector);
}


double 
T2Vector::deviatorLength() const
{
  return sqrt(theDeviator && theDeviator);
}


double 
T2Vector::octahedralShear(int isEngrgStain) const
{
  if (isEngrgStain) 
    return 2.* sqrt(1. / 3.) * deviatorLength();
  else
    return sqrt(1. / 3.) * deviatorLength();
}


double 
T2Vector::deviatorRatio(double residualPress) const
{
  if ((fabs(theVolume)+fabs(residualPress)) <= LOW_LIMIT) {
	opserr << "FATAL:T2Vector::deviatorRatio(): volume <=" << LOW_LIMIT << endln;
	exit(-1);
  }
  return sqrt(3./2.* (theDeviator && theDeviator)) / (fabs(theVolume)+fabs(residualPress));
}


const Vector &
T2Vector::unitT2Vector() const
{
  engrgStrain = theT2Vector;	
  double length = this->t2VectorLength();
  if (length <= LOW_LIMIT) {
    opserr << "WARNING:T2Vector::unitT2Vector(): vector length <=" << LOW_LIMIT << endln;
    engrgStrain /= LOW_LIMIT;
  } else
    engrgStrain /= length;

  return engrgStrain;
}


const Vector &
T2Vector::unitDeviator() const
{

  engrgStrain = theDeviator;;	
  double length = this->deviatorLength();
  if (length <= LOW_LIMIT) {
    opserr << "WARNING:T2Vector::unitT2Vector(): vector length <=" << LOW_LIMIT << endln;
    engrgStrain /= LOW_LIMIT;
  } else
    engrgStrain /= length;

  return engrgStrain;
}


double 
T2Vector::angleBetweenT2Vector(const T2Vector & a) const
{
  if (t2VectorLength() <= LOW_LIMIT || a.t2VectorLength() <= LOW_LIMIT) {
    opserr << "FATAL:T2Vector::angleBetweenT2Vector(T2Vector &): vector length <=" << LOW_LIMIT << endln;
    exit(-1);
  }

  double angle = (theT2Vector && a.theT2Vector) / (t2VectorLength() * a.t2VectorLength());
  if(angle > 1.) angle = 1.;
  if(angle < -1.) angle = -1.;

  return acos(angle);
}


double 
T2Vector::angleBetweenDeviator(const T2Vector & a) const
{
  if (deviatorLength() <= LOW_LIMIT || a.deviatorLength() <= LOW_LIMIT) {
    opserr << "FATAL:T2Vector::angleBetweenDeviator(T2Vector &): vector length <=" << LOW_LIMIT << endln;
    exit(-1);
  }

  double angle = (theDeviator && a.theDeviator) / (deviatorLength() * a.deviatorLength());
  if(angle > 1.) angle = 1.;
  if(angle < -1.) angle = -1.;

  return acos(angle);
}


int 
T2Vector::operator == (const T2Vector & a) const
{
  for(int i=0; i<6; i++)
    if(theT2Vector[i] != a.theT2Vector[i]) return 0;

  return 1;
}

int 
T2Vector::isZero(void) const
{
  for(int i=0; i<6; i++)
    if(theT2Vector[i] != 0.0) return 0;

  return 1;
}
int T2Vector::Zero(void)
{
  theT2Vector.Zero();
  theDeviator.Zero();
  theVolume=0.0;
return 1;
}


/*********************
ostream & operator<< (ostream & os, const T2Vector & a)
{
  os.precision(16);
  os.setf(ios::showpoint);

  os << "theT2Vector = " << a.t2Vector() << endln;
  os << "theDeviator = " << a.deviator() << endln;
  os << "theVolume = " << a.volume() << endln;

  return os;
}


istream & operator>> (istream & is, T2Vector & a)
{
  Vector temp;

  is >> temp;
  a = T2Vector(temp);

  return is;
}
*/








