//<<<<<<< T2Vector.h
//<<<<<<< T2Vector.h
// $Revision: 1.7 $
// $Date: 2002-05-16 00:07:47 $
//=======
// $Revision: 1.7 $
// $Date: 2002-05-16 00:07:47 $
//>>>>>>> 1.4
//=======
// $Revision: 1.7 $
// $Date: 2002-05-16 00:07:47 $
//>>>>>>> 1.6
// $Source: /usr/local/cvs/OpenSees/SRC/material/nD/soil/T2Vector.h,v $
                                                                        
// Written: ZHY
// Created: August 2000

//
// T2Vector.h
// --------
//

#ifndef _T2Vector_H_
#define _T2Vector_H_

#include <Vector.h>
#include <Channel.h>
#include <float.h>

#define UP_LIMIT    1.0e+30
#define LOW_LIMIT   20.*DBL_EPSILON

// global function: scalar product of two second order tensor vectors
double operator && (const Vector &, const Vector &);

// define second order tensor vector class
class T2Vector 
{

public:
  // constructors
  T2Vector();
  T2Vector(const Vector & T2Vector_init, int isEngrgStrain=0);
  T2Vector(const Vector & deviat_init, double volume_init);
  
  ~T2Vector();

  void setData(const Vector &init, int isEngrgStrain =0);
  void setData(const Vector &deviat, double volume);

  const Vector & t2Vector(int isEngrgStrain=0) const; 
  const Vector & deviator(int isEngrgStrain=0) const;
  double volume() const {return theVolume; }
  const Vector &unitT2Vector() const;
  const Vector &unitDeviator() const;
  double t2VectorLength() const;
  double deviatorLength() const;
  double octahedralShear(int isEngrgStrain=0) const;

  // = -sqrt(3/2*(S:S))/(p+residualPress)
  double deviatorRatio(double residualPress=0.) const; 

  //next function return the angle between two T2Vectors in radians (-PI to PI)
  double angleBetweenT2Vector(const T2Vector &) const; 

  //next function return the angle between deviatoric components of
  //two vectors in radians (-PI to PI)
  double angleBetweenDeviator(const T2Vector &) const; 

  int operator == (const T2Vector & a) const;
  int isZero(void) const;

protected:

private:
  Vector theT2Vector;
  Vector theDeviator;
  double theVolume;
  static Vector engrgStrain;
};


#endif
