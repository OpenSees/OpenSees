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

                                                                        

// $Revision: 1.2 $                                                              

// $Date: 2009-01-08 22:00:17 $                                                                  

// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/backbone/ManderBackbone.cpp,v $                                                                



// Written: MHS

// Created: Mar 2001

//

// Description: This file contains the implementation of 

// ManderBackbone, which the concrete backbone function given

// by Mander, Priestly, and Park (1988)



#include <ManderBackbone.h>

#include <Vector.h>

#include <Channel.h>



#include <math.h>

#include <elementAPI.h>

void * OPS_ADD_RUNTIME_VPV(OPS_ManderBackbone)
{
  HystereticBackbone *theBackbone = 0;

  if (OPS_GetNumRemainingInputArgs() < 4) {
    opserr << "Invalid number of args, want: hystereticBackbone Mander tag? fc? epsc? E?" << endln;
    return 0;
  }

  int iData[1];
  double dData[3];
  
  int numData = 1;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid tag for hystereticBackbone Mander" << endln;
    return 0;
  }

  numData = 3;
  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "WARNING invalid data for hystereticBackbone Mander" << endln;
    return 0;
  }

  theBackbone = new ManderBackbone(iData[0], dData[0], dData[1], dData[2]);
  if (theBackbone == 0) {
    opserr << "WARNING could not create ManderBackbone\n";
    return 0;
  }

  return theBackbone;
}

ManderBackbone::ManderBackbone(int tag, double f, double e, double E):
  HystereticBackbone(tag,BACKBONE_TAG_Mander),
  fpc(f), epsc(e), Ec(E)
{
  fpc = fabs(fpc);
  epsc = fabs(epsc);
  Ec = fabs(Ec);

  /*

  if (Ec <= fpc/epsc) {

    opserr << "ManderBackbone::ManderBackbone -- Ec <= Esec, setting Ec = 2*fpc/epsc" << endln;

    Ec = 2*fpc/epsc;

  }

  */

}



ManderBackbone::ManderBackbone():

  HystereticBackbone(0,BACKBONE_TAG_Mander),

  fpc(0.0), epsc(0.0), Ec(0.0)

{



}



ManderBackbone::~ManderBackbone()

{



}



double

ManderBackbone::getTangent (double strain)

{

  if (strain > 0.0)

    return 0.0;

  

  strain *= -1;



  double oneOverepsc = 1.0/epsc;

  

  double x = strain*oneOverepsc;

  double Esec = fpc*oneOverepsc;

  

  double r = Ec/(Ec-Esec);

  

  double xr = pow(x,r);

  double denom = r-1.0+xr;

  

  return Esec*r*(r-1.0)*(1.0-xr)/(denom*denom);

}



double

ManderBackbone::getStress (double strain)

{

  if (strain > 0.0)

    return 0.0;

  

  strain *= -1;



  double oneOverepsc = 1.0/epsc;

  

  double x = strain*oneOverepsc;

  double Esec = fpc*oneOverepsc;

  

  double r = Ec/(Ec-Esec);

  

  return -fpc*(x*r)/(r-1.0+pow(x,r));

}



double

ManderBackbone::getEnergy (double strain)

{

  return 0.0;

}



double

ManderBackbone::getYieldStrain(void)

{

  return epsc;

}



HystereticBackbone*

ManderBackbone::getCopy(void)

{

  ManderBackbone *theCopy =

    new ManderBackbone (this->getTag(), fpc, epsc, Ec);

  

  return theCopy;

}



void

ManderBackbone::Print(OPS_Stream &s, int flag)

{

  s << "ManderBackbone, tag: " << this->getTag() << endln;

  s << "\tfpc: " << fpc << endln;

  s << "\tepsc: " << epsc << endln;

  s << "\tEc: " << Ec << endln;

}



int

ManderBackbone::setVariable (char *argv)

{

  return -1;

}



int

ManderBackbone::getVariable (int varID, double &theValue)

{

  return -1;

}



int

ManderBackbone::sendSelf(int commitTag, Channel &theChannel)

{

  int res = 0;

  

  static Vector data(4);

  

  data(0) = this->getTag();

  data(1) = fpc;

  data(2) = epsc;

  data(3) = Ec;

  

  res += theChannel.sendVector(this->getDbTag(), commitTag, data);

  if (res < 0) {

    opserr << "ManderBackbone::sendSelf -- could not send Vector" << endln;

    return res;

  }

  

  return res;

}



int

ManderBackbone::recvSelf(int commitTag, Channel &theChannel, 

			 FEM_ObjectBroker &theBroker)

{

  int res = 0;

  

  static Vector data(4);

  

  res += theChannel.recvVector(this->getDbTag(), commitTag, data);

  if (res < 0) {

    opserr << "ManderBackbone::recvSelf -- could not receive Vector" << endln;

    return res;

  }

  

  this->setTag(int(data(0)));

  fpc = data(1);

  epsc = data(2);

  Ec = data(3);

  

  return res;

}

