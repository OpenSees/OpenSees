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

// $Revision: 1.0 $
// $Date: 2012-06-04 00:00:00 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/AxialSp.cpp,v $

// Written: Kazuki Tanimoto
// Created: June 2012
//
// Top and bottom series of axial springs for three-dimensonal multi-spring mechanical model.
//
// Description: This file contains the class definition for AxialSp.
//
//   strain  --> disp.
//   stress  --> force
//   modulus --> stiffness
//
// Description: This file contains the function to parse the TCL input
// uniaxialMaterial AxialSp matTag? sce? fty? fcy? < bte? bty? bcy? fcr? >


#include <Vector.h>
#include <string.h>

#include <AxialSp.h>
#include <Channel.h>
#include <math.h>
#include <float.h>

#include <elementAPI.h>

void* OPS_AxialSp()
{
    int numdata = OPS_GetNumRemainingInputArgs();
    if (numdata < 4) {
	opserr << "WARNING invalid number of arguments\n";
	return 0;
    }

    int tag;
    numdata = 1;
    if (OPS_GetIntInput(&numdata, &tag) < 0) {
	opserr << "WARNING invalid AxialSp tag\n";
	return 0;
    }

    double data[3];
    numdata = 3;
    if (OPS_GetDoubleInput(&numdata, data) < 0) {
	opserr << "WARNING invalid double inputs\n";
	return 0;
    }

    double opt[4] = {0,0,0,0};
    numdata = OPS_GetNumRemainingInputArgs();
    if (numdata > 4) numdata = 4;
    if (OPS_GetDoubleInput(&numdata, opt) < 0) {
	opserr << "WARNING invalid double inputs\n";
	return 0;
    }

    return new AxialSp(tag,data[0],data[1],data[2],opt[0],opt[1],opt[2],opt[3]);
}


AxialSp::AxialSp(int tag, double sce, double fty, double fcy,
		 double bte, double bty, double bcy, double fcr)
  :UniaxialMaterial(tag,MAT_TAG_AxialSp),sce(sce),fty(fty),fcy(fcy),
   bte(bte),bty(bty),bcy(bcy),fcr(fcr)
{

  if ( fty < 0.0 ) {
    opserr << "WARNING invalid fty\n";
	opserr << "fty>=0\n";
    opserr << "uniaxialMaterial AxialSp: " << tag << endln;
  }

  if ( fcy > 0.0 ) {
    opserr << "WARNING invalid fcy\n";
	opserr << "fcy<=0\n";
    opserr << "uniaxialMaterial AxialSp: " << tag << endln;
  }

  if ( !(bte >= 0.0 && bte <= 1.0) ) {
    opserr << "WARNING invalid bte\n";
	opserr << "0<=bte<=1\n";
    opserr << "uniaxialMaterial AxialSp: " << tag << endln;
  }

  if ( !(bty >= 0.0 && bty <= 1.0) ) {
    opserr << "WARNING invalid bty\n";
	opserr << "0<=bty<=1\n";
	opserr << "uniaxialMaterial AxialSp: " << tag << endln;
  }

  if ( !(bcy >= 0.0 && bcy <= 1.0) ) {
    opserr << "WARNING invalid bcy\n";
	opserr << "0<=bcy<=1\n";
    opserr << "uniaxialMaterial AxialSp: " << tag << endln;
  }

  if ( !(fcr <= 0.0 && fcr >= fcy) ) {
    opserr << "WARNING invalid fcr\n";
	opserr << "0<=fcr<=1\n";
    opserr << "uniaxialMaterial AxialSp: " << tag << endln;
  }

  this->revertToStart();
}

AxialSp::AxialSp()
  :UniaxialMaterial(0,MAT_TAG_AxialSp)
{
  trialDeformation   = 0.0;
  trialForce   = 0.0;
  trialStiffness  = 0.0;
  commitDeformation  = 0.0;
  commitForce  = 0.0;
  commitStiffness = 0.0;
}

AxialSp::~AxialSp()
{

}

int
AxialSp::setTrialStrain(double strain, double strainRate)
{
  //
  trialDeformation = strain;

  // Stg index
  switch(trialStg) {

  //1
  case 1:

    //
    if ( trialDeformation < ucy )
      trialStg = 6;

    //
    else if ( trialDeformation >= ucy && trialDeformation <= 0.0 )
      trialStg = 1;

    //
    else if ( trialDeformation > 0.0 && trialDeformation <= uty )
      trialStg = 2;

    //
    else if ( trialDeformation > uty )
      trialStg = 3;

    break;

  //2
  case 2:

    //
    if ( trialDeformation < ucy )
      trialStg = 6;

    //
    else if ( trialDeformation >= ucy && trialDeformation <= 0.0 )
      trialStg = 1;

    //
    else if ( trialDeformation > 0.0 && trialDeformation <= uty )
      trialStg = 2;

    //
    else if ( trialDeformation > uty )
      trialStg = 3;

    break;

  //3
  case 3:

    //
    if ( trialDeformation >= commitDeformation )
      trialStg = 3;

    //
    else {
      ur1 = commitDeformation;
      fr1 = commitForce;

      if ( trialDeformation > ucr )
	trialStg = 4;

      else if ( trialDeformation > ucy && trialDeformation <= ucr )
	trialStg = 6;

    }

    break;

  //4
  case 4:

    //
    if ( trialDeformation <= commitDeformation ) {
      if ( trialDeformation > ucr )
	trialStg = 4;

      else if ( trialDeformation > ucy && trialDeformation <= ucr )
	trialStg = 1;

      else if ( trialDeformation <= ucy )
	trialStg = 6;
    }

    //
    else {
      ur2 = commitDeformation;
      fr2 = commitForce;
      ur3 = (ste*ur2-sty*uty+fty-fr2)/(ste-sty);
      fr3 = sty*(ur3-uty) + fty;

      if ( trialDeformation > ur3 )
	trialStg = 3;
      else
	trialStg = 5;
    }

    break;

  //5
  case 5:

    if ( trialDeformation > ur3 )
      trialStg = 3;
    else if ( trialDeformation > ur2 && trialDeformation <= ur3 )
      trialStg = 5;
    else if ( trialDeformation > ucr && trialDeformation <= ur2 )
      trialStg = 4;
	else if ( trialDeformation > ucy && trialDeformation <= ucr )
	  trialStg = 1;
    else if ( trialDeformation <= ucy )
      trialStg = 6;

    break;

  //6
  case 6:

    //
    if ( trialDeformation <= commitDeformation )
      trialStg = 6;

    //
    else {
      ur4 = commitDeformation;
      fr4 = commitForce;

      uc0 = -fr4/sce + ur4;
      ur5 = (ste*uc0-sty*uty+fty)/(ste-sty);
      fr5 = ste*(ur5-uc0);

      if ( trialDeformation <= uc0 )
	trialStg = 7;
      else if ( trialDeformation > uc0 && trialDeformation <= ur5 )
	trialStg = 8;
      else if ( trialDeformation > ur5 && trialDeformation <= uty )
	trialStg = 9;
      else if ( trialDeformation > uty )
	trialStg = 3;
    }

    break;

  //7
  case 7:

    if ( trialDeformation <= ur4 )
      trialStg = 6;
    else if ( trialDeformation > ur4 && trialDeformation <= uc0 )
      trialStg = 7;
    else if ( trialDeformation > uc0 && trialDeformation <= ur5 )
      trialStg = 8;
    else if ( trialDeformation > ur5 && trialDeformation <= uty )
      trialStg = 9;
    else if ( trialDeformation > uty )
      trialStg =3;

    break;

  //8
  case 8:

    if ( trialDeformation <= ur4 )
      trialStg = 6;
    else if ( trialDeformation > ur4 && trialDeformation <= uc0 )
      trialStg = 7;
    else if ( trialDeformation > uc0 && trialDeformation <= ur5 )
      trialStg = 8;
    else if ( trialDeformation > ur5 && trialDeformation <= uty )
      trialStg = 9;
    else if ( trialDeformation > uty )
      trialStg = 3;

    break;

  //9
  case 9:

    //
    if ( trialDeformation >= commitDeformation ) {
      if ( trialDeformation <= uty )
	trialStg = 9;
      else
	trialStg = 3;
    }

    //
    else {
      ur5 = commitDeformation;
      fr5 = commitForce;
      uc0 = -fr5/ste + ur5;
      ur4 = (sce*uc0-scy*ucy+fcy)/(sce-scy);
      fr4 = scy*(ur4-ucy) + fcy;

      if ( trialDeformation <= ur4 )
	trialStg = 6;
      else if ( trialDeformation > ur4 && trialDeformation <= uc0 )
	trialStg = 7;
      else if ( trialDeformation > uc0 )
	trialStg = 8;
    }

    break;

  }

  // stiffness, force
  switch(trialStg) {

  // 1: compression, elastic
  case 1:

    trialStiffness = sce;
    trialForce = sce*trialDeformation;
    break;

  // 2: tension, elastic
  case 2:

    trialStiffness = ste;
    trialForce = ste*trialDeformation;

    break;

  // 3: after tensile yield point (skeleton curve)
  case 3:

    trialStiffness = sty;
    trialForce = sty*(trialDeformation-uty) + fty;

    break;

  // 4: after tensile yield point (unload)
  case 4:

    trialStiffness = (fcr-fr1)/(ucr-ur1);
    trialForce = trialStiffness*(trialDeformation-ur1) + fr1;

    break;

  // 5: after tensile yield point (unload, approach reversal point)
  case 5:

    trialStiffness = ste;
    trialForce = ste*(trialDeformation-ur2) + fr2;

    break;

  // 6: after compressive yield point (skeleton curve)
  case 6:

    trialStiffness = scy;
    trialForce = scy*(trialDeformation-ucy) + fcy;

    break;

  // 7: after compressive yield point (unload)
  case 7:

    trialStiffness = sce;
    trialForce = sce*(trialDeformation-ur4) + fr4;

    break;

  // 8: after compressive yield point (unload, turn into tensile stress)
  case 8:

    trialStiffness = ste;
    trialForce = ste*(trialDeformation-uc0);

    break;

  // 8: after compressive yield point (unload, turn into tensile stress, ur5<=strain)
  case 9:

    trialStiffness = sty;
    trialForce = sty*(trialDeformation-ur5) + fr5;

    break;

  }

  return 0;
}

double
AxialSp::getStress(void)
{
  return trialForce;
}

double
AxialSp::getTangent(void)
{
  return trialStiffness;
}

double
AxialSp::getInitialTangent(void)
{
  return sce;
}

double
AxialSp::getStrain(void)
{
  return trialDeformation;
}

int
AxialSp::commitState(void)
{
  commitDeformation  = trialDeformation;
  commitForce  = trialForce;
  commitStiffness = trialStiffness;

  commitStg = trialStg;

  return 0;
}

int
AxialSp::revertToLastCommit(void)
{
  trialDeformation  = commitDeformation;
  trialForce  = commitForce;
  trialStiffness = commitStiffness;

  trialStg = commitStg;

  return 0;
}

int
AxialSp::revertToStart(void)
{
  trialDeformation = 0.0;
  trialForce = 0.0;
  trialStiffness = sce;
  commitDeformation = 0.0;
  commitForce = 0.0;
  commitStiffness = sce;

  trialStg = 1;
  commitStg = 1;

  ste  = bte * sce;
  sty  = bty * sce;
  scy  = bcy * sce;
  uty  = fty / ste;
  ucy  = fcy / sce;
  ucr  = fcr / sce;
  uc0  = 0.0;
  ur1  = 0.0;
  fr1  = 0.0;
  ur2  = 0.0;
  fr2  = 0.0;
  ur3  = 0.0;
  fr3  = 0.0;
  ur4  = 0.0;
  fr4  = 0.0;
  ur5  = 0.0;
  fr5  = 0.0;

  return 0;
}

UniaxialMaterial *
AxialSp::getCopy(void)
{
  AxialSp *theCopy = new AxialSp(this->getTag(), sce, fty, fcy,
				 bte, bty, bcy, fcr);

  // Copy committed history variables
  theCopy->commitStg = commitStg;

  // Copy trial history variables
  theCopy->trialStg = trialStg;

  // Copy trial state variables
  theCopy->trialDeformation = trialDeformation;
  theCopy->trialForce = trialForce;
  theCopy->trialStiffness = trialStiffness;

  // Copy commit state variables
  theCopy->commitDeformation = commitDeformation;
  theCopy->commitForce = commitForce;
  theCopy->commitStiffness = commitStiffness;

  return theCopy;
}

int
AxialSp::sendSelf(int cTag, Channel &theChannel)
{
  int res = 0;

  static Vector data(16+11);

  data(0)  = this->getTag();
  data(1)  = sce;
  data(2)  = fty;
  data(3)  = fcy;
  data(4)  = bte;
  data(5)  = bty;
  data(6)  = bcy;
  data(7)  = fcr;
  data(8)  = commitDeformation;
  data(9)  = commitForce;
  data(10) = commitStiffness;
  data(11) = commitStg;
  data(12) = trialDeformation;
  data(13) = trialForce;
  data(14) = trialStiffness;
  data(15) = trialStg;

  data(16) = uc0;
  data(17) = ur1;
  data(18) = fr1;
  data(19) = ur2;
  data(20) = fr2;
  data(21) = ur3;
  data(22) = fr3;
  data(23) = ur4;
  data(24) = fr4;
  data(25) = ur5;
  data(26) = fr5;
  
  res = theChannel.sendVector(this->getDbTag(), cTag, data);
  if ( res < 0 )
    opserr << "AxialSp::sendSelf() - failed to send data\n";

  return res;
}

int
AxialSp::recvSelf(int cTag, Channel &theChannel,
		  FEM_ObjectBroker &theBroker)
{
  int res = 0;

  static Vector data(16+11);
  res = theChannel.recvVector(this->getDbTag(), cTag, data);

  if ( res < 0 ) {
    opserr << "AxialSp::recvSelf() - failed to receive data\n";
    this->setTag(0);
  }

  else {
    this->setTag((int)data(0));
    sce = data(1);
    fty = data(2);
    fcy = data(3);
    bte = data(4);
    bty = data(5);
    bcy = data(6);
    fcr = data(7);
    commitDeformation = data(8);
    commitForce = data(9);
    commitStiffness = data(10);
    commitStg = (int)data(11);
    trialDeformation = data(12);
    trialForce = data(13);
    trialStiffness = data(14);
    trialStg = (int)data(15);

    ste  = bte * sce;
    sty  = bty * sce;
    scy  = bcy * sce;
    uty  = fty / ste;
    ucy  = fcy / sce;
    ucr  = fcr / sce;

    uc0 = data(16);
    ur1 = data(17);
    fr1 = data(18);
    ur2 = data(19);
    fr2 = data(20);
    ur3 = data(21);
    fr3 = data(22);
    ur4 = data(23);
    fr4 = data(24);
    ur5 = data(25);
    fr5 = data(26);
  }

  return res;
}

void
AxialSp::Print(OPS_Stream &s, int flag)
{
  s << "AxialSp : " << this->getTag() << endln;
  s << "    sce : " << endln;
  s << "    fty : " << endln;
  s << "    fcy : " << endln;
  s << "    bte : " << endln;
  s << "    bty : " << endln;
  s << "    bcy : " << endln;
  s << "    fcr : " << endln;

  return;

}


