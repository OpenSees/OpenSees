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
// $Date: 2013-06-03 00:00:00 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/AxialSpHD.cpp,v $

// Written: Kazuki Tanimoto
// Created: June 2012
//
// Top and bottom series of axial springs for three-dimensonal multi-spring mechanical model.
// This model is considering hardening at tension side.
//
// Description: This file contains the class definition for AxialSpHD.
//              This file contains the function to parse the TCL input
//
//   strain  --> disp.
//   stress  --> force
//   modulus --> stiffness
//
// uniaxialMaterial AxialSpHD matTag? sce? fty? fcy? < bte? bty? bth? bcy? fcr? ath? >


#include <string.h>
#include <AxialSpHD.h>
#include <Vector.h>
#include <Channel.h>
#include <math.h>
#include <float.h>

#include <elementAPI.h>

void* OPS_AxialSpHD()
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

    double opt[6] = {1,1,1,1,0,1};
    numdata = OPS_GetNumRemainingInputArgs();
    if (numdata > 6) numdata = 6;
    if (OPS_GetDoubleInput(&numdata, opt) < 0) {
	opserr << "WARNING invalid double inputs\n";
	return 0;
    }

    return new AxialSpHD(tag,data[0],data[1],data[2],opt[0],opt[1],opt[2],opt[3],opt[4],opt[5]);
}



AxialSpHD::AxialSpHD(int tag, double sce, double fty, double fcy, double bte,
		 double bty, double bth, double bcy, double fcr, double ath)
  :UniaxialMaterial(tag,MAT_TAG_AxialSpHD),sce(sce),fty(fty),fcy(fcy),bte(bte),
   bty(bty),bth(bth),bcy(bcy),fcr(fcr),ath(ath)
{

  if ( fty < 0.0 ) {
    opserr << "WARNING invalid fty\n";
	opserr << "fty>=0\n";
    opserr << "uniaxialMaterial AxialSpHD: " << tag << endln;
  }

  if ( fcy > 0.0 ) {
    opserr << "WARNING invalid fcy\n";
	opserr << "fcy<=0\n";
    opserr << "uniaxialMaterial AxialSpHD: " << tag << endln;
  }

  if ( !(bte >= 0.0 && bte <= 1.0) ) {
    opserr << "WARNING invalid bte\n";
	opserr << "0<=bte<=1\n";
    opserr << "uniaxialMaterial AxialSpHD: " << tag << endln;
  }

  if ( !(bty >= 0.0 && bty <= 1.0) ) {
    opserr << "WARNING invalid bty\n";
	opserr << "0<=bty<=1\n";
    opserr << "uniaxialMaterial AxialSpHD: " << tag << endln;
  }

  if ( !(bth >= 0.0 && bth <= 1.0 && bth > bty && bth < bte) ) {
    opserr << "WARNING invalid bth\n";
	opserr << "0<=bth<=1 and bty<bth<bte\n";
    opserr << "uniaxialMaterial AxialSpHD: " << tag << endln;
  }

  if ( !(bcy >= 0.0 && bcy <= 1.0) ) {
    opserr << "WARNING invalid bcy\n";
	opserr << "0<=bcy<=1\n";
    opserr << "uniaxialMaterial AxialSpHD: " << tag << endln;
  }

  if ( !(fcr <= 0.0 && fcr >= fcy) ) {
    opserr << "WARNING invalid fcr\n";
	opserr << "0<=fcr<=fcy\n";
    opserr << "uniaxialMaterial AxialSpHD: " << tag << endln;
  }

  if ( ath < 1.0 ) {
    opserr << "WARNING invalid ath\n";
	opserr << "ath>=1\n";
    opserr << "uniaxialMaterial AxialSpHD: " << tag << endln;
  }

  this->revertToStart();
}

AxialSpHD::AxialSpHD()
  :UniaxialMaterial(0,MAT_TAG_AxialSpHD)
{
  trialDeformation = 0.0;
  trialForce   = 0.0;
  trialStiffness  = 0.0;
  commitDeformation  = 0.0;
  commitForce  = 0.0;
  commitStiffness = 0.0;
}

AxialSpHD::~AxialSpHD()
{

}

int
AxialSpHD::setTrialStrain(double strain, double strainRate)
{
  //
  trialDeformation = strain;

  // Stg index
  switch(trialStg) {

  // 1: compression, elastic
  case 1:

    // after compressive yield point
    if ( trialDeformation < ucy )
      trialStg = 5;

    // compression, elastic
    else if ( trialDeformation >= ucy && trialDeformation <= 0.0 )
      trialStg = 1;

    // tension, elastic
    else if ( trialDeformation > 0.0 && trialDeformation <= uty )
      trialStg = 2;

    // after tensile yield point
    else if ( trialDeformation > uty && trialDeformation <= uth )
      trialStg = 3;

    // after tensile hardening point
    else if ( trialDeformation > uth )
      trialStg = 4;

    break;

  // 2: tension, elastic
  case 2:

    // after compressive yield point
    if ( trialDeformation < ucy )
      trialStg = 5;

    // compression, elastic
    else if ( trialDeformation >= ucy && trialDeformation <= 0.0 )
      trialStg = 1;

    // tension, elastic
    else if ( trialDeformation > 0.0 && trialDeformation <= uty )
      trialStg = 2;

    // after tensile yield point
    else if ( trialDeformation > uty && trialDeformation <= uth )
      trialStg = 3;

    // after tensile hardening point
    else if ( trialDeformation > uth )
      trialStg = 4;

    break;

  // 3:  after tensile yield point (skeleton curve)
  case 3:

    // tension
    if ( trialDeformation >= commitDeformation ) {
      if ( trialDeformation <= uty )
	trialStg = 2; //  tension, elastic
      else if ( trialDeformation > uty && trialDeformation <= uth )
	trialStg = 3; // after tensile yield point
      else if ( trialDeformation > uth )
	trialStg = 4; // after tensile hardening point
    }

    // compression (unload)
    else {
      ur1 = commitDeformation; // strain of unload point
      fr1 = commitForce; // stress of unload point

      // trialStg=8
      if ( ur1 <= utr ) {
	ur7 = (-ste*ur1+fr1)/(sce-ste);
	fr7 = sce*ur7;

	if ( trialDeformation >= ur7 )
	  trialStg = 8;
	else if ( trialDeformation >= ucy && trialDeformation < ur7 )
	  trialStg = 1;
	else if ( trialDeformation < ucy )
	  trialStg = 5;
      }

      // trialStg=6
      else {
	ur2 = (ste*ur1-fr1)/(ste-sty);
	fr2 = sty*ur2;

	if ( trialDeformation >= ur2 )
	  trialStg = 6;
	else if ( trialDeformation >= ucr && trialDeformation < ur2 )
	  trialStg = 7;
	else if ( trialDeformation >= ucy && trialDeformation < ucr )
	  trialStg = 1;
	else if ( trialDeformation < ucy )
	  trialStg = 5;
      }
    }

    break;

  // 4: after tensile hardening point (skeleton curve)
  case 4:

    // tension
    if ( trialDeformation >= commitDeformation )
      trialStg = 4;

    // unload
    else {
      ur3 = commitDeformation; // strain of unload point
      fr3 = commitForce; // stress of unload point

      ur4 = (ste*ur3-fr3)/(ste-sty);
      fr4 = sty*ur4;

      stp = (fcr-fr4)/(ucr-ur4);
      uch = (fth-fcr-ste*uth+stp*ucr)/(stp-ste);
      fch = fcr + stp*(uch-ucr);
      ur2 = uch;
      fr2 = fch;

      if ( trialDeformation > ur4 )
	trialStg = 9;
      else if ( trialDeformation > uch && trialDeformation <= ur4 )
	trialStg = 13;
      else if ( trialDeformation > ucr && trialDeformation <= uch )
	trialStg = 7;
      else if ( trialDeformation >= ucy && trialDeformation < ucr )
	trialStg = 1;
      else if ( trialDeformation < ucy )
	trialStg = 5;
    }

    break;


  // 5: after compressive yield point (skeleton curve)
  case 5:

    // compression
    if ( trialDeformation <= commitDeformation )
      trialStg = 5;

    // tension (unload)
    else {
      ur5 = commitDeformation; // strain of unload point
      fr5 = commitForce; // stress of unload point

      uc0 = -fr5/sce + ur5;
      ur6 = (ste*uc0-sty*uty+fty)/(ste-sty);
      fr6 = ste*(ur6-uc0);

      if ( trialDeformation <= uc0 )
	trialStg = 10;
      else if ( trialDeformation > uc0 && trialDeformation <= ur6 )
	trialStg = 11;
      else if ( trialDeformation > ur6 && trialDeformation <= uty )
	trialStg = 12;
      else if ( trialDeformation > uty && trialDeformation <= uth )
	trialStg = 3;
      else if ( trialDeformation > uth )
	trialStg = 4;
    }

    break;


  // 6: after tensile yield point (unload)
  case 6:

    if ( trialDeformation > uth )
      trialStg = 4;
    else if ( trialDeformation > ur1 && trialDeformation <= uth )
      trialStg = 3;
    else if ( trialDeformation > ur2 && trialDeformation <= ur1 )
      trialStg = 6;
    else if ( trialDeformation > ucr && trialDeformation <= ur2 )
      trialStg = 7;
    else if ( trialDeformation > ucy && trialDeformation <= ucr )
      trialStg = 1;
    else if ( trialDeformation <= ucy )
      trialStg = 5;

    break;


  // 7: after tensile yield point (unload, approach reversal point)
  case 7:

    // compression
    if ( trialDeformation <= commitDeformation ) {
      if ( trialDeformation > ucr )
	trialStg = 7;
      else if ( trialDeformation > ucy && trialDeformation <= ucr )
	trialStg = 1;
      else if ( trialDeformation <= ucy )
	trialStg = 5;
    }

    // tensioniunloadj
    else {
      ur2 = commitDeformation; // strain of unload point
      fr2 = commitForce; // stress of unload point
      ur1 = (ste*ur2-sty*uty+fty-fr2)/(ste-sty);
      fr1 = sty*(ur1-uty) + fty;

      if ( trialDeformation > uth )
	trialStg = 4;
      else if ( trialDeformation > ur1 && trialDeformation <= uth )
	trialStg = 3;
      else if ( trialDeformation <= ur1 )
	trialStg = 6;
    }

    break;


  // 8: after tensile yield point (unload, strain<=utr)
  case 8:
    if ( trialDeformation > uth )
      trialStg = 4;
    else if ( trialDeformation > ur1 && trialDeformation <= uth )
      trialStg = 3;
    else if ( trialDeformation > ur7 && trialDeformation <= ur1 )
      trialStg = 8;
    else if ( trialDeformation > ucy && trialDeformation <= ur7 )
      trialStg = 1;
    else if ( trialDeformation <= ucy )
      trialStg = 5;

    break;


  // 9: after tensile hardening point (unload)
  case 9:
    if ( trialDeformation > ur3 )
      trialStg = 4;
    else if ( trialDeformation > ur4 && trialDeformation <= ur3 )
      trialStg = 9;
    else if ( trialDeformation > uch && trialDeformation <= ur4 )
      trialStg = 13;
    else if ( trialDeformation > ucr && trialDeformation <= uch )
      trialStg = 7;
    else if ( trialDeformation > ucy && trialDeformation <= ucr )
      trialStg = 1;
    else if ( trialDeformation <= ucy )
      trialStg = 5;

    break;


  // 10: after compressive yield point (unload)
  case 10:

    if ( trialDeformation <= ur5 )
      trialStg = 5;
    else if ( trialDeformation > ur5 && trialDeformation <= uc0 )
      trialStg = 10;
    else if ( trialDeformation > uc0 && trialDeformation <= ur6 )
      trialStg = 11;
    else if ( trialDeformation > ur6 && trialDeformation <= uty )
      trialStg = 12;
    else if ( trialDeformation > uty && trialDeformation <= uth )
      trialStg = 3;
    else if ( trialDeformation > uth )
      trialStg = 4;

    break;

  // 11: after compressive yield point (unload, turn into tensile stress)
  case 11:

    if ( trialDeformation <= ur5 )
      trialStg = 5;
    else if ( trialDeformation > ur5 && trialDeformation <= uc0 )
      trialStg = 10;
    else if ( trialDeformation > uc0 && trialDeformation <= ur6 )
      trialStg = 11;
    else if ( trialDeformation > ur6 && trialDeformation <= uty )
      trialStg = 12;
    else if ( trialDeformation > uty && trialDeformation <= uth )
      trialStg = 3;
    else if ( trialDeformation > uty )
      trialStg = 4;

    break;

  // 12: after compressive yield point (unload, turn into tensile stress, ur6<=strain)
  case 12:

    // tension
    if ( trialDeformation >= commitDeformation ) {
      if ( trialDeformation <= uty )
	trialStg = 12;
      else if ( trialDeformation > uty && trialDeformation <= uth )
	trialStg = 3;
      else if ( trialDeformation > uth )
	trialStg = 4;
    }

    // compression (unload)
    else {
      ur6 = commitDeformation; // strain of unload point
      fr6 = commitForce; // stress of unload point
      uc0 = -fr6/ste + ur6;
      ur5 = (sce*uc0-scy*ucy+fcy)/(sce-scy);
      fr5 = sce*(ur5-uc0);

      if ( trialDeformation <= ur5 )
	trialStg = 5;
      else if ( trialDeformation > ur5 && trialDeformation <= uc0 )
	trialStg = 10;
      else if ( trialDeformation > uc0 )
	trialStg = 11;
    }

    break;


  // 13: after tensile hardening point (unload, turn into compressive stress)
  case 13:
    // compression
    if ( trialDeformation <= commitDeformation ) {
      if ( trialDeformation > uch )
	trialStg = 13;
      else if ( trialDeformation > ucr && trialDeformation <= uch )
	trialStg = 7;
      else if ( trialDeformation > ucy && trialDeformation <= ucr )
	trialStg = 1;
      else if ( trialDeformation <= ucy )
	trialStg = 5;
    }

    // tension (unload)
    else {
      ur4 = commitDeformation; // strain of unload point
      fr4 = commitForce; // stress of unload point
      ur3 = (ste*ur4-sth*uth+fth-fr4)/(ste-sth);
      fr3 = sth*(ur3-uth) + fth;

      if ( trialDeformation > ur3 )
	trialStg = 4;
      else
	trialStg = 9;
    }

    break;

  }


  // stiffness,force
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

  // 3:  after tensile yield point (skeleton curve)
  case 3:

    trialStiffness = sty;
    trialForce = sty*(trialDeformation-uty) + fty;

    break;

  // 4: after tensile hardening point (skeleton curve)
  case 4:

    trialStiffness = sth;
    trialForce = sth*(trialDeformation-uth) + fth;

    break;

  // 5: after compressive yield point (skeleton curve)
  case 5:

    trialStiffness = scy;
    trialForce = scy*(trialDeformation-ucy) + fcy;

    break;

  // 6: after tensile yield point (unload)
  case 6:

    trialStiffness = ste;
    trialForce = ste*(trialDeformation-ur1) + fr1;

    break;

  // 7: after tensile yield point (unload, approach reversal point)
  case 7:

    trialStiffness = (fcr-fr2)/(ucr-ur2);
    trialForce = trialStiffness*(trialDeformation-ucr) + fcr;

    break;

  // 8: after tensile yield point (unload, strain<=utr)
  case 8:

    trialStiffness = ste;
    trialForce = ste*(trialDeformation-ur7) + fr7;

    break;

  // 9: after tensile hardening point (unload)
  case 9:

    trialStiffness = ste;
    trialForce = ste*(trialDeformation-ur3) + fr3;

    break;

  // 10: after compressive yield point (unload)
  case 10:

    trialStiffness = sce;
    trialForce = sce*(trialDeformation-ur5) + fr5;

    break;

  // 11: after compressive yield point (unload, turn into tensile stress)
  case 11:

    trialStiffness = ste;
    trialForce = ste*(trialDeformation-uc0);

    break;

  // 12: after compressive yield point (unload, turn into tensile stress, ur6<=strain)
  case 12:

    trialStiffness = sty;
    trialForce = sty*(trialDeformation-ur6) + fr6;

    break;

  // 13: after tensile hardening point (unload, turn into compressive stress)
  case 13:

    trialStiffness = (fcr-fr4)/(ucr-ur4);
    trialForce = trialStiffness*(trialDeformation-ucr) + fcr;

    break;

  }

  return 0;
}

double
AxialSpHD::getStress(void)
{
  return trialForce;
}

double
AxialSpHD::getTangent(void)
{
  return trialStiffness;
}

double
AxialSpHD::getInitialTangent(void)
{
  return sce;
}

double
AxialSpHD::getStrain(void)
{
  return trialDeformation;
}

int
AxialSpHD::commitState(void)
{
  commitDeformation  = trialDeformation;
  commitForce  = trialForce;
  commitStiffness = trialStiffness;

  commitStg = trialStg;

  return 0;
}

int
AxialSpHD::revertToLastCommit(void)
{
  trialDeformation  = commitDeformation;
  trialForce  = commitForce;
  trialStiffness = commitStiffness;

  trialStg = commitStg;

  return 0;
}

int
AxialSpHD::revertToStart(void)
{
  trialDeformation = 0.0;
  trialForce = 0.0;
  trialStiffness = sce;
  commitDeformation = 0.0;
  commitForce = 0.0;
  commitStiffness = sce;

  trialStg = 1;
  commitStg = 1;

  ste  = bte*sce;
  sty  = bty*sce;
  sth  = bth*sce;
  scy  = bcy*sce;
  uty  = fty/ste;
  ucy  = fcy/sce;
  ucr  = fcr/sce;
  utr  = (ste*ucr-sty*uty+fty-fcr)/(ste-sty);
  ftr  = sty*(utr-uty) + fty;
  uth  = ath*uty;
  fth  = sty*(uth-uty) + fty;
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
  ur6  = 0.0;
  fr6  = 0.0;
  ur7  = 0.0;
  fr7  = 0.0;
  
  return 0;
}

UniaxialMaterial *
AxialSpHD::getCopy(void)
{
  AxialSpHD *theCopy = new AxialSpHD(this->getTag(), sce, fty, fcy, bte,
				 bty, bth, bcy, fcr, ath);

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
AxialSpHD::sendSelf(int cTag, Channel &theChannel)
{
  int res = 0;

  static Vector data(18+15);

  data(0)  = this->getTag();
  data(1)  = sce;
  data(2)  = fty;
  data(3)  = fcy;
  data(4)  = bte;
  data(5)  = bty;
  data(6)  = bth;
  data(7)  = bcy;
  data(8)  = fcr;
  data(9)  = ath;
  data(10) = commitDeformation;
  data(11) = commitForce;
  data(12) = commitStiffness;
  data(13) = commitStg;
  data(14) = trialDeformation;
  data(15) = trialForce;
  data(16) = trialStiffness;
  data(17) = trialStg;

  data(18) = uc0;
  data(19) = ur1;
  data(20) = fr1;
  data(21) = ur2;
  data(22) = fr2;
  data(23) = ur3;
  data(24) = fr3;
  data(25) = ur4;
  data(26) = fr4;
  data(27) = ur5;
  data(28) = fr5;
  data(29) = ur6;
  data(30) = fr6;
  data(31) = ur7;
  data(32) = fr7;  
  
  res = theChannel.sendVector(this->getDbTag(), cTag, data);
  if ( res < 0 )
    opserr << "AxialSpHD::sendSelf() - failed to send data\n";

  return res;
}

int
AxialSpHD::recvSelf(int cTag, Channel &theChannel,
		  FEM_ObjectBroker &theBroker)
{
  int res = 0;

  static Vector data(18+15);
  res = theChannel.recvVector(this->getDbTag(), cTag, data);

  if ( res < 0 ) {
    opserr << "AxialSpHD::recvSelf() - failed to receive data\n";
    this->setTag(0);
  }

  else {
    this->setTag((int)data(0));
    sce = data(1);
    fty = data(2);
    fcy = data(3);
    bte = data(4);
    bty = data(5);
    bth = data(6);
    bcy = data(7);
    fcr = data(8);
    ath = data(9);
    commitDeformation = data(10);
    commitForce = data(11);
    commitStiffness = data(12);
    commitStg = (int)data(13);
    trialDeformation = data(14);
    trialForce = data(15);
    trialStiffness = data(16);
    trialStg = (int)data(17);

    ste  = bte*sce;
    sty  = bty*sce;
    sth  = bth*sce;
    scy  = bcy*sce;
    uty  = fty/ste;
    ucy  = fcy/sce;
    ucr  = fcr/sce;
    utr  = (ste*ucr-sty*uty+fty-fcr)/(ste-sty);
    ftr  = sty*(utr-uty) + fty;
    uth  = ath*uty;
    fth  = sty*(uth-uty) + fty;

    uc0 = data(18);
    ur1 = data(19);
    fr1 = data(20);
    ur2 = data(21);
    fr2 = data(22);
    ur3 = data(23);
    fr3 = data(24);
    ur4 = data(25);
    fr4 = data(26);
    ur5 = data(27);
    fr5 = data(28);
    ur6 = data(29);
    fr6 = data(30);
    ur7 = data(31);
    fr7 = data(32);
  }

  return res;
}

void
AxialSpHD::Print(OPS_Stream &s, int flag)
{
  s << "AxialSpHD : " << this->getTag() << endln;
  s << "    sce : " << endln;
  s << "    fty : " << endln;
  s << "    fcy : " << endln;
  s << "    bte : " << endln;
  s << "    bty : " << endln;
  s << "    bth : " << endln;
  s << "    bcy : " << endln;
  s << "    fcr : " << endln;
  s << "    ath : " << endln;

  return;

}

