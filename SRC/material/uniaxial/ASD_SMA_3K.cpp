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
// $Date: 2020-09-01 11:29:01 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/ASD_SMA_3K.cpp,v $

// Written: Luca Aceto
// Created: Aug
// Revision: A
//
// Description: This file contains the class implementation for ASD_SMA_3K. 

/*
ASD_SMA_3K is written by Eng. Luca Aceto (luca.aceto@unich.it), University of Chieti-Pescara, InGeo department in collaboration with ASDEA Software Technology: https://asdeasoft.net


This material is a modified version of Self Centering Material written by JAE at Oct 2007
With ASD_SMA_3K it is possible to replicate the behavior of SMA material with a different unloading stiffness (k3)

        k1 = Load stiffness
        k2 = Post-Activation Stiffness (0<$k2<$k1)
        k3 = Unload stiffness
        sigAct = Forward Activation Stress/Force
        beta= Ratio of Forward to Reverse Activation Stress/Force


ASD_SMA_3K matTag? k1? k2? k3? sigF? beta?

*/


#include <ASD_SMA_3K.h>
#include <Vector.h>
#include <Channel.h>
#include <Matrix.h>
#include <Information.h>
#include <Parameter.h>

#include <cmath>
#include <float.h>
#include <elementAPI.h>
#include <iostream>

/* static class instance counter */
static int ASD_SMA3K_counter = 0;

void* OPS_ASD_SMA_3K()
{



    int numdata = OPS_GetNumRemainingInputArgs();
    if (numdata < 5) {
	opserr << "WARNING: Insufficient arguments\n";
	opserr << "Want: uniaxialMaterial ASD_SMA_3K matTag? k1? k2? k3? sigF? beta?";
	return 0;
    }

    int tag;
    numdata = 1;
    if (OPS_GetIntInput(&numdata,&tag) < 0) {
	opserr << "WARNING invalid tag\n";
	return 0;
    }

    double data[5] = {0,0,0,0,0};
    numdata = OPS_GetNumRemainingInputArgs();
    if (numdata > 5) {
	numdata = 5;
    }
    if (OPS_GetDoubleInput(&numdata,data)) {
	opserr << "WARNING invalid double inputs\n";
	return 0;
    }

    UniaxialMaterial* mat = new ASD_SMA_3K(tag,data[0],data[1],data[2],data[3],data[4]);
    if (mat == 0) {
	opserr << "WARNING: failed to create ASD_SMA_3K material\n";
	return 0;
    }

    return mat;
}

ASD_SMA_3K::ASD_SMA_3K(int tag, double K1, double K2, double K3,
					     double fa, double b)
  : UniaxialMaterial(tag,MAT_TAG_ASD_SMA_3K),
    k1(K1), k2(K2), k3(K3), ActF(fa), beta(b)
{
    // on first call
    if (ASD_SMA3K_counter == 0) {
        ASD_SMA3K_counter++;
        opserr <<
            "\n*******************************************************************************************\n"
            "* ASD_SMA_3K - Written by Eng. Luca Aceto, University of Chieti-Pescara, InGeo department *\n"
            "* in collaboration with ASDEA Software Technology                                         *\n"
            "* Eng. Luca Aceto luca.aceto@unich.it                                                     *\n"
            "* ASDEA Software Technology: https://asdeasoft.net                                        *\n"
            "* STKO (Scientific ToolKit for OpenSees): https://asdeasoft.net/stko/                     *\n"
            "*******************************************************************************************\n"
            ;
    }
  // Find Equivalent Slip Force
  ActDef = ActF / k1;

  
  // Initialize variables
  this->revertToStart();
  
}

ASD_SMA_3K::ASD_SMA_3K()
  : UniaxialMaterial(0,MAT_TAG_ASD_SMA_3K),
    k1(0.0), k2(0.0), k3(0.0), ActF(0.0), beta(0.0)
{
  // Initialize variables
  ActDef = 0;
  
  this->revertToStart();
  
}

ASD_SMA_3K::~ASD_SMA_3K()
{
  
}

int 

ASD_SMA_3K::setTrialStrain(double strain, double strainRate)
{


    diffStrain = strain - Cstrain;

    if (fabs(diffStrain) < DBL_EPSILON)
        return 0;

    // Set total strain
    Tstrain = strain;
    // Tstrain = Tstrain - CslipStrain;

    // Middle Elastic Portion (outside any upper or lower activation)
    //     Entirely elastic response
    if (fabs(Tstrain) <= ((1 - beta) * ActF / k1)) {

        Tstress = k1 * Tstrain;
        Ttangent = k1;

        

        TupperStrainPos = ActDef;
        TupperStressPos = ActF;

        TupperStrainNeg = -ActDef;
        TupperStressNeg = -ActF;

        No_Y_Pos = 0;
        No_Y_Neg = 0;

    }
    else {

        // Positive Quadrant (Top Right) where strain >= 0
        if (Tstrain >= 0) {

            TupperStrainNeg = -ActDef;
            TupperStressNeg = -ActF;

                 // Linear range movement (no upper or lower activation)  
            if ((Tstrain >= ClowerStrainPos) &&
                (Tstrain <= CupperStrainPos)) {

                if (diffStrain < DBL_EPSILON) {
                    if (No_Y_Pos == 1) {
                        Tstress = Cstress + diffStrain * k3;
                        Ttangent = k3;

                        TupperStrainPos = (k1 * Tstrain - Tstress - k2 * ActDef + ActF) / (k1 - k2);
                        TupperStressPos = k2 * TupperStrainPos - k2 * ActDef + ActF;
                    }
                    else {
                        Tstress = fmin(Cstress + diffStrain * k1, Tstrain * k1);
                        Ttangent = k1;
                    }
                }
                else if (diffStrain > DBL_EPSILON) {

                    if (No_Y_Pos == 1) {
                        Tstress = fmin(Cstress + diffStrain * k1,Tstrain * k1);
                        Ttangent = k1;
                        No_k2_Pos = 0;

                        TupperStrainPos = (k1 * Tstrain - Tstress - k2 * ActDef + ActF) / (k1 - k2);
                        TupperStressPos = k2 * TupperStrainPos - k2 * ActDef + ActF;

                    } else if (No_Y_Pos == 0) {

                        Tstress = fmin((Tstrain - CactivStrainPos) * k1,Tstrain*k1);
                        Ttangent = k1;
                        No_k2_Pos = 0;
                    }


                }

            }
            // Upper Activation 
            else if (Tstrain > CupperStrainPos) {

                No_Y_Pos = 1;

                TupperStressPos = CupperStressPos + (Tstrain - CupperStrainPos) * k2;
                TupperStrainPos = Tstrain;

                X = (-k3 * TupperStrainPos + TupperStressPos) / (k1 - k3);
                Y = k1 * ((-k3 * TupperStrainPos + TupperStressPos) / (k1 - k3));

                TlowerStrainPos = (TupperStressPos - TupperStrainPos * k3 - ActF * (1 - beta) * (1 - (k2 / k1))) / (k2 - k3);
                TlowerStressPos = k2 * ((TupperStressPos - TupperStrainPos * k3 - ActF * (1 - beta) * (1 - (k2 / k1))) / (k2 - k3)) + ActF * (1 - beta) * (1 - (k2 / k1));
                No_k2_Pos = 0;

                if (TlowerStrainPos < X) {
                    TlowerStrainPos = X;
                    TlowerStressPos = Y;
                    No_k2_Pos = 1;
                }

                Tstress = TupperStressPos;
                TactivStrainPos = TupperStrainPos - Tstress / k3;

                Ttangent = k2;


            }

            else { // Tstrain < ClowerStrainPos

                if (k1 * Tstrain <= ClowerStressPos && No_k2_Pos == 1) {
                    Tstress = fmin((Tstrain - CactivStrainPos) * k1, Tstrain * k1);
                    Ttangent = k1;

                    TupperStrainPos = ActDef;
                    TupperStressPos = ActF;

                    TactivStrainPos = Tstrain - Tstress / k1;


                }
                else {
                    TlowerStressPos = ClowerStressPos +
                        (Tstrain - ClowerStrainPos) * k2;
                    TlowerStrainPos = Tstrain;
                    TupperStrainPos = Tstrain + beta * ActF / k1;
                    TupperStressPos = TlowerStressPos + beta * ActF;
                    Tstress = TlowerStressPos;
                    TactivStrainPos = TlowerStrainPos - Tstress / k1;

                    Ttangent = k2;
                }
            }
        }

        // Negative Quadrant (Bottom Left) where strain < 0
        else { // Tstrain < 0)

            TupperStrainPos = ActDef;
            TupperStressPos = ActF;

            // Linear range movement (no upper or
            //     lower activation)
            if ((Tstrain <= ClowerStrainNeg) &&
                (Tstrain >= CupperStrainNeg)) {

                if (diffStrain > DBL_EPSILON) {
                    if (No_Y_Neg == 1) {
                        Tstress = Cstress+ diffStrain * k3;
                        Ttangent = k3;

                        TupperStrainNeg = (k1 * Tstrain - Tstress + k2 * ActDef - ActF) / (k1 - k2);
                        TupperStressNeg = k2 * TupperStrainNeg + k2 * ActDef - ActF;
                    }
                    else {
                        Tstress = fmax(Cstress + diffStrain * k1, Tstrain * k1);
                        Ttangent = k1;
                    }
                }
                else if (diffStrain < DBL_EPSILON) {

                    if (No_Y_Neg == 1) {
                        Tstress = fmax(Cstress + diffStrain * k1, Tstrain * k1);
                        Ttangent = k1;
                        No_k2_Neg = 0;

                        TupperStrainNeg = (k1 * Tstrain - Tstress + k2 * ActDef - ActF) / (k1 - k2);
                        TupperStressNeg = k2 * TupperStrainNeg + k2 * ActDef - ActF;

                    }
                    else if (No_Y_Neg == 0) {
                        Tstress = fmax((Tstrain - CactivStrainNeg) * k1, Tstrain * k1);
                        Ttangent = k1;
                        No_k2_Neg = 0;


                    }
                }

            }
            // Upper Activation
            else if (Tstrain < CupperStrainNeg) {

                No_Y_Neg = 1;

                TupperStressNeg = CupperStressNeg + (Tstrain - CupperStrainNeg) * k2;
                TupperStrainNeg = Tstrain;

                X = (-k3 * TupperStrainNeg + TupperStressNeg) / (k1 - k3);
                Y = k1 * ((-k3 * TupperStrainNeg + TupperStressNeg) / (k1 - k3));

                TlowerStrainNeg = (TupperStressNeg - TupperStrainNeg * k3 + ActF * (1 - beta) * (1 - (k2 / k1))) / (k2 - k3);
                TlowerStressNeg = k2 * ((TupperStressNeg - TupperStrainNeg * k3 + ActF * (1 - beta) * (1 - (k2 / k1))) / (k2 - k3)) - ActF * (1 - beta) * (1 - (k2 / k1));


                No_k2_Neg = 0;

                if (TlowerStrainNeg > X) {
                    TlowerStrainNeg = X;
                    TlowerStressNeg = Y;
                    No_k2_Neg = 1;
                }

                Tstress = TupperStressNeg;
                TactivStrainNeg = TupperStrainNeg - Tstress / k3;

                Ttangent = k2;
            }
            // Lower Activation
            else { // Tstrain > ClowerStrainNeg


                if (k1 * Tstrain >= ClowerStressNeg && No_k2_Neg == 1) {
                    Tstress = fmax((Tstrain - CactivStrainNeg) * k1, Tstrain * k1);
                    Ttangent = k1;

                    TupperStrainNeg = -ActDef;
                    TupperStressNeg = -ActF;

                    TactivStrainNeg = Tstrain - Tstress / k1;


                }
                else {
                    TlowerStressNeg = ClowerStressNeg + (Tstrain - ClowerStrainNeg) * k2;
                    TlowerStrainNeg = Tstrain;
                    TupperStrainNeg = Tstrain - beta * ActF / k1;
                    TupperStressNeg = TlowerStressNeg - beta * ActF;
                    Tstress = TlowerStressNeg;
                    TactivStrainNeg = TlowerStrainNeg - Tstress / k1;

                    Ttangent = k2;
                }



            }
        }

    }



    return 0;
}

double 
ASD_SMA_3K::getStress(void)
{
  return Tstress;
}

double 
ASD_SMA_3K::getTangent(void)
{
  return Ttangent;
}

double 
ASD_SMA_3K::getStrain(void)
{
  return Tstrain;
}

int 
ASD_SMA_3K::commitState(void)
{

      // Commit trial history variables
      CactivStrainPos = TactivStrainPos;
      CactivStrainNeg = TactivStrainNeg;
      CupperStrainPos = TupperStrainPos;
      ClowerStrainPos = TlowerStrainPos;	
      CupperStressPos = TupperStressPos;
      ClowerStressPos = TlowerStressPos;
      CupperStrainNeg = TupperStrainNeg;
      ClowerStrainNeg = TlowerStrainNeg;
      CupperStressNeg = TupperStressNeg;
      ClowerStressNeg = TlowerStressNeg;
      Cstrain = Tstrain;
      Cstress = Tstress;
      Ctangent = Ttangent;

      CLastStrain = TLastStrain;
  
      return 0;
}

int 
ASD_SMA_3K::revertToLastCommit(void)
{
  Tstrain = Cstrain;
  Tstress = Cstress;
  Ttangent = Ctangent;

  return 0;
}

int 
ASD_SMA_3K::revertToStart(void)
{
  // Reset committed history variables
  CactivStrainPos = 0.0;
  CactivStrainNeg = 0.0;
  CupperStrainPos = ActDef;
  ClowerStrainPos = (1-beta) * ActDef;	
  CupperStressPos = ActF;
  ClowerStrainPos = (1 - beta) * ActDef;
  CupperStrainNeg = -CupperStrainPos;
  ClowerStrainNeg = -ClowerStrainPos;
  CupperStressNeg = -CupperStressPos;
  ClowerStressNeg = -ClowerStressPos;

  CLastStrain = 0.0;
  
  // Reset trial history variables
  TactivStrainPos = 0.0;
  TactivStrainNeg = 0.0;
  TupperStrainPos = ActDef;
  TlowerStrainPos = (1-beta) * ActDef;	
  TupperStressPos = ActF;
  TlowerStressPos = (1-beta) * ActF;
  TupperStrainNeg = -CupperStrainPos;
  TlowerStrainNeg = -ClowerStrainPos;
  TupperStressNeg = -CupperStressPos;
  TlowerStressNeg = -ClowerStressPos;

  TLastStrain = 0.0;
  
  // Initialize state variables
  Tstrain = 0.0;
  Tstress = 0.0;
  Ttangent = k1;
  
  Cstrain = 0.0;

  No_k2_Pos = 0;
  No_k2_Neg = 0;
  No_Y_Pos = 0;
  No_Y_Neg = 0;
  
  return 0;
}

UniaxialMaterial *
ASD_SMA_3K::getCopy(void)
{
  ASD_SMA_3K *theCopy =
    new ASD_SMA_3K(this->getTag(), k1, k2, k3, ActF, beta);
  
  // Copy committed history variables
  theCopy->CactivStrainPos = CactivStrainPos;
  theCopy->CactivStrainNeg = CactivStrainNeg;
  theCopy->CupperStrainPos = CupperStrainPos;
  theCopy->ClowerStrainPos = ClowerStrainPos;	
  theCopy->CupperStressPos = CupperStressPos;
  theCopy->ClowerStressPos = ClowerStressPos;
  theCopy->CupperStrainNeg = CupperStrainNeg;
  theCopy->ClowerStrainNeg = ClowerStrainNeg;
  theCopy->CupperStressNeg = CupperStressNeg;
  theCopy->ClowerStressNeg = ClowerStressNeg;

  theCopy->CLastStrain = CLastStrain;
  
  // Copy trial history variables
  theCopy->TactivStrainPos = TactivStrainPos;
  theCopy->TactivStrainNeg = TactivStrainNeg;
  theCopy->TupperStrainPos = TupperStrainPos;
  theCopy->TlowerStrainPos = TlowerStrainPos;	
  theCopy->TupperStressPos = TupperStressPos;
  theCopy->TlowerStressPos = TlowerStressPos;
  theCopy->TupperStrainNeg = TupperStrainNeg;
  theCopy->TlowerStrainNeg = TlowerStrainNeg;
  theCopy->TupperStressNeg = TupperStressNeg;
  theCopy->TlowerStressNeg = TlowerStressNeg;

  theCopy->TLastStrain = TLastStrain;
  
  // Copy trial state variables
  theCopy->Tstrain = Tstrain;
  theCopy->Tstress = Tstress;
  theCopy->Ttangent = Ttangent;
  
  theCopy->Cstrain = Cstrain;
  
  return theCopy;
}

int 
ASD_SMA_3K::sendSelf(int cTag, Channel &theChannel)
{
  int res = 0;
  
  static Vector data(23);
  
  data(0) = this->getTag();
  data(1) = k1;
  data(2) = k2;
  data(3) = k3;
  data(4) = ActF;
  data(5) = beta;
  data(6) = ActDef;
  data(7) = CactivStrainPos;
  data(8) = CactivStrainNeg;
  data(9) = CupperStrainPos;
  data(10) = ClowerStrainPos;
  data(11) = CupperStressPos;
  data(12) = ClowerStressPos;
  data(13) = CupperStrainNeg;
  data(14) = ClowerStrainNeg;
  data(15) = CupperStressNeg;
  data(16) = ClowerStressNeg;
  data(17) = Tstrain;
  data(18) = Tstress;
  data(19) = Ttangent;
  data(20) = Cstrain;
  data(21) = TLastStrain;
  data(22) = CLastStrain;
  
  res = theChannel.sendVector(this->getDbTag(), cTag, data);
  if (res < 0) 
    opserr << "ASD_SMA_3K::sendSelf() - failed to send data\n";

  return res;
}

int 
ASD_SMA_3K::recvSelf(int cTag, Channel &theChannel,
			       FEM_ObjectBroker &theBroker)
{
  int res = 0;
  
  static Vector data(23);
  res = theChannel.recvVector(this->getDbTag(), cTag, data);
  
  if (res < 0) {
      opserr << "ASD_SMA_3K::recvSelf() - failed to receive data\n";
      this->setTag(0);      
  }
  else {
    this->setTag((int)data(0));
    k1 = data(1);
	k2 = data(2);
    k3 = data(3);
	ActF = data(4);
	beta = data(5);
	ActDef = data(6);
	CactivStrainPos = data(7);
	CactivStrainNeg = data(8);
	CupperStrainPos = data(9);
	ClowerStrainPos = data(10);
	CupperStressPos = data(11);
	ClowerStressPos = data(12);
	CupperStrainNeg = data(13);
	ClowerStrainNeg = data(14);
	CupperStressNeg = data(15);
	ClowerStressNeg = data(16);
    Tstrain = data(17);
	Tstress = data(18);
	Ttangent = data(19);
	Cstrain = data(20);
    TLastStrain = data(21) ;
    CLastStrain=data(22) ;
  }
    
  return res;
}

void 
ASD_SMA_3K::Print(OPS_Stream &s, int flag)
{
    s << "ASD_SMA_3K, tag: " << this->getTag() << endln;
    s << "  k1: " << k1 << endln;
    s << "  k2: " << k2 << endln;
    s << "  k3: " << k3 << endln;
    s << "  ActF: " << ActF << endln;
    s << "  beta: " << beta << endln;
}

