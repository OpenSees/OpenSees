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
// $Date: 2005-11-09 00:30:14 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/Steel03.cpp,v $
                                                                        
// Written: mackie
// Created: 06/2005
// Revision: A
//
// Description: This file contains the class implementation for 
// Steel03. Steel03 is Steel01 verbatim but with added Giuffre-Menegotto-Pinto 
// transitions on the loading and unloading loops.  
// references:
// 1.) 	Menegotto, M., and Pinto, P.E. (1973). Method of analysis of cyclically loaded 
//	RC plane frames including changes in geometry and non-elastic behavior of 
//	elements under normal force and bending. Preliminary Report IABSE, vol 13. 
// 2.)	Dhakal, R.J., and Maekawa, K. (2002). Path-dependent cyclic stress-strain relationship
//	of reinforcing bar including buckling. Engineering Structures, 24(11): 1383-96.
// 3.)	Gomes, A., and Appleton, J. (1997). Nonlinear cyclic stress-strain relationship of 
//	reinforcing bars including buckling. Engineering Structures, 19(10): 822-6.
//
// What: "@(#) Steel03.C, revA"


#include <Steel03.h>
#include <Vector.h>
#include <Matrix.h>
#include <Channel.h>
#include <Information.h>
#include <math.h>
#include <float.h>

Steel03::Steel03
(int tag, double FY, double E, double B, double R, int RTYP, 
 double A1, double A2, double A3, double A4):
   UniaxialMaterial(tag,MAT_TAG_Steel03),
   fy(FY), E0(E), b(B), r(R), rtype(RTYP), a1(A1), a2(A2), a3(A3), a4(A4)
{
   // Sets all history and state variables to initial values
   // History variables
   CminStrain = 0.0;
   CmaxStrain = 0.0;
   CshiftP = 1.0;
   CshiftN = 1.0;
   Cloading = 0;
   CbStrain = 0.0;
   CbStress = 0.0;
   CrStrain = 0.0;
   CrStress = 0.0;
   Cplastic = 0.0;

   TminStrain = 0.0;
   TmaxStrain = 0.0;
   TshiftP = 1.0;
   TshiftN = 1.0;
   Tloading = 0;
   TbStrain = 0.0;
   TbStress = 0.0;
   TrStrain = 0.0;
   TrStress = 0.0;
   Tplastic = 0.0;

   // State variables
   Cstrain = 0.0;
   Cstress = 0.0;
   Ctangent = E0;
   CcurR = getR(0);

   Tstrain = 0.0;
   Tstress = 0.0;
   Ttangent = E0;
   TcurR = getR(0);
}

Steel03::Steel03():UniaxialMaterial(0,MAT_TAG_Steel03),
 fy(0.0), E0(0.0), b(0.0), r(0.0), rtype(1), a1(0.0), a2(0.0), a3(0.0), a4(0.0)
{

}

Steel03::~Steel03 ()
{

}

int Steel03::setTrialStrain (double strain, double strainRate)
{
   // Reset history variables to last converged state
   TminStrain = CminStrain;
   TmaxStrain = CmaxStrain;
   TshiftP = CshiftP;
   TshiftN = CshiftN;
   Tloading = Cloading;
   TbStrain = CbStrain;
   TbStress = CbStress;
   TrStrain = CrStrain;
   TrStress = CrStress;
   Tplastic = Cplastic;
   TcurR = CcurR;

   // Set trial strain
   Tstrain = strain;

   // Determine change in strain from last converged state
   double dStrain = Tstrain - Cstrain;

   // Calculate the trial state given the trial strain
   if (fabs(dStrain) > DBL_EPSILON)   
     determineTrialState (dStrain);

   return 0;
}

int Steel03::setTrial (double strain, double &stress, double &tangent, double strainRate)
{
   // Reset history variables to last converged state
   TminStrain = CminStrain;
   TmaxStrain = CmaxStrain;
   TshiftP = CshiftP;
   TshiftN = CshiftN;
   Tloading = Cloading;
   TbStrain = CbStrain;
   TbStress = CbStress;
   TrStrain = CrStrain;
   TrStress = CrStress;
   Tplastic = Cplastic;
   TcurR = CcurR;

   // Set trial strain
   Tstrain = strain;

   // Determine change in strain from last converged state
   double dStrain;
   dStrain = Tstrain - Cstrain;

   // Calculate the trial state given the trial strain
   if (fabs(dStrain) > DBL_EPSILON) 
     determineTrialState (dStrain);

   stress = Tstress;
   tangent = Ttangent;

   return 0;
}

double Steel03::getR (double x_in)
{
    // maybe modify this later, but it gives us better degradation at smaller strains
    x_in = fabs(x_in);
    double temp_r = r;
    
    if (rtype == 2) {
        // Dhakal and Maekawa (see ref above)
        temp_r = r - 18.5*x_in/(0.15+x_in);
        if (temp_r < 0)
            temp_r = 1.0e-8;
    } else if (rtype == 3) {
        // Gomes and Appleton (see ref above)
        temp_r = r - 19.0*x_in/(0.3+x_in);
        if (temp_r < 0)
            temp_r = 1.0e-8;
    } else {
        // Mackie, rough trilinear fit to the tangent to the x_in-r first 
        // quadrant circle.  Try using with small values of R like 2 and 3
        temp_r = r;
        double t1 = -x_in/7+15/7*r;
        double t2 = -4*x_in+6*r;
        if (t1 > temp_r)
            temp_r = t1;
        if (t2 > temp_r)
            temp_r = t2;
        //opserr << "xin = " << x_in << " rout = " << temp_r << endln;
    }
    
    return temp_r;
}

void Steel03::determineTrialState (double dStrain)
{
      double fyOneMinusB = fy * (1.0 - b);

      double Esh = b*E0;
      double epsy = fy/E0;
      
      double c1 = Esh*Tstrain;
      double c2 = TshiftN*fyOneMinusB;
      double c3 = TshiftP*fyOneMinusB;
      double c = Cstress + E0*dStrain;
      
      //
      // Determine if a load reversal has occurred due to the trial strain
      //

      // Determine initial loading condition
      if (Tloading == 0 && dStrain != 0.0) {
          TmaxStrain = epsy;
          TminStrain = -epsy;
	  if (dStrain > 0.0) {
	    Tloading = 1;
            TbStrain = TmaxStrain;
            TbStress = fy;
            Tplastic = TmaxStrain;
          }
	  else {
	    Tloading = -1;
            TbStrain = TminStrain;
            TbStress = -fy;
            Tplastic = TminStrain;
          }

          double intval = 1+pow(fabs(Tstrain/epsy),TcurR);
          Tstress = c1+(1-b)*E0*Tstrain/pow(intval,1/TcurR);
          Ttangent = Esh+E0*(1-b)/pow(intval,1+1/TcurR);
      }
          
      // Transition from loading to unloading, i.e. positive strain increment
      // to negative strain increment
      if (Tloading == 1 && dStrain < 0.0) {
	  Tloading = -1;
	  if (Cstrain > TmaxStrain)
	    TmaxStrain = Cstrain;
          Tplastic = TminStrain;
	  TshiftN = 1 + a1*pow((TmaxStrain-TminStrain)/(2.0*a2*epsy),0.8);
          TrStrain = Cstrain;
          TrStress = Cstress;
          TbStrain = (c2+c)/E0/(b-1)+Tstrain/(1-b);
          TbStress = 1/(b-1)*(b*c2+b*c-c1)-c2;
          TcurR = getR((TbStrain-TminStrain)/epsy);
      }

      // Transition from unloading to loading, i.e. negative strain increment
      // to positive strain increment
      if (Tloading == -1 && dStrain > 0.0) {
	  Tloading = 1;
	  if (Cstrain < TminStrain)
	    TminStrain = Cstrain;
          Tplastic = TmaxStrain;
	  TshiftP = 1 + a3*pow((TmaxStrain-TminStrain)/(2.0*a4*epsy),0.8);
          TrStrain = Cstrain;
          TrStress = Cstress;
          TbStrain = (c3-c)/E0/(1-b)+Tstrain/(1-b);
          TbStress = 1/(1-b)*(b*c3-b*c+c1)+c3;
          TcurR = getR((TmaxStrain-TbStrain)/epsy);
      }
      
      if (Cloading != 0) {
          double c4 = TbStrain - TrStrain;
          double c5 = TbStress - TrStress;
          double c6 = Tstrain - TrStrain;
          double c4c5 = c5/c4;
          double intval = 1+pow(fabs(c6/c4),TcurR);
          
          Tstress = TrStress+b*c4c5*c6+(1-b)*c4c5*c6/pow(intval,1/TcurR);
          Ttangent = c4c5*b+c4c5*(1-b)/pow(intval,1+1/TcurR);
      }
}

double Steel03::getStrain ()
{
   return Tstrain;
}

double Steel03::getStress ()
{
   return Tstress;
}

double Steel03::getTangent ()
{
   return Ttangent;
}

int Steel03::commitState ()
{
   // History variables
   CminStrain = TminStrain;
   CmaxStrain = TmaxStrain;
   CshiftP = TshiftP;
   CshiftN = TshiftN;
   Cloading = Tloading;
   CbStrain = TbStrain;
   CbStress = TbStress;
   CrStrain = TrStrain;
   CrStress = TrStress;
   Cplastic = Tplastic;
   CcurR = TcurR;

   // State variables
   Cstrain = Tstrain;
   Cstress = Tstress;
   Ctangent = Ttangent;

   return 0;
}

int Steel03::revertToLastCommit ()
{
   // Reset trial history variables to last committed state
   TminStrain = CminStrain;
   TmaxStrain = CmaxStrain;
   TshiftP = CshiftP;
   TshiftN = CshiftN;
   Tloading = Cloading;
   TbStrain = CbStrain;
   TbStress = CbStress;
   TrStrain = CrStrain;
   TrStress = CrStress;
   Tplastic = Cplastic;

   // Reset trial state variables to last committed state
   Tstrain = Cstrain;
   Tstress = Cstress;
   Ttangent = Ctangent;
   TcurR = CcurR;

   return 0;
}

int Steel03::revertToStart ()
{
   // History variables
   CminStrain = 0.0;
   CmaxStrain = 0.0;
   CshiftP = 1.0;
   CshiftN = 1.0;
   Cloading = 0;
   CbStrain = 0.0;
   CbStress = 0.0;
   CrStrain = 0.0;
   CrStress = 0.0;
   Cplastic = 0.0;

   TminStrain = 0.0;
   TmaxStrain = 0.0;
   TshiftP = 1.0;
   TshiftN = 1.0;
   Tloading = 0;
   TbStrain = 0.0;
   TbStress = 0.0;
   TrStrain = 0.0;
   TrStress = 0.0;
   Tplastic = 0.0;

   // State variables
   Cstrain = 0.0;
   Cstress = 0.0;
   Ctangent = E0;
   CcurR = getR(0);

   Tstrain = 0.0;
   Tstress = 0.0;
   Ttangent = E0;
   TcurR = getR(0);

   return 0;
}

UniaxialMaterial* Steel03::getCopy ()
{
   Steel03* theCopy = new Steel03(this->getTag(), fy, E0, b, r, rtype, 
				  a1, a2, a3, a4);

   // Converged history variables
   theCopy->CminStrain = CminStrain;
   theCopy->CmaxStrain = CmaxStrain;
   theCopy->CshiftP = CshiftP;
   theCopy->CshiftN = CshiftN;
   theCopy->Cloading = Cloading;
   theCopy->CbStrain = CbStrain;
   theCopy->CbStress = CbStress;
   theCopy->CrStrain = CrStrain;
   theCopy->CrStress = CrStress;
   theCopy->Cplastic = Cplastic;

   // Trial history variables
   theCopy->TminStrain = TminStrain;
   theCopy->TmaxStrain = TmaxStrain;
   theCopy->TshiftP = TshiftP;
   theCopy->TshiftN = TshiftN;
   theCopy->Tloading = Tloading;
   theCopy->TbStrain = TbStrain;
   theCopy->TbStress = TbStress;
   theCopy->TrStrain = TrStrain;
   theCopy->TrStress = TrStress;
   theCopy->Tplastic = Tplastic;

   // Converged state variables
   theCopy->Cstrain = Cstrain;
   theCopy->Cstress = Cstress;
   theCopy->Ctangent = Ctangent;
   theCopy->CcurR = CcurR;

   // Trial state variables
   theCopy->Tstrain = Tstrain;
   theCopy->Tstress = Tstress;
   theCopy->Ttangent = Ttangent;
   theCopy->TcurR = TcurR;

   return theCopy;
}

int Steel03::sendSelf (int commitTag, Channel& theChannel)
{
   int res = 0;
   static Vector data(24);
   data(0) = this->getTag();

   // Material properties
   data(1) = fy;
   data(2) = E0;
   data(3) = b;
   data(4) = r;
   data(5) = rtype;
   data(6) = a1;
   data(7) = a2;
   data(8) = a3;
   data(9) = a4;

   // History variables from last converged state
   data(10) = CminStrain;
   data(11) = CmaxStrain;
   data(12) = CshiftP;
   data(13) = CshiftN;
   data(14) = Cloading;
   data(15) = CbStrain;
   data(16) = CbStress;
   data(17) = CrStrain;
   data(18) = CrStress;
   data(19) = Cplastic;

   // State variables from last converged state
   data(20) = Cstrain;
   data(21) = Cstress;
   data(22) = Ctangent;
   data(23) = CcurR;

   // Data is only sent after convergence, so no trial variables
   // need to be sent through data vector

   res = theChannel.sendVector(this->getDbTag(), commitTag, data);
   if (res < 0) 
      opserr << "Steel03::sendSelf() - failed to send data\n";

   return res;
}

int Steel03::recvSelf (int commitTag, Channel& theChannel,
                                FEM_ObjectBroker& theBroker)
{
   int res = 0;
   static Vector data(24);
   res = theChannel.recvVector(this->getDbTag(), commitTag, data);
  
   if (res < 0) {
      opserr << "Steel03::recvSelf() - failed to receive data\n";
      this->setTag(0);      
   }
   else {
      this->setTag(int(data(0)));

      // Material properties
      fy = data(1);
      E0 = data(2);
      b = data(3);
      r = data(4);
      rtype = int(data(5));
      a1 = data(6);
      a2 = data(7);
      a3 = data(8);
      a4 = data(9);

      // History variables from last converged state
      CminStrain = data(10);
      CmaxStrain = data(11);
      CshiftP = data(12);
      CshiftN = data(13);
      Cloading = int(data(14));
      CbStrain = data(15);
      CbStress = data(16);
      CrStrain = data(17);
      CrStress = data(18);
      Cplastic = data(19);

      // Copy converged history values into trial values since data is only
      // sent (received) after convergence
      TminStrain = CminStrain;
      TmaxStrain = CmaxStrain;
      TshiftP = CshiftP;
      TshiftN = CshiftN;
      Tloading = Cloading;
      TbStrain = CbStrain;
      TbStress = CbStress;
      TrStrain = CrStrain;
      TrStress = CrStress;
      Tplastic = Cplastic;

      // State variables from last converged state
      Cstrain = data(20);
      Cstress = data(21);
      Ctangent = data(22);
      CcurR = data(23);

      // Copy converged state values into trial values
      Tstrain = Cstrain;
      Tstress = Cstress;
      Ttangent = Ctangent;
      TcurR = CcurR;
   }
    
   return res;
}

void Steel03::Print (OPS_Stream& s, int flag)
{
   s << "Steel03 tag: " << this->getTag() << endln;
   s << "  type: " << rtype << " fy: " << fy << " ";
   s << "  E0: " << E0 << " ";
   s << "  b: " << b << " ";
   s << "  r:  " << r << endln;
   s << "  a1: " << a1 << " ";
   s << "  a2: " << a2 << " ";
   s << "  a3: " << a3 << " ";
   s << "  a4: " << a4 << " ";
}

