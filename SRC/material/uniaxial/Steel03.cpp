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
                                                                        
// $Revision: 1.4 $
// $Date: 2006-05-24 21:44:31 $
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
#include <elementAPI.h>

void * OPS_ADD_RUNTIME_VPV(OPS_Steel03)
{

    int argc = OPS_GetNumRemainingInputArgs() + 2;
    
    // Check that there is the minimum number of arguments
    if (argc < 9) {
        opserr << "WARNING insufficient arguments\n";
        opserr << "Want: uniaxialMaterial Steel03 tag? fy? E0? b? r? cR1 cR2?";
        opserr << " <a1? a2? a3? a4?>\n";
        return 0;
    }
      
    int tag;
    int numdata = 1;
    if (OPS_GetIntInput(&numdata, &tag) < 0) {
        opserr << "WARNING invalid uniaxialMaterial Steel03 tag\n";
        return 0;
    }
      
    // Read required Steel01 material parameters
    // fy, E, b, r, cR1, cR2;
    double data[6];
    numdata = 6;
    if (OPS_GetDoubleInput(&numdata, data) < 0) {
        opserr << "WARNING invalid double inputs\n";
        return 0;
    }
      
    // Read optional Steel01 material parameters
    // a1, a2, a3, a4
    if (argc > 9) {
	double opt[4];
	numdata = 4;

	if (argc < 13) {
	    opserr << "WARNING insufficient number of hardening parameters\n";
	    opserr << "uniaxialMaterial Steel03: " << tag << "\n";
	    return 0;
	}
	if (OPS_GetDoubleInput(&numdata, opt) < 0) {
	    opserr << "WARNING invalid double inputs\n";
	    return 0;
	}
	
        // Parsing was successful, allocate the material
	return new Steel03(tag, data[0], data[1], data[2], data[3],
			   data[4], data[5], opt[0], opt[1], opt[2], opt[3]);
    }
    else
	// Parsing was successful, allocate the material
	return new Steel03(tag, data[0], data[1], data[2], data[3],
			   data[4], data[5]);
     
}

//uniaxialMaterial Steel02 $matTag $Fy $E $b $R0 $cR1 $cR2 $a1 $a2 $a3 $a4

Steel03::Steel03
(int tag, double FY, double E, double B, double R, double R1, double R2, 
 double A1, double A2, double A3, double A4):
   UniaxialMaterial(tag,MAT_TAG_Steel03),
   fy(FY), E0(E), b(B), r(R), cR1(R1), cR2(R2), a1(A1), a2(A2), a3(A3), a4(A4)
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
 fy(0.0), E0(0.0), b(0.0), r(0.0), cR1(0.0), cR2(0.0), a1(0.0), a2(0.0), a3(0.0), a4(0.0)
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

   // Determine change in strain from last converged state
   double dStrain = strain - Cstrain;

   if (fabs(dStrain) > DBL_EPSILON) {
     // Set trial strain
     Tstrain = strain;

     // Calculate the trial state given the trial strain
     determineTrialState (dStrain);
   }

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

   // Determine change in strain from last converged state
   double dStrain = strain - Cstrain;

   if (fabs(dStrain) > DBL_EPSILON) {
     // Set trial strain
     Tstrain = strain;

     // Calculate the trial state given the trial strain
     determineTrialState (dStrain);
   }

   stress = Tstress;
   tangent = Ttangent;

   return 0;
}

double Steel03::getR (double x_in)
{
    // maybe modify this later, but it gives us better degradation at smaller strains
    x_in = fabs(x_in);
    double temp_r = r;
    
    // new input parameters are supposed to match Steel02 which look like 
    // Dhakal and Maekawa values: cr1 = 0.925, cr2 = 0.15
    // where 0.925 comes from 18.5/20.0 where R0=20 and 18.5 is the old coefficient provided
    //
    // so for old Dhakal and Maekawa R0 = 20, cr1 = 18.5/R0 = 0.925, cr2 = 0.15
    // so for old Gomes and Appleton R0 = 20, cr1 = 19.0/R0 = 0.95, cR2 = 0.3
    // so for my old model, now just use R0 = 20, cR1 = 0, cR2 = 0
    if (cR1 < 0.1 && cR2 < 0.1) {
        // Mackie, rough trilinear fit to the tangent to the x_in-r first 
        // quadrant circle.  Try using with values of R0 like 20 to 30
        temp_r = r*2.0/20.0;
        double t1 = -x_in/7+15/7*temp_r;
        double t2 = -4*x_in+6*temp_r;
        if (t1 > temp_r)
            temp_r = t1;
        if (t2 > temp_r)
            temp_r = t2;
        //opserr << "xin = " << x_in << " rout = " << temp_r << endln;
    } else {
    	temp_r = r * (1.0 - cR1*x_in/(cR2+x_in));
        if (temp_r < 0)
            temp_r = 1.0e-8;
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
   Steel03* theCopy = new Steel03(this->getTag(), fy, E0, b, r, cR1, cR2, 
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
   static Vector data(25);
   data(0) = this->getTag();

   // Material properties
   data(1) = fy;
   data(2) = E0;
   data(3) = b;
   data(4) = r;
   data(5) = cR1;
   data(6) = cR2;
   data(7) = a1;
   data(8) = a2;
   data(9) = a3;
   data(10) = a4;

   // History variables from last converged state
   data(11) = CminStrain;
   data(12) = CmaxStrain;
   data(13) = CshiftP;
   data(14) = CshiftN;
   data(15) = Cloading;
   data(16) = CbStrain;
   data(17) = CbStress;
   data(18) = CrStrain;
   data(19) = CrStress;
   data(20) = Cplastic;

   // State variables from last converged state
   data(21) = Cstrain;
   data(22) = Cstress;
   data(23) = Ctangent;
   data(24) = CcurR;

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
   static Vector data(25);
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
      cR1 = data(5);
      cR2 = data(6);
      a1 = data(7);
      a2 = data(8);
      a3 = data(9);
      a4 = data(10);

      // History variables from last converged state
      CminStrain = data(11);
      CmaxStrain = data(12);
      CshiftP = data(13);
      CshiftN = data(14);
      Cloading = int(data(15));
      CbStrain = data(16);
      CbStress = data(17);
      CrStrain = data(18);
      CrStress = data(19);
      Cplastic = data(20);

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
      Cstrain = data(21);
      Cstress = data(22);
      Ctangent = data(23);
      CcurR = data(24);

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
	if (flag == OPS_PRINT_PRINTMODEL_MATERIAL) {
		s << "Steel03 tag: " << this->getTag() << endln;
		s << " fy: " << fy << " ";
		s << "  E0: " << E0 << " ";
		s << "  b: " << b << " ";
		s << "  r:  " << r << " cR1: " << cR1 << " cR2: " << cR2 << endln;
		s << "  a1: " << a1 << " ";
		s << "  a2: " << a2 << " ";
		s << "  a3: " << a3 << " ";
		s << "  a4: " << a4 << " ";
	}

	if (flag == OPS_PRINT_PRINTMODEL_JSON) {
		s << "\t\t\t{";
		s << "\"name\": \"" << this->getTag() << "\", ";
		s << "\"type\": \"Steel03\", ";
		s << "\"E\": " << E0 << ", ";
		s << "\"fy\": " << fy << ", ";
		s << "\"b\": " << b << ", ";
		s << "\"R0\": " << r << ", ";
		s << "\"cR1\": " << cR1 << ", ";
		s << "\"cR2\": " << cR2 << ", ";
		s << "\"a1\": " << a1 << ", ";
		s << "\"a2\": " << a2 << ", ";
		s << "\"a3\": " << a3 << ", ";
		s << "\"a4\": " << a4 << "}";
	}
}

