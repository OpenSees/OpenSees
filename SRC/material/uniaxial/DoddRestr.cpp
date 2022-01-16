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
                                                                        
// $Revision: 1.20 $
// $Date: 2008-08-26 16:35:21 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/Steel01.cpp,v $
                                                                        
// Written: Dr. Ioannis Koutromanos, Virginia Tech, ikoutrom@vt.edu 
// Created: 01/14
// Revision: A
//
// Description: This file contains the C++ "wrapper" for the Dodd-Restrepo material routine 
// DoddRestr. 
//
// What: "@(#) DoddRestr.cpp, revA"


#include <DoddRestr.h>
#include <Vector.h>
#include <Matrix.h>
#include <Channel.h>
#include <Information.h>
#include <Parameter.h>

#include <string.h>

#include <math.h>
#include <float.h>


#include <elementAPI.h>
#include <OPS_Globals.h>


//#include <iostream>
//#include <fstream>
//using namespace std;

void * OPS_ADD_RUNTIME_VPV(OPS_DoddRestr)
{
  // Pointer to a uniaxial material that will be returned
  UniaxialMaterial *theMaterial = 0;

  int    iData[1];
  double dData[12];
  int numData = 1;

  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid uniaxialMaterial DoddRestr tag" << endln;
    return 0;
  }

  numData = OPS_GetNumRemainingInputArgs();

  if (numData != 9 && numData != 12) {
    opserr << "Invalid #args, want: uniaxialMaterial DoddRestr " << iData[0] << " Eo? fy? esh? esh1? fsh1? esu? fsu? Pmajor? Pminor? <slcf? tlcf? Dcrit?>>" << endln;
    return 0;
  }

  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "Invalid #args, want: uniaxialMaterial DoddRestr " << iData[0] << " Eo? fy? esh? esh1? fsh1? esu? fsu? Pmajor? Pminor? <slcf? tlcf? Dcrit?>>" << endln;
    return 0;
  }

  if (numData == 9) {
    dData[9] = 0.;
    dData[10] = 0.;
    dData[11] = 0.;
  }
  
  // Parsing was successful, allocate the material
  theMaterial = new DoddRestr(iData[0], dData[0], dData[1], dData[2], 
			      dData[3], dData[4], dData[5], dData[6], 
			      dData[7], dData[8], dData[9], dData[10], dData[11]);
  
  
  if (theMaterial == 0) {
    opserr << "WARNING could not create uniaxialMaterial of type DoddRestr Material\n";
    return 0;
  }

  return theMaterial;
}



DoddRestr::DoddRestr
(int tag, double Eo1, double fy1, double esh1b,
 double esh11, double fsh11, double esu1, double fsu1, double Pmajor1,
 double Pminor1, double slcf1, double tlcf1, double Dcrit1):
   UniaxialMaterial(tag,MAT_TAG_DoddRestr),
   Eo(Eo1), fy(fy1), esh(esh1b), esh1(esh11), fsh1(fsh11), esu(esu1), fsu(fsu1), Pmajor(Pmajor1),
   Pminor(Pminor1), slcf(slcf1), tlcf(tlcf1), Dcrit(Dcrit1)
{
   // Sets all history and state variables to initial values
   // History variables

	Ce_so=0;	
	Cf_so=0;	
	Cyield1=0;
	Cregion=0;
	Cpoint11=0;
	Cpoint12=0;
	Cpoint13=0;
	Cpoint21=0;
	Cpoint22=0;
	Cpoint23=0;
	Cpoint31=0;
	Cpoint32=0;
	Cpoint33=0;
	Cpoint41=0;
	Cpoint42=0;
	Cpoint43=0;
	Cpoint51=0;
	Cpoint52=0;
	Cpoint53=0;
	Cep_o1=0;
	Cep_o2=0;
	Cep_M=0;
	Cfps_so=Eo;
	Chist11=0;
	Chist12=0;
	Cpoint61=0;
	Cpoint62=0;
	Cpoint63=0;
	Csim1=0;
    CDam=0;



	Te_so=0;	
	Tf_so=0;	
	Tyield1=0;
	Tregion=0;
	Tpoint11=0;
	Tpoint12=0;
	Tpoint13=0;
	Tpoint21=0;
	Tpoint22=0;
	Tpoint23=0;
	Tpoint31=0;
	Tpoint32=0;
	Tpoint33=0;
	Tpoint41=0;
	Tpoint42=0;
	Tpoint43=0;
	Tpoint51=0;
	Tpoint52=0;
	Tpoint53=0;
	Tep_o1=0;
	Tep_o2=0;
	Tep_M=0;
	Tfps_so=Eo;
	Thist11=0;
	Thist12=0;
	Tpoint61=0;
	Tpoint62=0;
	Tpoint63=0;
	Tsim1=0;
    TDam=0;

   // State variables
   Cstrain = 0.0;
   Cstress = 0.0;
   Ctangent = Eo;

   Tstrain = 0.0;
   Tstress = 0.0;
   Ttangent = Eo;


}

DoddRestr::DoddRestr():UniaxialMaterial(0,MAT_TAG_DoddRestr),
   Eo(0.0), fy(0.0), esh(0.0), esh1(0.0), fsh1(0.0), esu(0.0), fsu(0.0), Pmajor(0.0),
   Pminor(0.0), slcf(0.0), tlcf(0.0), Dcrit(0.0)
{



}

DoddRestr::~DoddRestr ()
{

}


#ifdef _WIN32

#define steeldr_ STEELDR

extern "C" int  steeldr_ (double *strn1, double *strs1, double *Etangm, double *dmatr, double *hrv1);

// Add more declarations as needed


#else

extern "C" int  steeldr_ (double *strn1, double *strs1, double *Etangm, double *dmatr, double *hrv1); 


#endif


int DoddRestr::setTrialStrain (double strain, double strainRate)
{
   // Reset history variables to last converged state
	Te_so	= Ce_so;
	Tf_so	= Cf_so;
	Tyield1	= Cyield1;
	Tregion	= Cregion;
	Tpoint11= Cpoint11;
	Tpoint12= Cpoint12;
	Tpoint13= Cpoint13;
	Tpoint21= Cpoint21;
	Tpoint22= Cpoint22;
	Tpoint23= Cpoint23;
	Tpoint31= Cpoint31;
	Tpoint32= Cpoint32;
	Tpoint33= Cpoint33;
	Tpoint41= Cpoint41;
	Tpoint42= Cpoint42;
	Tpoint43= Cpoint43;
	Tpoint51= Cpoint51;
	Tpoint52= Cpoint52;
	Tpoint53= Cpoint53;
	Tep_o1	= Cep_o1;
	Tep_o2	= Cep_o2;
	Tep_M	= Cep_M;
	Tfps_so	= Cfps_so;
	Thist11	= Chist11;
	Thist12	= Chist12;
	Tpoint61= Cpoint61;
	Tpoint62= Cpoint62;
	Tpoint63= Cpoint63;
	Tsim1	= Csim1;
	TDam	= CDam;

   Tstrain = Cstrain;
   Tstress = Cstress;
   Ttangent = Ctangent;

   // Determine change in strain from last converged state
   double dStrain = strain - Cstrain;

//   if (fabs(dStrain) > DBL_EPSILON) {
     // Set trial strain
     Tstrain = strain;

     // Calculate the trial state given the trial strain
     determineTrialState (dStrain);

//   }

   return 0;
}

int DoddRestr::setTrial (double strain, double &stress, double &tangent, double strainRate)
{
   
	// Reset history variables to last converged state
	Te_so	= Ce_so;
	Tf_so	= Cf_so;
	Tyield1	= Cyield1;
	Tregion	= Cregion;
	Tpoint11= Cpoint11;
	Tpoint12= Cpoint12;
	Tpoint13= Cpoint13;
	Tpoint21= Cpoint21;
	Tpoint22= Cpoint22;
	Tpoint23= Cpoint23;
	Tpoint31= Cpoint31;
	Tpoint32= Cpoint32;
	Tpoint33= Cpoint33;
	Tpoint41= Cpoint41;
	Tpoint42= Cpoint42;
	Tpoint43= Cpoint43;
	Tpoint51= Cpoint51;
	Tpoint52= Cpoint52;
	Tpoint53= Cpoint53;
	Tep_o1	= Cep_o1;
	Tep_o2	= Cep_o2;
	Tep_M	= Cep_M;
	Tfps_so	= Cfps_so;
	Thist11	= Chist11;
	Thist12	= Chist12;
	Tpoint61= Cpoint61;
	Tpoint62= Cpoint62;
	Tpoint63= Cpoint63;
	Tsim1	= Csim1;
	TDam	= CDam;

   Tstrain = Cstrain;
   Tstress = Cstress;
   Ttangent = Ctangent;

   // Determine change in strain from last converged state
   double dStrain = strain - Cstrain;

//   if (fabs(dStrain) > DBL_EPSILON) {
     // Set trial strain
     Tstrain = strain;

     // Calculate the trial state given the trial strain
     determineTrialState (dStrain);

//   }


   stress = Tstress;
   tangent = Ttangent;

   return 0;
}

void DoddRestr::determineTrialState (double dStrain)
{


// ALL THE CODE FOR STRESS UPDATE SHOULD BE ADDED HERE!!!


	// ********************************************************************************************************************

	strn1=Cstrain+dStrain;



		// Put committed history variables in vector hrv1:
	hrv1[0]	= Ce_so;
	hrv1[1]	= Cf_so;
	hrv1[2]	= Cyield1;
	hrv1[3]	= Cregion;
	hrv1[4]= Cpoint11;
	hrv1[5]= Cpoint12;
	hrv1[6]= Cpoint13;
	hrv1[7]= Cpoint21;
	hrv1[8]= Cpoint22;
	hrv1[9]= Cpoint23;
	hrv1[10]= Cpoint31;
	hrv1[11]= Cpoint32;
	hrv1[12]= Cpoint33;
	hrv1[13]= Cpoint41;
	hrv1[14]= Cpoint42;
	hrv1[15]= Cpoint43;
	hrv1[16]= Cpoint51;
	hrv1[17]= Cpoint52;
	hrv1[18]= Cpoint53;
	hrv1[19]	= Cep_o1;
	hrv1[20]	= Cep_o2;
	hrv1[21]	= Cep_M;
	hrv1[22]	= Cfps_so;
	hrv1[23]	= Chist11;
	hrv1[24]	= Chist12;
	hrv1[25]= Cpoint61;
	hrv1[26]= Cpoint62;
	hrv1[27]= Cpoint63;
	hrv1[28]	= Csim1;
	hrv1[29]	= CDam;


	// Put significant material parameters in vector dmatr:

	dmatr[0] = Eo;
	dmatr[1] = 0.;
	dmatr[2] = fy;
	dmatr[3] = esh;
	dmatr[4] = esh1;
	dmatr[5] = fsh1;
	dmatr[6] = esu;
	dmatr[7] = fsu;
	dmatr[8] = Pmajor;
	dmatr[9] = Pminor;
	dmatr[10] = 0.;
	dmatr[11] = slcf;
	dmatr[12] = tlcf;
	dmatr[13] = Dcrit;

	// Call the external (fortran) procedure to conduct the stress update:

	steeldr_ (&strn1, &strs1, &Etangm, dmatr, hrv1);


	Tstress=strs1;
	Ttangent=Etangm;



	// Write trial values of history variables which have been returned in hrv1:

	Te_so	= hrv1[0];
	Tf_so	= hrv1[1];
	Tyield1	= hrv1[2];
	Tregion	= hrv1[3];
	Tpoint11= hrv1[4];
	Tpoint12= hrv1[5];
	Tpoint13= hrv1[6];
	Tpoint21= hrv1[7];
	Tpoint22= hrv1[8];
	Tpoint23= hrv1[9];
	Tpoint31= hrv1[10];
	Tpoint32= hrv1[11];
	Tpoint33= hrv1[12];
	Tpoint41= hrv1[13];
	Tpoint42= hrv1[14];
	Tpoint43= hrv1[15];
	Tpoint51= hrv1[16];
	Tpoint52= hrv1[17];
	Tpoint53= hrv1[18];
	Tep_o1	= hrv1[19];
	Tep_o2	= hrv1[20];
	Tep_M	= hrv1[21];
	Tfps_so	= hrv1[22];
	Thist11	= hrv1[23];
	Thist12	= hrv1[24];
	Tpoint61= hrv1[25];
	Tpoint62= hrv1[26];
	Tpoint63= hrv1[27];
	Tsim1	= hrv1[28];
	TDam	= hrv1[29];




}

//	****	END OF STRESS UPDATE ALGORITHM	*** =================================



double DoddRestr::getStrain ()
{
   return Tstrain;
}

double DoddRestr::getStress ()
{
   return Tstress;
}

double DoddRestr::getTangent ()
{
   return Ttangent;
}

int DoddRestr::commitState ()
{
   // History variables
   	Ce_so	= Te_so;
	Cf_so	= Tf_so;
	Cyield1	= Tyield1;
	Cregion	= Tregion;
	Cpoint11= Tpoint11;
	Cpoint12= Tpoint12;
	Cpoint13= Tpoint13;
	Cpoint21= Tpoint21;
	Cpoint22= Tpoint22;
	Cpoint23= Tpoint23;
	Cpoint31= Tpoint31;
	Cpoint32= Tpoint32;
	Cpoint33= Tpoint33;
	Cpoint41= Tpoint41;
	Cpoint42= Tpoint42;
	Cpoint43= Tpoint43;
	Cpoint51= Tpoint51;
	Cpoint52= Tpoint52;
	Cpoint53= Tpoint53;
	Cep_o1	= Tep_o1;
	Cep_o2	= Tep_o2;
	Cep_M	= Tep_M;
	Cfps_so	= Tfps_so;
	Chist11	= Thist11;
	Chist12	= Thist12;
	Cpoint61= Tpoint61;
	Cpoint62= Tpoint62;
	Cpoint63= Tpoint63;
	Csim1	= Tsim1;
	CDam	= TDam;


   // State variables
   Cstrain = Tstrain;
   Cstress = Tstress;
   Ctangent = Ttangent;

   return 0;
}

int DoddRestr::revertToLastCommit ()
{
   // Reset trial history variables to last committed state
   	Te_so	= Ce_so;
	Tf_so	= Cf_so;
	Tyield1	= Cyield1;
	Tregion	= Cregion;
	Tpoint11= Cpoint11;
	Tpoint12= Cpoint12;
	Tpoint13= Cpoint13;
	Tpoint21= Cpoint21;
	Tpoint22= Cpoint22;
	Tpoint23= Cpoint23;
	Tpoint31= Cpoint31;
	Tpoint32= Cpoint32;
	Tpoint33= Cpoint33;
	Tpoint41= Cpoint41;
	Tpoint42= Cpoint42;
	Tpoint43= Cpoint43;
	Tpoint51= Cpoint51;
	Tpoint52= Cpoint52;
	Tpoint53= Cpoint53;
	Tep_o1	= Cep_o1;
	Tep_o2	= Cep_o2;
	Tep_M	= Cep_M;
	Tfps_so	= Cfps_so;
	Thist11	= Chist11;
	Thist12	= Chist12;
	Tpoint61= Cpoint61;
	Tpoint62= Cpoint62;
	Tpoint63= Cpoint63;
	Tsim1	= Csim1;
	TDam	= CDam;


   // Reset trial state variables to last committed state
   Tstrain = Cstrain;
   Tstress = Cstress;
   Ttangent = Ctangent;

   return 0;
}

int DoddRestr::revertToStart ()
{
   // History variables
	Ce_so=0;	
	Cf_so=0;	
	Cyield1=0;
	Cregion=0;
	Cpoint11=0;
	Cpoint12=0;
	Cpoint13=0;
	Cpoint21=0;
	Cpoint22=0;
	Cpoint23=0;
	Cpoint31=0;
	Cpoint32=0;
	Cpoint33=0;
	Cpoint41=0;
	Cpoint42=0;
	Cpoint43=0;
	Cpoint51=0;
	Cpoint52=0;
	Cpoint53=0;
	Cep_o1=0;
	Cep_o2=0;
	Cep_M=0;
	Cfps_so=Eo;
	Chist11=0;
	Chist12=0;
	Cpoint61=0;
	Cpoint62=0;
	Cpoint63=0;
	Csim1=0;
    CDam=0;


	Te_so=0;	
	Tf_so=0;	
	Tyield1=0;
	Tregion=0;
	Tpoint11=0;
	Tpoint12=0;
	Tpoint13=0;
	Tpoint21=0;
	Tpoint22=0;
	Tpoint23=0;
	Tpoint31=0;
	Tpoint32=0;
	Tpoint33=0;
	Tpoint41=0;
	Tpoint42=0;
	Tpoint43=0;
	Tpoint51=0;
	Tpoint52=0;
	Tpoint53=0;
	Tep_o1=0;
	Tep_o2=0;
	Tep_M=0;
	Tfps_so=Eo;
	Thist11=0;
	Thist12=0;
	Tpoint61=0;
	Tpoint62=0;
	Tpoint63=0;
	Tsim1=0;
    TDam=0;

   // State variables
   Cstrain = 0.0;
   Cstress = 0.0;
   Ctangent = Eo;

   Tstrain = 0.0;
   Tstress = 0.0;
   Ttangent = Eo;



   return 0;
}



UniaxialMaterial* DoddRestr::getCopy ()
{
   DoddRestr* theCopy = new DoddRestr(this->getTag(), Eo, fy, esh, esh1,
                                  fsh1, esu, fsu, Pmajor, Pminor, slcf, tlcf, Dcrit);

   // Converged history variables
   	theCopy->Ce_so	= Ce_so;
	theCopy->Cf_so	= Cf_so;
	theCopy->Cyield1	= Cyield1;
	theCopy->Cregion	= Cregion;
	theCopy->Cpoint11= Cpoint11;
	theCopy->Cpoint12= Cpoint12;
	theCopy->Cpoint13= Cpoint13;
	theCopy->Cpoint21= Cpoint21;
	theCopy->Cpoint22= Cpoint22;
	theCopy->Cpoint23= Cpoint23;
	theCopy->Cpoint31= Cpoint31;
	theCopy->Cpoint32= Cpoint32;
	theCopy->Cpoint33= Cpoint33;
	theCopy->Cpoint41= Cpoint41;
	theCopy->Cpoint42= Cpoint42;
	theCopy->Cpoint43= Cpoint43;
	theCopy->Cpoint51= Cpoint51;
	theCopy->Cpoint52= Cpoint52;
	theCopy->Cpoint53= Cpoint53;
	theCopy->Cep_o1	= Cep_o1;
	theCopy->Cep_o2	= Cep_o2;
	theCopy->Cep_M	= Cep_M;
	theCopy->Cfps_so	= Cfps_so;
	theCopy->Chist11	= Chist11;
	theCopy->Chist12	= Chist12;
	theCopy->Cpoint61= Cpoint61;
	theCopy->Cpoint62= Cpoint62;
	theCopy->Cpoint63= Cpoint63;
	theCopy->Csim1	= Csim1;
	theCopy->CDam	= CDam;

   // Trial history variables
   	theCopy->Te_so	= Te_so;
	theCopy->Tf_so	= Tf_so;
	theCopy->Tyield1	= Tyield1;
	theCopy->Tregion	= Tregion;
	theCopy->Tpoint11= Tpoint11;
	theCopy->Tpoint12= Tpoint12;
	theCopy->Tpoint13= Tpoint13;
	theCopy->Tpoint21= Tpoint21;
	theCopy->Tpoint22= Tpoint22;
	theCopy->Tpoint23= Tpoint23;
	theCopy->Tpoint31= Tpoint31;
	theCopy->Tpoint32= Tpoint32;
	theCopy->Tpoint33= Tpoint33;
	theCopy->Tpoint41= Tpoint41;
	theCopy->Tpoint42= Tpoint42;
	theCopy->Tpoint43= Tpoint43;
	theCopy->Tpoint51= Tpoint51;
	theCopy->Tpoint52= Tpoint52;
	theCopy->Tpoint53= Tpoint53;
	theCopy->Tep_o1	= Tep_o1;
	theCopy->Tep_o2	= Tep_o2;
	theCopy->Tep_M	= Tep_M;
	theCopy->Tfps_so	= Tfps_so;
	theCopy->Thist11	= Thist11;
	theCopy->Thist12	= Thist12;
	theCopy->Tpoint61= Tpoint61;
	theCopy->Tpoint62= Tpoint62;
	theCopy->Tpoint63= Tpoint63;
	theCopy->Tsim1	= Tsim1;
	theCopy->TDam	= TDam;

   // Converged state variables
   theCopy->Cstrain = Cstrain;
   theCopy->Cstress = Cstress;
   theCopy->Ctangent = Ctangent;

   // Trial state variables
   theCopy->Tstrain = Tstrain;
   theCopy->Tstress = Tstress;
   theCopy->Ttangent = Ttangent;

   return theCopy;
}

int DoddRestr::sendSelf (int commitTag, Channel& theChannel)
{
   int res = 0;
   static Vector data(46);
   data(0) = this->getTag();

   // Material properties
   data(1) = Eo;
   data(2) = fy;
   data(3) = esh;
   data(4) = esh1;
   data(5) = fsh1;
   data(6) = esu;
   data(7) = fsu;
   data(8) = Pmajor;
   data(9) = Pminor;
   data(10) = slcf;
   data(11) = tlcf;
   data(12) = Dcrit;


   // History variables from last converged state

    data(13)	= Ce_so;
    data(14)	= Cf_so;
	data(15)	= Cyield1;
	data(16)	= Cregion;
	data(17)	= Cpoint11;
	data(18)	= Cpoint12;
	data(19)	= Cpoint13;
	data(20)	= Cpoint21;
	data(21)	= Cpoint22;
	data(22)	= Cpoint23;
	data(23)	= Cpoint31;
	data(24)	= Cpoint32;
	data(25)	= Cpoint33;
	data(26)	= Cpoint41;
	data(27)	= Cpoint42;
	data(28)	= Cpoint43;
	data(29)	= Cpoint51;
	data(30)	= Cpoint52;
	data(31)	= Cpoint53;
	data(32)	= Cep_o1;
	data(33)	= Cep_o2;
	data(34)	= Cep_M;
	data(35)	= Cfps_so;
	data(36)	= Chist11;
	data(37)	= Chist12;
	data(38)	= Cpoint61;
	data(39)	= Cpoint62;
	data(40)	= Cpoint63;
	data(41)	= Csim1;
	data(42)	= CDam;

   // State variables from last converged state
   data(43) = Cstrain;
   data(44) = Cstress;
   data(45) = Ctangent;

   // Data is only sent after convergence, so no trial variables
   // need to be sent through data vector

   res = theChannel.sendVector(this->getDbTag(), commitTag, data);
   if (res < 0) 
      opserr << "DoddRestr::sendSelf() - failed to send data\n";

   return res;
}

int DoddRestr::recvSelf (int commitTag, Channel& theChannel,
                                FEM_ObjectBroker& theBroker)
{
   int res = 0;
   static Vector data(46);
   res = theChannel.recvVector(this->getDbTag(), commitTag, data);
  
   if (res < 0) {
      opserr << "DoddRestr::recvSelf() - failed to receive data\n";
      this->setTag(0);      
   }
   else {
      this->setTag(int(data(0)));

      // Material properties
      Eo = data(1);
      fy = data(2);
      esh = data(3);
      esh1 = data(4);
      fsh1 = data(5);
	  esu = data(6);
      fsu = data(7);
      Pmajor = data(8);
	  Pminor = data(9);
	  slcf = data(10);
	  tlcf = data(11);
	  Dcrit = data(12);


      // History variables from last converged state
	  Ce_so = data(13);	
	  Cf_so = data(14);	
	  Cyield1 = data(15);
	  Cregion = data(16);
	  Cpoint11 = data(17);
	  Cpoint12 = data(18);
	  Cpoint13 = data(19);
	  Cpoint21 = data(20);
	  Cpoint22 = data(21);
	  Cpoint23 = data(22);
	  Cpoint31 = data(23);
	  Cpoint32 = data(24);
	  Cpoint33 = data(25);
	  Cpoint41 = data(26);
	  Cpoint42 = data(27);
	  Cpoint43 = data(28);
	  Cpoint51 = data(29);
	  Cpoint52 = data(30);
	  Cpoint53 = data(31);
	  Cep_o1 = data(32);
	  Cep_o2 = data(33);
	  Cep_M = data(34);
	  Cfps_so = data(35);
	  Chist11 = data(36);
	  Chist12 = data(37);
	  Cpoint61 = data(38);
	  Cpoint62 = data(39);
	  Cpoint63 = data(40);
	  Csim1 = data(41);
      CDam = data(42);


      // Copy converged history values into trial values since data is only
      // sent (received) after convergence
   	Te_so	= Ce_so;
	Tf_so	= Cf_so;
	Tyield1	= Cyield1;
	Tregion	= Cregion;
	Tpoint11= Cpoint11;
	Tpoint12= Cpoint12;
	Tpoint13= Cpoint13;
	Tpoint21= Cpoint21;
	Tpoint22= Cpoint22;
	Tpoint23= Cpoint23;
	Tpoint31= Cpoint31;
	Tpoint32= Cpoint32;
	Tpoint33= Cpoint33;
	Tpoint41= Cpoint41;
	Tpoint42= Cpoint42;
	Tpoint43= Cpoint43;
	Tpoint51= Cpoint51;
	Tpoint52= Cpoint52;
	Tpoint53= Cpoint53;
	Tep_o1	= Cep_o1;
	Tep_o2	= Cep_o2;
	Tep_M	= Cep_M;
	Tfps_so	= Cfps_so;
	Thist11	= Chist11;
	Thist12	= Chist12;
	Tpoint61= Cpoint61;
	Tpoint62= Cpoint62;
	Tpoint63= Cpoint63;
	Tsim1	= Csim1;
	TDam	= CDam;


      // State variables from last converged state
      Cstrain = data(43);
      Cstress = data(44);
      Ctangent = data(45);      

      // Copy converged state values into trial values
      Tstrain = Cstrain;
      Tstress = Cstress;
      Ttangent = Ctangent;
   }
    
   return res;
}

void DoddRestr::Print (OPS_Stream& s, int flag)
{
   s << "Concshcr tag: " << this->getTag() << endln;
   s << "  Eo: " << Eo << " ";
   s << "  fy: " << fy << " ";
   s << "  esh:  " << esh << " ";
   s << "  esh1: " << esh1 << " ";
   s << "  fsh1: " << fsh1 << " ";
   s << "  esu: " << esu << " ";
   s << "  fsu: " << fsu << " ";
   s << "  Pmajor: " << Pmajor << " ";
   s << "  Pminor: " << Pminor << " ";
   s << "  slcf: " << slcf << " ";
   s << "  tlcf: " << tlcf << " ";
   s << "  Dcrit: " << Dcrit << " ";
}




