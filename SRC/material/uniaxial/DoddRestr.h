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
                                                                        
// $Revision: 1.6 $
// $Date: 2014-01-16 23:42:19 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/DoddRestr.h,v $
                                                             
#ifndef DoddRestr_h
#define DoddRestr_h

// Written: Dr. Ioannis Koutromanos, Virginia Tech, ikoutrom@vt.edu
// Created: 02/15
// Revision: A
//
// Description: This file contains the class definition for 
// DoddRestr. DoddRestr provides the abstraction
// of a Dodd-Restrepo material which can also account for 
// fracture due to low-cycle fatigue 
//
// What: "@(#) DoddRestr.h, revA"

#include <UniaxialMaterial.h>

class DoddRestr : public UniaxialMaterial
{
  public:
    DoddRestr(int tag, double Eo, double fy, double esh, double esh1,  double fsh1,  double esu, double fsu, double Pmajor, double Pminor);    
    DoddRestr(int tag, double Eo, double fy, double esh, double esh1,  double fsh1,  double esu, double fsu, double Pmajor, double Pminor, double slcf, double tlcf, double Dcrit);    
    DoddRestr();    

    ~DoddRestr();

    const char *getClassType(void) const {return "DoddRestr";};

    int setTrialStrain(double strain, double strainRate = 0.0);
	int setTrial (double strain, double &stress, double &tangent, double strainRate = 0.0);
    double getStrain(void);          
    double getStress(void);
    double getTangent(void);

    double getInitialTangent(void) {return Eo;};

    int commitState(void);
    int revertToLastCommit(void);    
    int revertToStart(void);    

    UniaxialMaterial *getCopy(void);
    
    int sendSelf(int commitTag, Channel &theChannel);  
    int recvSelf(int commitTag, Channel &theChannel, 
                 FEM_ObjectBroker &theBroker);    
    
    void Print(OPS_Stream &s, int flag =0);
    
  protected:
    
  private:
    double Eo;		 // initial Young's modulus
    double fy;       // yield stress
    double esh;           // strain at end of yield plateau
    double esh1;			 // strain at point of strain hardening region 
	double fsh1;				 // stress at point of strain hardening region
    double esu;				 // ultimate strain
	double fsu;				 //	ultimate tensile strength
	double Pmajor;			 // coefficient P for shape of major unloading-reloading curves (if negative, it is calculated automatically)
	double Pminor;		     // coefficient P for shape of minor unloading-reloading curves
	double slcf;		     // parameter S of low-cycle fatigue criterion (optional)
	double tlcf;			 // parameter t of low-cycle fatigue criterion (optional)
	double Dcrit;			 // critical damage parameter for fracture induced by low-cycle fatigue
	
	
	double hrv1[30];		// temporary vector with history variables which is passed to and from the fortran routine with the DoddRestr model
	double dmatr[14];		// temporary vector material parameters passed to the fortran routine with the DoddRestr model
	double strn1;
	double strs1;
	double Etangm;
	//  Trial Response quantities


    double Tstrain;      // current trial strain
    double Tstress;      // current trial stress
    double Ttangent;     // current trial tangent

//	Committed Response quantities

    double Cstrain;     // last committed strain
    double Cstress;     // last committed stress
    double Ctangent;    // last committed  tangent


//	Committed History Variables

    double Cepso;     // last committed shrinkage strain
    double Ch1;     // last committed h1
    double Ch2;    // last committed h2
	double Ctn;	   // last committed "time for shrinkage calculation"!		





	double Ce_so;	//last committed strain
	double Cf_so;	//last committed stress
	double Cyield1;
	double Cregion;
	double Cpoint11;
	double Cpoint12;
	double Cpoint13;
	double Cpoint21;
	double Cpoint22;
	double Cpoint23;
	double Cpoint31;
	double Cpoint32;
	double Cpoint33;
	double Cpoint41;
	double Cpoint42;
	double Cpoint43;
	double Cpoint51;
	double Cpoint52;
	double Cpoint53;
	double Cep_o1;
	double Cep_o2;
	double Cep_M;
	double Cfps_so;
	double Chist11;
	double Chist12;
	double Cpoint61;
	double Cpoint62;
	double Cpoint63;
	double Csim1;
    double CDam;





//	Trial History Variables

	double Te_so;	
	double Tf_so;	
	double Tyield1;
	double Tregion;
	double Tpoint11;
	double Tpoint12;
	double Tpoint13;
	double Tpoint21;
	double Tpoint22;
	double Tpoint23;
	double Tpoint31;
	double Tpoint32;
	double Tpoint33;
	double Tpoint41;
	double Tpoint42;
	double Tpoint43;
	double Tpoint51;
	double Tpoint52;
	double Tpoint53;
	double Tep_o1;
	double Tep_o2;
	double Tep_M;
	double Tfps_so;
	double Thist11;
	double Thist12;
	double Tpoint61;
	double Tpoint62;
	double Tpoint63;
	double Tsim1;
    double TDam;

	void determineTrialState (double dStrain);


};


#endif



 
