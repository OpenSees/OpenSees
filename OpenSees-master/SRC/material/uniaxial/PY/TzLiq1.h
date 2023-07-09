/* *********************************************************************
**    Module:	TzLiq1.h 
**
**    Purpose:	Provide a t-z material that gets pore pressure from a
**				specified element that contains a PorousFluidSolid.
**              
**
**    Developed by Ross W. Boulanger
**    (C) Copyright 2002, All Rights Reserved.
**
** ****************************************************************** */

// $Revision: 1.0
// $Date: 2002/2/5
// $Source: /OpenSees/SRC/material/uniaxial/TzLiq1.h

#ifndef TZLIQ1_H
#define TZLIQ1_H

// Written: RWB
// Created: Feb 2002
//
// Description: This file contains the class definition for TzLiq1.
// 

#include "TzSimple1.h"
#include <UniaxialMaterial.h>
#include <Domain.h>

#include <string>

class FourNodeQuad;
class TimeSeries;

//#include <SSPquadUP.h>
//#include <SSPquad.h>

class TzLiq1 : public TzSimple1 {
  public:
    TzLiq1(int tag, int classtag, int tzType, double tult, double z50,
		      double dashpot, int solidElem1, int solidElem2, Domain *theDomain);
	TzLiq1(int tag, int classtag, int tzType, double tult, double z50,
			  double dashpot, Domain *theDomain, TimeSeries *theSeries);
    TzLiq1();
    ~TzLiq1();

    const char *getClassType(void) const {return "TzLiq1";};

    int setTrialStrain(double y, double yRate); 
    double getStrain(void);          
    double getStress(void);
    double getTangent(void);
    double getStrainRate(void);
    double getDampTangent(void);
    double getInitialTangent(void);

    int commitState(void);
    int revertToLastCommit(void);    
    int revertToStart(void);        

    UniaxialMaterial *getCopy(void);
    
    int sendSelf(int commitTag, Channel &theChannel);  
    int recvSelf(int commitTag, Channel &theChannel, 
		 FEM_ObjectBroker &theBroker);

    //  Command for initiating vertConsolStress from TclUpdateMaterialStageCommand
     int setParameter(const char **argv, int argc, Parameter &param);
     int updateParameter(int snum, Information &eleInformation);
    
    void Print(OPS_Stream &s, int flag =0);

   
  protected:
    
  private:

	// Committed and trial values for t, z, and ru
	double Tz;
	double Cz;
	double Tt;
	double Ct;
	double Tangent;
	double maxTangent;
	double Tru;
	double Cru;
	double Hru;

	// Solid element from which pore pressures are obtained, domain pointer
	// and stage information to get the initial vertical effective stress.
	int solidElem1;
	int solidElem2;
	int theSeriesTag;

	double meanConsolStress;
	double ru;

    static int loadStage;
	int    lastLoadStage;
	std::string elemFlag; 
	Domain *theDomain;
	TimeSeries *theSeries;
	FourNodeQuad *theQuad1;
	FourNodeQuad *theQuad2;

	// Initial tangent
	double initialTangent;
	
	// Function for obtaining effective stresses from adjoining solid soil elements
	double getEffectiveStress(void);
	double getEffectiveStress(TimeSeries *theSeries);
	static Vector stressV3;
	
};

#endif // TZLIQ1_H



