/* *********************************************************************
**    Module:	PyLiq1.h 
**
**    Purpose:	Provide a p-y material that gets pore pressure from a
**				specified element that contains a PorousFluidSolid.
**              
**
**    Developed by Ross W. Boulanger
**    (C) Copyright 2002, All Rights Reserved.
**
** ****************************************************************** */

// $Revision: 1.0
// $Date: 2002/5/15
// $Source: /OpenSees/SRC/material/uniaxial/PyLiq1.h

#ifndef PYLIQ1_H
#define PYLIQ1_H

// Written: RWB
// Created: May 2002
//
// Description: This file contains the class definition for PyLiq1.
// 

#include <UniaxialMaterial.h>
#include <Domain.h>
#include <FourNodeQuad.h>
#include <FluidSolidPorousMaterial.h>
#include <PySimple1.h>
#include <iostream>
#include <string>

class PyLiq1 : public PySimple1
{
  public:
    PyLiq1(int tag, int classtag, int soilType, double pult, double y50, double drag,
	   double dashpot, double pRes, int solidElem1, int solidElem2, Domain *theDomain);
    PyLiq1();
    ~PyLiq1();

    const char *getClassType(void) const {return "PyLiq1";};

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
    int updateParameter(int snum, Information &eleInformation);
    
    void Print(OPS_Stream &s, int flag =0);
    
  protected:
    
  private:

	// Residual p (other parameters in PySimple1 base class)
	double pRes;

	// Committed and trial values for p, y, and ru
	double Ty;
	double Cy;
	double Tp;
	double Cp;
	double Tangent;
	double maxTangent;
	double Tru;
	double Cru;
	double Hru;

	// Solid element from which pore pressures are obtained, domain pointer
	// and stage information to get the initial vertical effective stress.
	int solidElem1;
	int solidElem2;
	double meanConsolStress;
    static int loadStage;
	int lastLoadStage;
	std::string elemFlag; 
	Domain *theDomain;
	FourNodeQuad *theQuad1;
	FourNodeQuad *theQuad2;

	// Initial tangent
	double initialTangent;
	
	// Function for obtaining effective stresses from adjoining solid soil elements
	double getEffectiveStress(void);
	static Vector stressV3;
	
};

#endif // PYLIQ1_H
