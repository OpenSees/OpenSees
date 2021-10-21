/* *********************************************************************
**    Module:	QzLiq1.h 
**
**    Purpose:	Provide a q-z material that gets pore pressure from a
**				specified element that contains a PorousFluidSolid.
**              
**
**    Developed by Sumeet Kumar Sinha
**    (C) Copyright 2002, All Rights Reserved.
**
** ****************************************************************** */

// $Revision: 1.0
// $Date: 2020/6/12
// $Source: /OpenSees/SRC/material/uniaxial/QzLiq1.h

#ifndef QZLIQ1_H
#define QZLIQ1_H

// Written: RWB
// Created: Feb 2002
//
// Description: This file contains the class definition for QzLiq1.
// 

#include <UniaxialMaterial.h>
#include <Domain.h>
#include <FourNodeQuad.h>
#include <FluidSolidPorousMaterial.h>
#include <QzSimple1.h>
#include <iostream>
#include <FourNodeQuadUP.h>
#include <Nine_Four_Node_QuadUP.h>
#include <TimeSeries.h>


class QzLiq1 : public QzSimple1
{
  public:
    QzLiq1(int tag, int qzType, double tult, double z50, double suction,
		      double dashpot, double alpha, int solidElem1, int solidElem2, Domain *theDomain);
	QzLiq1(int tag, int qzType, double tult, double z50, double suction, 
			  double dashpot, double alpha, Domain *theDomain, TimeSeries *theSeries);
    QzLiq1();
    ~QzLiq1();

    const char *getClassType(void) const {return "QzLiq1";};

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

	// Committed and trial values for q, z, and ru
	double Tz; // Trial z (displacement)
	double Cz; // Commit z (displacement)
	double Tt; // Trial load
	double Ct; // Commit load
	double Tangent; // Tangent
	double maxTangent; // maximum tangent
	double Tru; // Trial ru
	double Cru; // Commit ru
	double Hru; // 
	double alpha; // factor (1-ru)^alpha

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

#endif // QZLIQ1_H



