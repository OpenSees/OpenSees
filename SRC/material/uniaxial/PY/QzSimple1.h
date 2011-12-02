/* *********************************************************************
**    Module:	QzSimple1.h 
**
**    Purpose:	Provide a simple Q-z material for OpenSees
**
**
**    Developed by Ross W. Boulanger
**    (C) Copyright 2002, All Rights Reserved.
**
** ****************************************************************** */

// $Revision: 1.0
// $Date: 2002/1/22
// $Source: /OpenSees/SRC/material/uniaxial/QzSimple1.h

#ifndef QZSIMPLE1_H
#define QZSIMPLE1_H

// Written: RWB
// Created: Jan 2002
// tested and checked: Boris Jeremic (jeremic@ucdavis.edu) Spring 2002
//
// Description: This file contains the class definition for QzSimple1.
// 

#include <UniaxialMaterial.h>


	// Controls on internal iteration between components
	static int QZmaxIterations = 20;      
	static double QZtolerance  = 1.0e-12; 


class QzSimple1 : public UniaxialMaterial
{
  public:
    QzSimple1(int tag, int qzType, double Qult, double z50, double suction,
		      double dashpot);
    QzSimple1();
    ~QzSimple1();

    int setTrialStrain(double z, double zRate); 
    double getStrain(void);          
    double getStress(void);
    double getTangent(void);
    double getInitialTangent(void);
    double getStrainRate(void);
    double getDampTangent(void);

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

	// Functions to get Q & z for each component individually
	void getGap(double zlast, double dz, double dz_old);
	void getClosure(double zlast, double dz);
	void getSuction(double zlast, double zy);
	void getNearField(double zlast, double dz, double dz_old);
	void getFarField(double z);

    // Material parameters
	int    QzType;		// Q-z relation selection
    double Qult;		// Material capacity
    double z50;			// z at 50% of Qult in compression
    double suction;		// ratio of max suction force to Qult
	double zref;		// reference point for Near Field component
	double np;			// exponent for hardening shape of Near Field component
	double Elast;		// Q/Qult when yielding first occurs in virgin compression
	double maxElast;	// max size of elastic range (in terms of dQ/Qult)
	double nd;			// exponent for hardening shape of suction component
	double dashpot;     // dashpot on the far-field (elastic) component

	// Generated parameters or constants (not user input)
	double NFkrig;		// stiffness of the "rigid" portion of Near Field 
	
    // Committed history variables for entire Q-z material
    double Cz;			// Committed z
    double CQ;			// Committed Q
    double Ctangent;	// Committed tangent

	// Trial history variables for entire Q-z material
    double Tz;			// Trial z
    double TQ;			// Trial Q
    double Ttangent;	// Trial tangent
	double TzRate;      // Trial velocity

	// Committed internal parameters for the NearField rigid-plastic component
	double CNF_Qinr;		//  Q at start of current plastic loading cycle - right
	double CNF_Qinl;		//  Q at start of current plastic loading cycle - left
	double CNF_zinr;		//  z at start of current plastic loading cycle - right
	double CNF_zinl;		//  z at start of current plastic loading cycle - left
	double CNF_Q;			//  current Q
	double CNF_z;			//  current z
	double CNF_tang;		//  tangent

	// Trial internal parameters for the NearField plastic component
	double TNF_Qinr;		//  Q at start of current plastic loading cycle - right
	double TNF_Qinl;		//  Q at start of current plastic loading cycle - left
	double TNF_zinr;		//  z at start of current plastic loading cycle - right
	double TNF_zinl;		//  z at start of current plastic loading cycle - left
	double TNF_Q;			//  current Q
	double TNF_z;			//  current z
	double TNF_tang;		//  tangent

	// Committed internal parameters for the Suction component
	double CSuction_Qin;	//  Q at start of current plastic loading cycle
	double CSuction_zin;	//  z at start of current plastic loading cycle
	double CSuction_Q;		//  current Q
	double CSuction_z;		//  current z
	double CSuction_tang;	//  tangent

	// Trial internal parameters for the Suction component
	double TSuction_Qin;	//  Q at start of current plastic loading cycle
	double TSuction_zin;	//  z at start of current plastic loading cycle
	double TSuction_Q;		//  current Q
	double TSuction_z;		//  current z
	double TSuction_tang;	//  tangent

	// Committed internal parameters for the Closure component
	double CClose_Q;		//  current Q
	double CClose_z;		//  current z
	double CClose_tang;		//  tangent

	// Trial internal parameters for the Closure component
	double TClose_Q;		//  current Q
	double TClose_z;		//  current z
	double TClose_tang;		//  tangent

	// Committed internal parameters for the Gap (Suction + Closure)
	double CGap_z;			//	z
	double CGap_Q;			//  combined Q
	double CGap_tang;		//  combined tangent

	// Trial internal parameters for the Gap (Suction + Closure)
	double TGap_z;			//	z
	double TGap_Q;			//  combined Q
	double TGap_tang;		//  combined tangent

	// Committed internal parameters for the Far Field component
	double CFar_z;			//  z
	double CFar_Q;			//  current Q
	double CFar_tang;       //  tangent

	// Trial internal parameters for the Far Field component
	double TFar_z;			//  z
	double TFar_Q;			//  current Q
	double TFar_tang;       //  tangent

	double initialTangent;
};

#endif
