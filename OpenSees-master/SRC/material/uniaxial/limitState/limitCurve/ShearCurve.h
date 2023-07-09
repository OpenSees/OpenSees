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
// $Date: 2006-09-05 22:39:58 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/limitState/limitCurve/ShearCurve.h,v $
                                                                        
// Written: kje
// Created: 08/01
// Modified: 09/02
// Revision: A
//
// Description: This file contains the class definition for 
// ShearCurve. Defines the curve used by LimitStateMaterial  
// to determine if shear failure has occurred according to empirical 
// drift capacity model by Elwood (2002).


#ifndef ShearCurve_h
#define ShearCurve_h

#define TAG_ShearCurve	1972 //???????? 

#include <LimitCurve.h>

class Element;
class Domain;

class ShearCurve : public LimitCurve
{
  public:
	ShearCurve(int tag, int eleTag, Domain *theDomain, 
		double rho, double fc, double b, double h, double d, double Fsw, //SDK
		double Kdeg, double Fres, 
		int defType, int forType,
		int ndI = 0, int ndJ = 0, int dof = 0, int perpDirn = 0,
		double delta = 0.0);
	ShearCurve();
    ~ShearCurve();

    LimitCurve *getCopy (void);

	int checkElementState(double springForce);

	double getDegSlope(void);
	double getResForce(void);
	double getUnbalanceForce(void);

	int sendSelf(int commitTag, Channel &theChannel);  
	int recvSelf(int commitTag, Channel &theChannel, 
		FEM_ObjectBroker &theBroker);    
    
	void Print(OPS_Stream &s, int flag =0);

	double findLimit(double input);

	int revertToStart(void);        
	int setParameter(const char **argv, int argc, Parameter &param);
	int    updateParameter          (int parameterID, Information &info);


  protected:
    
  private:
    void setDegSlope(double V, double Dshear); // sets degrading slope upon shear failure
											   // based on calc drift at axial failure
											   // if Kdeg >= 0.0
    
	// Associated beam-colum element information
	int eleTag;			// tag for associated beam column element

	Domain *theDomain;	// needed to retrieve element pointer
	Element *theElement;// element pointer


	int stateFlag;	// state of limitstate material
					// stateFlag = 0: prior to failure
					// stateFlag = 1: first time limit surface is reached
					// stateFlag = 2: on limit surface after first failure
					// stateFlag = 3: off limit surface after first failure

	double Kdeg;	// degrading slope for LimitState material after failure
	double Fres;	// residual capacity for LimitState material after failure (assumed to be positive)
	int defType, forType; // flags indicating the axes of the limit state surface
						  //	defType = 1 for max chord rotation
						  //		    = 2 for drift based on displ of nodes ndI and ndJ
	                      //    forType = 0 for force in associated LimitState material
	                      //			= 1 for shear from beam-column element
						  //            = 2 for axial load from beam-column element
	
	double rho;     // transverse reinforcement ratio Ast/b/s
	double fc;		// concrete strength in psi
	double b;		// column width
	double h;       // column cross section height
	double d;		// column effective depth

	int ndI;		// nodes for determining interstory drift
	int ndJ;
	int dof;		// degree of freedom in which drift is desired
	int perpDirn;	// perpendicular direction to dof to get distance between nodes
	double oneOverL;// 1/dist between nodes

	double P;		// axial load in associated beam-column element
	double Fsw;		// Asw*fy*dcore/s



	double delta;   // drift ratio used to shift failure surface

	double theta1;
	double theta4;
	double theta5;
	double sigma;
	double eps_normal;
	
};


#endif
