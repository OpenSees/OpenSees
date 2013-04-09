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
// $Date: 2012/05/01 01:00:00 $
// $Source: /usr/local/cvs/OpenSees/PACKAGES/NewElement/cpp/RotationShearCurve.h,v $
                                                                        
#ifndef RotationShearCurve_h
#define RotationShearCurve_h

// Written: MRL 
//
// Description: This file contains the implementation for the rotation based shear curve.
//
// What: "RotationShearCurve.h Rev A"

#include <LimitCurve.h>

class Element;
class Domain;
class DomainComponent;
class Node;

#define LIMCRV_TAG_RotShearCurve 998

class RotationShearCurve : public LimitCurve
{
  public:
	RotationShearCurve(int crvTag, int eleTag,
	int ndI, int ndj, int rotAxis, double Vn, double Vr, double Kdeg, double rotLim, int defTyp,
	double b, double d, double h, double L, double st, double As, double Acc, double ld, double db, double rhot, double fc,
	double fy, double fyt, double delta, Domain *theDom, Element *theEle, Node *theNdI, Node *theNdJ);

	RotationShearCurve();	
    ~RotationShearCurve();	

    LimitCurve *getCopy (void);

	int checkElementState(double springForce);

	double getDegSlope(void);
	double getResForce(void);
	double getUnbalanceForce(void);

	int sendSelf(int commitTag, Channel &theChannel);  
	int recvSelf(int commitTag, Channel &theChannel, 
		FEM_ObjectBroker &theBroker);    
    
	void Print(OPS_Stream &s, int flag =0);

	double findLimit(double V);
	double findCritLimit(double V, double M);

	int revertToStart(void); 
	void setDegSlope(double V);
	void getElemForces(void);		

  protected:
    
  private:
	// Associated beam-colum element information
	int curveTag;		// tag for shear curve
	int eleTag;			// tag for associated beam column element

	Element *theElement;// element pointer
	Domain *theDomain;	// domain pointer
	Node *theNodeI;		// I node pointer
	Node *theNodeJ;		// J node pointer
	Vector *forceVec;	//vector of basic forces from beam column

	double thetaMin;	// lowest rotation limit where shear failure is allowed
	double P;			// axial load in associated beam-column element (kips)
	double M;			// moment in associated beam-column element (kips)
	int stateFlag;		// state of limitstate material
						// stateFlag = 0: prior to failure
						// stateFlag = 1: first time limit surface is reached
						// stateFlag = 2: on limit surface after first failure
						// stateFlag = 3: off limit surface after first failure

	int ndI;			// node I for determining rotation (column base)
	int ndJ;			// node J for determining rotation (should be located distance h from column base)
	int rotAxis;		// direction to indicate axis of rotation to be monifored
	double Vn;			// nominal shear critical capacity (kips)
	double Vr;			// residual capacity of column after failure (kips)
	double Kdeg;		// degrading slope for LimitState material after failure (kips/in)
	double rotLim;		// user defined rotation limit if calibrated SF model is not used (radians)
	int defType;		// 0 Forces user input rotatation limit
						// 1 Flexure-Shear capacity based on theta_total rotation capacity
						// 2 Flexure-Shear capacity based on theta_flexural rotation capacity
						// 3 Flexure-Shear capacity based on theta_total-plastic rotation capacity
						// 4 Flexure-Shear capacity based on theta_flexural-plastic rotation capacity
	double b;			// column width (in)
	double d;			// column effective depth (in)
	double h;			// column cross section height (in)
	double L;			// column clear span length (in)
	double st;			// transverse reinforcement spacing (in)
	double As;			// area of longitudinal steel bars in column section (in^2)
	double Acc;			// gross confined concrete area bounded by the transverse reinforcement in column section (in^2)
	double ld;			// development length of longitudinal bars using ACI 318-08 Equations 12-1 & 12-2 (in)
	double db;			// diameter of longitudinial bars (in)
	double rhot;		// transverse reinforcement ratio Ast/(st*b)
	double fc;			// concrete strength (ksi)
	double fy;			// longitudinal steel yield strength (ksi)
	double fyt;			// transverse steel yield strength (ksi)
	double delta;		// rotational offset to modify the calibrated shear failure model (radians)
};
#endif

