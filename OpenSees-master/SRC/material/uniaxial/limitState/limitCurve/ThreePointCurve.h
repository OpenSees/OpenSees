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
                                                                        
// $Revision: 1.1 $
// $Date: 2006-02-07 23:15:55 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/limitState/limitCurve/ThreePointCurve.h,v $


// Written: kje
// Created: 08/01
// Modified: 07/02
// Revision: A
//
// Description: This file contains the class definition for 
// ThreePointCurve. Defines the curve used by HystereticMaterial  
// to determine if a limit state has been reached.


#ifndef ThreePointCurve_h
#define ThreePointCurve_h

#define TAG_ThreePointCurve	1972 //???????? 

#include <LimitCurve.h>

class Element;
class Domain;

class ThreePointCurve : public LimitCurve
{
  public:
	ThreePointCurve(int tag, int eleTag, Domain *theDomain, 
		double x1, double y1,
		double x2, double y2,
		double x3, double y3,
		double Kdeg, double Fres,
		int defType, int forType,
		int ndI = 0, int ndJ = 0, int dof = 0, int perpDirn = 0);
	ThreePointCurve();
    ~ThreePointCurve();

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

  protected:
    
  private:
	// Associated beam-colum element information
	int eleTag;			// tag for associated beam column element
	Element *theElement;// element pointer
	Domain *theDomain;	// needed to retrieve element pointer

	int stateFlag;	// state of hysteretic material
					// stateFlag = 0: prior to failure
					// stateFlag = 1: first time limit surface is reached
					// stateFlag = 2: on limit surface after first failure
					// stateFlag = 3: off limit surface after first failure

	double Kdeg;	// degrading slope for hysteretic material after failure
	double Fres;	// residual capacity for hysteretic material after failure (assumed to be positive)
	int defType, forType; // flags indicating the axes of the limit state surface
						  //	defType = 1 for max chord rotation
						  //		    = 2 for drift based on displ of nodes ndI and ndJ
	                      //    forType = 0 for force in associated hysteretic material
						  //            = 1 for shear from beam-column element
						  //            = 2 for axial load from beam-column element
	double x1;		// corner points of limit surface
	double x2;
	double x3;
	double y1;
	double y2;
	double y3;

	int ndI;		// nodes for determining interstory drift
	int ndJ;
	int dof;		// degree of freedom in which drift is desired
	int perpDirn;	// perpendicular direction to dof to get distance between nodes
	double oneOverL;// 1/dist between nodes

	int count; //for debugging

};


#endif
