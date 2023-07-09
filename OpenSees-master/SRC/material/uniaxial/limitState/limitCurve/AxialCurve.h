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
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/limitState/limitCurve/AxialCurve.h,v $

// Written: kje
// Created: 08/01
// Modified: 07/02
// Revision: A
//
// Description: This file contains the class definition for 
// AxialCurve. Defines the curve used by LimitStateMaterial  
// to determine if a limit state has been reached.


#ifndef AxialCurve_h
#define AxialCurve_h

#define TAG_AxialCurve	1973 //???????? 

#include <LimitCurve.h>
#include <tcl.h>

class Element;
class Domain;

class AxialCurve : public LimitCurve
{
  public:
	AxialCurve (Tcl_Interp *theTclInterp, int tag, int eleTag, Domain *theDomain, 
		    double Fsw, double Kdeg, double Fres, //SDK
		    int defType, int forType, 
		    int ndI = 0, int ndJ = 0, int dof = 0, int perpDirn = 0,
		    double delta = 0.0,	int eleRemove = 0);
	AxialCurve();
	~AxialCurve();

    LimitCurve *getCopy (void);

	int checkElementState(double springForce);

	double getDegSlope(void);
	double getResForce(void);
	double getUnbalanceForce(void);	// get change in axial load for next time step

	int sendSelf(int commitTag, Channel &theChannel);  
	int recvSelf(int commitTag, Channel &theChannel, 
		FEM_ObjectBroker &theBroker);    
    
	void Print(OPS_Stream &s, int flag =0);

	double findLimit(double input);

	int revertToStart(void);  
	

	int    setParameter (const char **argv, int argc, Parameter &param);
	int    updateParameter          (int parameterID, Information &info);



  protected:
    
  private:

	Tcl_Interp *theTclInterp;

	// Associated beam-colum element information
	int eleTag;			// tag for associated beam column element

	Domain *theDomain;	// needed to retrieve element pointer
	Element *theElement;// element pointer
	

	int stateFlag;	// state of limitstate material
					// stateFlag = 0: prior to failure
					// stateFlag = 1: first time limit surface is reached
					// stateFlag = 2: on limit surface after first failure
					// stateFlag = 3: off limit surface after first failure
					// stateFlag = 4: at residual axial capacity

	double dP;		   // change in axial load
	double dP_old;     // History variable SDK
	double deform_old; // History variable SDK
	double failDrift;   // drift at failure SDK

	double Fsw;		// Asw*fy*dcore/s
	




	double Kdeg;	// degrading slope for LimitState material after failure
	double Fres;	// residual capacity for LimitState material after failure (assumed to be positive)
	int defType, forType; // flags indicating the axes of the limit state surface
						  //	defType = 1 for max chord rotation
						  //		    = 2 for drift based on displ of nodes ndI and ndJ
	                      //    forType = 0 for force in associated LimitState material
	                      //            = 1 for shear from beam-column element
						  //            = 2 for axial load from beam-column element

	int ndI;		// nodes for determining interstory drift
	int ndJ;
	int dof;		// degree of freedom in which drift is desired
	int perpDirn;	// perpendicular direction to dof to get distance between nodes
	double oneOverL;// 1/dist between nodes

	int eleRemove;  // option to remove element when failure is detected (not fully implemented)
					//   eleRemove = 0 (default) do not remove element
					//   eleRemove = 1 remove element after failure surface is exceeded
					//   eleRemove = 2 value after element has been removed
	double delta;   // drift ratio used to shift failure surface


	double theta2; //SDK
	double sigma;  //SDK
	double eps_normal; //SDK

	int stepCounter; //Terje





};


#endif

