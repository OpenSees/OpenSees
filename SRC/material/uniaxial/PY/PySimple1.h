/* *********************************************************************
**    Module:	PySimple1.h 
**
**    Purpose:	Provide a simple p-y spring for OpenSees
**              
**
**    Developed by Ross W. Boulanger
**
** Copyright @ 2002 The Regents of the University of California (The Regents). All Rights Reserved.
**
** The Regents grants permission, without fee and without a written license agreement, for (a) use, 
** reproduction, modification, and distribution of this software and its documentation by educational, 
** research, and non-profit entities for noncommercial purposes only; and (b) use, reproduction and 
** modification of this software by other entities for internal purposes only. The above copyright 
** notice, this paragraph and the following three paragraphs must appear in all copies and modifications 
** of the software and/or documentation.
**
** Permission to incorporate this software into products for commercial distribution may be obtained 
** by contacting the University of California 
** Office of Technology Licensing 
** 2150 Shattuck Avenue #510, 
** Berkeley, CA 94720-1620, 
** (510) 643-7201.
**
** This software program and documentation are copyrighted by The Regents of the University of California. 
** The Regents does not warrant that the operation of the program will be uninterrupted or error-free. The 
** end-user understands that the program was developed for research purposes and is advised not to rely 
** exclusively on the program for any reason.
**
** IN NO EVENT SHALL REGENTS BE LIABLE TO ANY PARTY FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR 
** CONSEQUENTIAL DAMAGES, INCLUDING LOST PROFITS, ARISING OUT OF THE USE OF THIS SOFTWARE AND ITS 
** DOCUMENTATION, EVEN IF REGENTS HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGE. REGENTS GRANTS 
** NO EXPRESS OR IMPLIED LICENSE IN ANY PATENT RIGHTS OF REGENTS BUT HAS IMPLEMENTED AN INDIVIDUAL 
** CONTRIBUTOR LICENSE AGREEMENT FOR THE OPENSEES PROJECT AT THE UNIVERISTY OF CALIFORNIA, BERKELEY 
** TO BENEFIT THE END USER.
**
** REGENTS SPECIFICALLY DISCLAIMS ANY WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES 
** OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE. THE SOFTWARE AND ACCOMPANYING DOCUMENTATION,
** IF ANY, PROVIDED HEREUNDER IS PROVIDED "AS IS". REGENTS HAS NO OBLIGATION TO PROVIDE MAINTENANCE, 
** SUPPORT, UPDATES, ENHANCEMENTS, OR MODIFICATIONS.
**
** ****************************************************************** */

// $Revision: 1.0
// $Date: 2001/10/15
// $Source: /OpenSees/SRC/material/uniaxial/PySimple1.h

#ifndef PYSIMPLE1_H
#define PYSIMPLE1_H

// Written: RWB
// Created: Oct 2001
// tested and checked: Boris Jeremic (jeremic@ucdavis.edu) Spring 2002
//
// Description: This file contains the class definition for PySimple1.
// 

#include <UniaxialMaterial.h>

class PySimple1 : public UniaxialMaterial
{
  public:
    PySimple1(int tag, int classtag, int soilType, double pult, double y50, 
	      double drag, double dashpot);
    PySimple1();
    ~PySimple1();


    const char *getClassType(void) const {return "PySimple1";};

    int setTrialStrain(double y, double yRate); 
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

    // Material parameters
	int    soilType;	// Soil type = 1 for soft clay 
    double pult;		// Spring capacity
    double y50;			// y at 50% of pult
    double drag;		// ratio of max gap drag force to pult
	double yref;		// reference point for Near Field component
	double np;			// exponent for hardening shape of Near Field component
	double Elast;		// p/pult when yielding first occurs in virgin loading
	double nd;			// exponent for hardening shape of drag component
	double dashpot;     // dashpot on the far-field (elastic) component


  private:

	// Functions to get p & y for each component individually
	void getGap(double ylast, double dy, double dy_old);
	void getClosure(double ylast, double dy);
	void getDrag(double ylast, double dy);
	void getNearField(double ylast, double dy, double dy_old);
	void getFarField(double y);

	// Generated parameters or constants (not user input)
	double NFkrig;		// stiffness of the "rigid" portion of Near Field spring
	
    // Committed history variables for entire p-y material
    double Cy;			// Committed p
    double Cp;			// Committed y
    double Ctangent;	// Committed tangent

	// Trial history variables for entire p-y material
    double Ty;			// Trial p
    double Tp;			// Trial y
    double Ttangent;	// Trial tangent
	double TyRate;      // Trial velocity

	// Committed internal parameters for the NearField rigid-plastic component
	double CNFpinr;		//  p at start of current plastic loading cycle - right side
	double CNFpinl;		//                                              - left side
	double CNFyinr;		//  y at start of current plastic loading cycle - right side
	double CNFyinl;		//                                              - left side
	double CNF_p;		//  current p
	double CNF_y;		//  current y
	double CNF_tang;	//  tangent

	// Trial internal parameters for the NearField rigid-plastic component
	double TNFpinr;		//  p at start of current plastic loading cycle - right side
	double TNFpinl;		//                                              - left side
	double TNFyinr;		//  y at start of current plastic loading cycle - right side
	double TNFyinl;		//                                              - left side
	double TNF_p;		//  current p
	double TNF_y;		//  current y
	double TNF_tang;	//  tangent

	// Committed internal parameters for the Drag component
	double CDrag_pin;		//  p at start of current plastic loading cycle
	double CDrag_yin;		//  y at start of current plastic loading cycle
	double CDrag_p;			//  current p
	double CDrag_y;			//  current y
	double CDrag_tang;		//  tangent

	// Trial internal parameters for the Drag component
	double TDrag_pin;		//  p at start of current plastic loading cycle
	double TDrag_yin;		//  y at start of current plastic loading cycle
	double TDrag_p;			//  current p
	double TDrag_y;			//  current y
	double TDrag_tang;		//  tangent

	// Committed internal parameters for the Closure component
	double CClose_yleft;	//  left reference point
	double CClose_yright;	//  right reference point
	double CClose_p;		//  current p
	double CClose_y;		//  current y
	double CClose_tang;		//  tangent

	// Trial internal parameters for the Closure component
	double TClose_yleft;	//  left reference point
	double TClose_yright;	//  right reference point
	double TClose_p;		//  current p
	double TClose_y;		//  current y
	double TClose_tang;		//  tangent

	// Committed internal parameters for the Gap (Drag + Closure)
	double CGap_y;			//	y
	double CGap_p;			//  combined p
	double CGap_tang;		//  combined tangent

	// Trial internal parameters for the Gap (Drag + Closure)
	double TGap_y;			//	y
	double TGap_p;			//  combined p
	double TGap_tang;		//  combined tangent

	// Committed internal parameters for the Far Field component
	double CFar_y;			//  y
	double CFar_p;			//  current p
	double CFar_tang;       //  tangent

	// Trial internal parameters for the Far Field component
	double TFar_y;			//  y
	double TFar_p;			//  current p
	double TFar_tang;       //  tangent

	double initialTangent;
};


#endif
