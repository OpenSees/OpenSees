/* *********************************************************************
**    Module:	QzSimple1.h 
**
**    Purpose:	Provide a simple Q-z material for OpenSees
**
**
**    Developed by Ross W. Boulanger
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


class QzSimple1 : public UniaxialMaterial
{
  public:
    QzSimple1(int tag, int qzType, double Qult, double z50, double suction,
		      double dashpot);
    QzSimple1();
    ~QzSimple1();

    const char *getClassType(void) const {return "QzSimple1";};

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
