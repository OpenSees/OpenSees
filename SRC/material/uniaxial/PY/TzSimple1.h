/* *********************************************************************
**    Module:	TzSimple1.h 
**
**    Purpose:	Provide a simple t-z spring for OpenSees
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
// $Date: 2002/1/19
// $Source: /OpenSees/SRC/material/uniaxial/TzSimple1.h

#ifndef TZSIMPLE1_H
#define TZSIMPLE1_H

// Written: RWB
// Created: Jan 2002
// tested and checked: Boris Jeremic (jeremic@ucdavis.edu) Spring 2002
//
// Description: This file contains the class definition for TzSimple1.
// 

#include <UniaxialMaterial.h>


class TzSimple1 : public UniaxialMaterial
{
  public:
    TzSimple1(int tag, int classtag, int tzType, double tult, double z50, double dashpot);
    TzSimple1();
    ~TzSimple1();

    const char *getClassType(void) const {return "TzSimple1";};

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
	int    tzType;		// tz relation type 
    double tult;		// Material capacity
    double z50;			// z at 50% of tult
	double zref;		// reference point for Near Field spring
	double np;			// exponent for hardening shape of Near Field Spring
	double dashpot;     // dashpot on the far-field (elastic) component

  private:

	// Functions to get t & z for elastic & plastic components
	void getNearField(double zlast, double dz, double dz_old);
	void getFarField(double z);

   // Committed history variables for entire t-z material
    double Cz;			// Committed t
    double Ct;			// Committed z
    double Ctangent;	// Committed tangent

	// Trial history variables for entire t-z material
    double Tz;			// Trial z
    double Tt;			// Trial t
    double Ttangent;	// Trial tangent
	double TzRate;      // Trial velocity

	// Committed internal parameters for the NearField plastic component
	double CNF_tin;			//  t at start of current plastic loading cycle
	double CNF_zin;			//  z at start of current plastic loading cycle
	double CNF_t;			//  current t
	double CNF_z;			//  current z
	double CNF_tang;		//  tangent

	// Trial internal parameters for the NearField plastic component
	double TNF_tin;			//  t at start of current plastic loading cycle
	double TNF_zin;			//  z at start of current plastic loading cycle
	double TNF_t;			//  current t
	double TNF_z;			//  current z
	double TNF_tang;		//  tangent

	// Committed internal parameters for the Far Field elastic component
	double CFar_z;			//  current z
	double CFar_t;			//  current t
	double CFar_tang;       //  tangent

	// Trial internal parameters for the Far Field component
	double TFar_z;			//  current z
	double TFar_t;			//  current t
	double TFar_tang;       //  tangent
	
	double initialTangent;

};


#endif
