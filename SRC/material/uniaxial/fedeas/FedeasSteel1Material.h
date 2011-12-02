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
// $Date: 2002-06-26 23:00:13 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/fedeas/FedeasSteel1Material.h,v $
                                                                      
// Written: MHS
// Created: Jan 2001
//
// Description: This file contains the class definition for 
// FedeasSteel1Material. FedeasSteel1Material wraps the FEDEAS
// 1d material subroutine Steel_1.

#ifndef FedeasSteel1Material_h
#define FedeasSteel1Material_h

#include <FedeasMaterial.h>

class FedeasSteel1Material : public FedeasMaterial
{
  public:
    FedeasSteel1Material(int tag,
		double fy, double E0, double b,
		double a1, double a2, double a3, double a4);

	// Constructor for no isotropic hardening
    FedeasSteel1Material(int tag,
		double fy, double E0, double b);

	FedeasSteel1Material(int tag, const Vector &d);
	FedeasSteel1Material(void);
    virtual ~FedeasSteel1Material();

	double getInitialTangent(void);
	UniaxialMaterial *getCopy(void);

  protected:

  private:

};


#endif

