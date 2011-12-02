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
// $Date: 2002-06-26 23:00:12 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/fedeas/FedeasConcr1Material.h,v $
                                                                      
// Written: MHS
// Created: Jan 2001
//
// Description: This file contains the class definition for 
// FedeasConcr1Material. FedeasConcr1Material wraps the FEDEAS
// 1d material subroutine Concr_1.

#ifndef FedeasConcr1Material_h
#define FedeasConcr1Material_h

#include <FedeasMaterial.h>

class FedeasConcr1Material : public FedeasMaterial
{
  public:
    FedeasConcr1Material(int tag,
		double fc, double ec, double fu, double eu);
	FedeasConcr1Material(int tag, const Vector &data);
	FedeasConcr1Material(void);
    virtual ~FedeasConcr1Material();

	double getInitialTangent(void);
	UniaxialMaterial *getCopy(void);

  protected:

  private:

};


#endif

