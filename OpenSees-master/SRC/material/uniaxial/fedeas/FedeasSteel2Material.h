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
                                                                        
// $Revision: 1.3 $
// $Date: 2006-08-03 23:45:48 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/fedeas/FedeasSteel2Material.h,v $
                                                                      
// Written: MHS
// Created: Jan 2001
//
// Description: This file contains the class definition for 
// FedeasSteel2Material. FedeasSteel2Material wraps the FEDEAS
// 1d material subroutine Steel_2.

#ifndef FedeasSteel2Material_h
#define FedeasSteel2Material_h

#include <FedeasMaterial.h>

class FedeasSteel2Material : public FedeasMaterial
{
  public:
    FedeasSteel2Material(int tag,
			 double fy, double E0, double b,
			 double R0, double cR1, double cR2,
			 double a1, double a2, double a3, double a4);
    
    // Constructor for no isotropic hardening
    FedeasSteel2Material(int tag,
			 double fy, double E0, double b,
			 double R0, double cR1, double cR2);
    
    // Constructor for no isotropic hardening
    // Also provides default values for R0, cR1, and cR2
    FedeasSteel2Material(int tag,
			 double fy, double E0, double b);

    FedeasSteel2Material(int tag, const Vector &d);
    FedeasSteel2Material(void);
    ~FedeasSteel2Material();


    const char *getClassType(void) const {return "FedeasSteel2Material";};
    double getInitialTangent(void);
    UniaxialMaterial *getCopy(void);

  protected:

  private:

};


#endif

