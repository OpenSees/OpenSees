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
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/fedeas/FedeasHyster1Material.h,v $
                                                                      
// Written: MHS
// Created: Jan 2001
//
// Description: This file contains the class definition for 
// FedeasHyster1Material. FedeasHyster1Material wraps the FEDEAS
// 1d material subroutine Hyster_1.

#ifndef FedeasHyster1Material_h
#define FedeasHyster1Material_h

#include <FedeasMaterial.h>

class FedeasHyster1Material : public FedeasMaterial
{
  public:
    FedeasHyster1Material(int tag,
			  double mom1p, double rot1p, double mom2p, double rot2p,
			  double mom1n, double rot1n, double mom2n, double rot2n,
			  double pinchX, double pinchY, double damfc1 = 0.0, double damfc2 = 0.0);
    FedeasHyster1Material(int tag, const Vector &d);
    FedeasHyster1Material(void);
    ~FedeasHyster1Material();

    const char *getClassType(void) const {return "FedeasHyster1Material";};

    double getInitialTangent(void);
    UniaxialMaterial *getCopy(void);

  protected:

  private:

};


#endif

