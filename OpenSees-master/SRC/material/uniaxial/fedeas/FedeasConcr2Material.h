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
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/fedeas/FedeasConcr2Material.h,v $
                                                                      
// Written: MHS
// Created: Jan 2001
//
// Description: This file contains the class definition for 
// FedeasConcr2Material. FedeasConcr2Material wraps the FEDEAS
// 1d material subroutine Concr_2.

#ifndef FedeasConcr2Material_h
#define FedeasConcr2Material_h

#include <FedeasMaterial.h>

class FedeasConcr2Material : public FedeasMaterial
{
  public:
    FedeasConcr2Material(int tag,
			 double fc, double ec, double fu, double eu,
			 double ratio, double ft, double Ets);
    FedeasConcr2Material(int tag, const Vector &data);
    FedeasConcr2Material(void);
    ~FedeasConcr2Material();

    const char *getClassType(void) const {return "FedeasConcr2Material";};

    double getInitialTangent(void);
    UniaxialMaterial *getCopy(void);

  protected:

  private:

};


#endif

