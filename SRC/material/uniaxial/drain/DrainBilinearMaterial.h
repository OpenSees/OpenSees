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
// $Date: 2008-04-14 21:28:21 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/drain/DrainBilinearMaterial.h,v $
                                                                      
// Written: MHS
// Created: June 2001
//
// Description: This file contains the class definition for 
// DrainBilinearMaterial.

#ifndef DrainBilinearMaterial_h
#define DrainBilinearMaterial_h

#include <DrainMaterial.h>

class DrainBilinearMaterial : public DrainMaterial
{
  public:
  DrainBilinearMaterial(int tag,
			double E, double fyp, double fyn, double alpha,
			double ecaps, double ecapk, double ecapa, double ecapd,
			double cs, double ck, double ca, double cd,
		double capSlope, double capDispP, double capDispN, double res, double beto = 0.0);
  DrainBilinearMaterial(int tag, const Vector &input, double beto = 0.0);
  DrainBilinearMaterial(void);
  virtual ~DrainBilinearMaterial();

  const char *getClassType(void) const {return "DrainBilinearMaterial";};

  int revertToStart(void);
  
  UniaxialMaterial *getCopy(void);

  protected:

  private:

};


#endif

