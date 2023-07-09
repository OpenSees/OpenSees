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
// $Date: 2006-08-04 18:11:35 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/nD/feap/FeapMaterial01.h,v $
                                                                      
// Written: MHS
// Created: June 2001
//
// Description: This file contains the class definition for 
// FeapMaterial01.

#ifndef FeapMaterial01_h
#define FeapMaterial01_h

#include <FeapMaterial.h>

class FeapMaterial01 : public FeapMaterial
{
 public:
  FeapMaterial01(int tag, double E, double nu, double rho = 0.0);
  FeapMaterial01(void);
  ~FeapMaterial01();

  const char *getClassType(void) const {return "FeapMaterial01";};
  
 protected:
  int fillDArray(void);
  
 private:
  
};

#endif
