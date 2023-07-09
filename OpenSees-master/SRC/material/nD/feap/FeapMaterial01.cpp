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
                                                                        
// $Revision: 1.1 $
// $Date: 2002-10-29 20:26:23 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/nD/feap/FeapMaterial01.cpp,v $
                                                                      
// Written: MHS
// Created: June 2001
//
// Description: This file contains the class definition for 
// FeapMaterial01.

#include <FeapMaterial01.h>
#include <Vector.h>

FeapMaterial01::FeapMaterial01(int tag, double E, double nu, double rho):
  // 0 history variables and 2 material parameters
  FeapMaterial(tag, ND_TAG_FeapMaterial01, 0, 2, rho)
{
  ud[0] = E;
  ud[1] = nu;
}

FeapMaterial01::FeapMaterial01(void):
  FeapMaterial(0, ND_TAG_FeapMaterial01, 0, 2)
{
  // Does nothing
}

FeapMaterial01::~FeapMaterial01(void)
{
  // Does nothing
}

int
FeapMaterial01::fillDArray(void)
{
  // Look in FEAP programmer's manual
  d[0] = ud[0];
  d[1] = ud[1];
  d[3] = rho;
  
  return 0;
}
