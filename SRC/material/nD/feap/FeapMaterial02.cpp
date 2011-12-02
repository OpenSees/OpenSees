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
// $Source: /usr/local/cvs/OpenSees/SRC/material/nD/feap/FeapMaterial02.cpp,v $
                                                                      
// Written: MHS
// Created: June 2001
//
// Description: This file contains the class definition for 
// FeapMaterial02.

#include <FeapMaterial02.h>
#include <Vector.h>

FeapMaterial02::FeapMaterial02(int tag,
			       double K, double G, double muK, double muG,
			       double lamK, double lamG, double theta):
  // 14 history variables and 7 material parameters
    FeapMaterial(tag, ND_TAG_FeapMaterial02, 14, 7)
{
  ud[0] = K;
  ud[1] = G;
  ud[2] = muK;
  ud[3] = muG;
  ud[4] = lamK;
  ud[5] = lamG;
  ud[6] = theta;
}

FeapMaterial02::FeapMaterial02(void):
  FeapMaterial(0, ND_TAG_FeapMaterial02, 14, 7)
{
  // Does nothing
}

FeapMaterial02::~FeapMaterial02(void)
{
  // Does nothing
}

int
FeapMaterial02::fillDArray(void)
{
  // Look in FEAP programmer's manual
  
  return 0;
}
