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
// $Date: 2003-03-04 00:48:09 $
// $Source: /usr/local/cvs/OpenSees/SRC/coordTransformation/CrdTransf2d.h,v $
                                                                        
                                                                        
// File: ~/CrdTransf/CrdTransf2d.h
//
// Written: Remo Magalhaes de Souza (rmsouza@ce.berkeley.edu)
// Created: 04/2000
// Revision: A
//
// Description: This file contains the class definition for 
// CrdTransf2d. CrdTransf2d provides the abstraction of spatial 
// coordinate transformation for a 2d frame. 

//
// What: "@(#) CrdTransf2d.h, revA"

#ifndef CrdTransf2d_h
#define CrdTransf2d_h

#include <CrdTransf.h>

class CrdTransf2d: public CrdTransf
{
  public:
    CrdTransf2d (int tag, int classTag);
    virtual ~CrdTransf2d();
   
    virtual CrdTransf2d *getCopy(void)=0;

  protected:
    
  private:
};

#endif
