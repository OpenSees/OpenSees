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
                                                                        
// $Revision: 1.1.1.1 $
// $Date: 2000-09-15 08:23:19 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/ElementalLoad.h,v $
                                                                        
                                                                        
#ifndef ElementalLoad_h
#define ElementalLoad_h

// File: ~/element/ElementalLoad.h
//
// Written: fmk 
// Created: 11/96
// Revision: A
//
// Purpose: This file contains the class definition for ElementalLoad.
// ElementalLoad is an abstract class.

#include <Load.h>

class ElementalLoad : public Load
{
  public:
    ElementalLoad(int elementTag, int tag, int classTag);
    ElementalLoad(int classTag);    
    ~ElementalLoad();

    virtual int getElementTag(void) const;

  protected:
	
  private:
    const int theElementTag; // tag indicating associated Element objects tag
};

#endif

