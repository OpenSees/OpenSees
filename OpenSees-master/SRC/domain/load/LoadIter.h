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
// $Source: /usr/local/cvs/OpenSees/SRC/domain/load/LoadIter.h,v $
                                                                        
                                                                        
// File: ~/model/LoadIter.h
//
// Written: fmk 11/95
// Revised:
//
// Purpose: This file contains the class definition for LoadCaseIter.
// LoadIter is an iter for going through a particular load case
// and returning all the loads in it.

#ifndef LoadIter_h
#define LoadIter_h

#include "LoadCase.h"

class LoadIter
{
  public:
    LoadIter(const LoadCase &loadCase);
    Load *operator()(void);
    
  private:
    const LoadCase &myLoadCase;
    int currIndexLoads;
    int currIndexLoadCases;
};

// Interface Functions:
//

#endif
