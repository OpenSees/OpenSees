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
// $Source: /usr/local/cvs/OpenSees/SRC/domain/load/LoadIter.cpp,v $
                                                                        
                                                                        
// File: ~/model/LoadIter.C
//
// Written: fmk 11/95
// Revised:
//
// Purpose: This file contains the method definitions for class LoadIter.
// LoadIter is an abstract class.


#include "LoadIter.h"

// LoadIter(LoadCase &loadcase):
//	constructor that takes the nodes associated with a load

LoadIter::LoadIter(const LoadCase &loadcase)
  :myLoadCase(loadcase), currIndexLoads(0), currIndexLoadCases(0)
{

}


Load *
LoadIter::operator()(void)
{
/** WARNING - HAVE TO CHANGE FOR LOAD CASES */

 // if we still have loads from myLoads to return
  if (currIndexLoads == myLoadCase.numLoads)
    return 0;
  else
    return myLoadCase.myLoads[currIndexLoads++];
}
  


    
    
