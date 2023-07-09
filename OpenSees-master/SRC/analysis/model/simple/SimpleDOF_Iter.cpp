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
// $Date: 2000-09-15 08:23:17 $
// $Source: /usr/local/cvs/OpenSees/SRC/analysis/model/simple/SimpleDOF_Iter.cpp,v $
                                                                        
                                                                        
// File: ~/analysis/model/simple/SimpleDOF_Iter.C
//
// Written: fmk 
// Created: Fri Sep 20 15:27:47: 1996
// Revision: A
//
// Description: This file contains the method definitions for class 
// SimpleDOF_Iter. SimpleDOF_Iter is a class for iterating through the 
// DOF_Groups of an simple analysis model. 

#include "SimpleDOF_Iter.h"
#include "AnalysisModel.h"


// SimpleDOF_Iter(AnalysisModel *theModel):
//	constructor that takes the model, just the basic iter

SimpleDOF_Iter::SimpleDOF_Iter(AnalysisModel &theModel)
  :myModel(theModel), currIndex(0), numDone(0)
{
}


SimpleDOF_Iter::~SimpleDOF_Iter()
{
}    

void
SimpleDOF_Iter::reset(void)
{
    currIndex = 0;
    numDone = 0;
}


DOF_Group *
SimpleDOF_Iter::operator()(void)
{
  // check if we still have elements in the model
  // if not return 0, indicating we are done

  if (numDone >= myModel.numDOF_Grp)
    return 0;

  // search through domains ele list till we find the next element
  while ((currIndex < myModel.sizeDOF) 
	 && (myModel.theDOFs[currIndex] == 0))
      currIndex++;

  // if not at the end of the list return the element
  if (currIndex < myModel.sizeDOF) {
      DOF_Group *result = myModel.theDOFs[currIndex];
      numDone++; currIndex++;
      return(result);
  }
  return (0);
}

    
    
