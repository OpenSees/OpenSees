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
// $Source: /usr/local/cvs/OpenSees/SRC/analysis/model/simple/SimpleFE_Iter.cpp,v $
                                                                        
                                                                        
// File: ~/analysis/model/simple/SimpleFE_Iter.C
//
// Written: fmk 
// Created: Fri Sep 20 15:27:47: 1996
// Revision: A
//
// Description: This file contains the method definitions for class 
// SimpleFE_Iter. SimpleFE_Iter is a class for iterating through the 
// elements of an analysis model. 


#include <AnalysisModel.h>
#include <SimpleFE_Iter.h>


// SimpleFE_Iter(AnalysisModel &theModel):
//	constructor that takes the model, just the basic iter

SimpleFE_Iter::SimpleFE_Iter(AnalysisModel &theModel)
  :myModel(theModel), currIndex(0), numDone(0)
{
}


SimpleFE_Iter::~SimpleFE_Iter()
{
}    

void
SimpleFE_Iter::reset(void)
{
    currIndex = 0;
    numDone = 0;
}

FE_Element *
SimpleFE_Iter::operator()(void)
{
  // check if we still have elements in the model
  // if not return 0, indicating we are done
  if (numDone >= myModel.numFE_Ele)
    return 0;

  // search through domains ele list till we find the next element
  while ((currIndex < myModel.sizeEle) 
	 && (myModel.theFEs[currIndex] == 0))
      currIndex++;

  // if not at the end of the list return the element
  if (currIndex < myModel.sizeEle) {
      FE_Element *result = myModel.theFEs[currIndex];
      numDone++; currIndex++;
      return(result);
  }
  return (0);
}

    
    
