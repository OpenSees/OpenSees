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
// $Date: 2005-11-28 22:07:24 $
// $Source: /usr/local/cvs/OpenSees/SRC/analysis/model/DOF_GrpIter.cpp,v $
                                                                        
// Written: fmk 
// Created: 10/05
// Revision: A
//

#include "DOF_GrpIter.h"

#include <DOF_Group.h>
#include <TaggedObjectIter.h>
#include <TaggedObjectStorage.h>


// DOF_GrpIter(SingleDomain &theDomain):
//	constructor that takes the model, just the basic iter

DOF_GrpIter::DOF_GrpIter(TaggedObjectStorage *theStorage)
  :myIter(&(theStorage->getComponents()))
{

}


DOF_GrpIter::~DOF_GrpIter()
{
}    

void
DOF_GrpIter::reset(void)
{
    myIter->reset();





}    


DOF_Group *
DOF_GrpIter::operator()(void)
{
    // check if we still have elements in the model
    // if not return 0, indicating we are done
    TaggedObject *theComponent = (*myIter)();
    if (theComponent == 0)
	return 0;
    else {
	DOF_Group *result = (DOF_Group *)theComponent;
	return result;
    }
}

    
    
