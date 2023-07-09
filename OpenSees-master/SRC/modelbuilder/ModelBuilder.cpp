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
// $Date: 2000-09-15 08:23:23 $
// $Source: /usr/local/cvs/OpenSees/SRC/modelbuilder/ModelBuilder.cpp,v $
                                                                        
                                                                        
// File: ~/model/ModelBuilder.C
//
// Written: fmk 
// Created: Mon Sept 15 14:47:47: 1996
// Revision: A
//
// Description: This file contains the class definition for ModelBuilder
// ModelBuilder is a class used to model a structure with a plane frame. 
// The object creates the components of the model and adds these to the
// Domain with which it is associated.



#include <ModelBuilder.h>
#include <Domain.h>

//  ModelBuilderModel(Domain &theDomain);
//	constructor
ModelBuilder::ModelBuilder(Domain &theDomain)
:myDomain(&theDomain)
{

}

ModelBuilder::~ModelBuilder()
{
    
}

Domain *
ModelBuilder::getDomainPtr(void) const 
{
    return myDomain;
}
    

