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
                                                                        

// Written: MHS 
//
// Description: This file contains the FORMStorage class implementation

#include <FORMStorage.h>
#include <Information.h>
#include <string.h>

FORMStorage::FORMStorage()
{
    alpha = 0;
}


FORMStorage::~FORMStorage()
{
    if (alpha != 0)
        delete alpha;
}


const char *
FORMStorage::getClassType(void) const
{
    return "FORM";
}


int 
FORMStorage::setVariable(const char *variable, Information &theInfo)
{
    if (strcmp(variable,"alphaFORM") == 0) {
        alpha = new Vector(*(theInfo.theVector));
    }
    else if (strcmp(variable,"gradientU") == 0) {
        gradient = new Vector(*(theInfo.theVector));
    }
    else if (strcmp(variable,"beta") == 0) {
        
    }
    else {
        opserr << "FORMStorage:: unknown variable " << variable << 
            " in setVariable()" << endln;
    }
        
    return 0;
}


int 
FORMStorage::getVariable(const char *variable, Information &theInfo)
{
    if (strcmp(variable,"alphaFORM") == 0) {
        if (alpha != 0) {
            theInfo.theType = VectorType;
            theInfo.setVector(*alpha);
        }
        else
            return -1;
    }
    else if (strcmp(variable,"gradient") == 0) {
        if (gradient != 0) {
            theInfo.theType = VectorType;
            theInfo.setVector(*gradient);
        }
        else
            return -1;
    }
    else if (strcmp(variable,"beta") == 0) {
        
    }
    else {
        opserr << "FORMStorage:: unknown variable " << variable << 
            " in getVariable()" << endln;
    }
    
    return 0;
}
