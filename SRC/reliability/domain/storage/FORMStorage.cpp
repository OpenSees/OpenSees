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
    gradientU = 0;
    gradientX = 0;
    beta = 0;
    firstCurvature = 0;
}


FORMStorage::~FORMStorage()
{
    if (alpha != 0)
        delete alpha;
    if (gradientU != 0)
        delete gradientU;
    if (gradientX != 0)
        delete gradientX;
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
    
    else if (strcmp(variable,"gradientUFORM") == 0) {
        gradientU = new Vector(*(theInfo.theVector));
    }
    
    else if (strcmp(variable,"gradientXFORM") == 0) {
        gradientX = new Vector(*(theInfo.theVector));
    }
    
    else if (strcmp(variable,"betaFORM") == 0) {

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
    
    else if (strcmp(variable,"gradientUFORM") == 0) {
        if (gradientU != 0) {
            theInfo.theType = VectorType;
            theInfo.setVector(*gradientU);
        }
        else
            return -1;
    }
    
    else if (strcmp(variable,"gradientXFORM") == 0) {
        if (gradientX != 0) {
            theInfo.theType = VectorType;
            theInfo.setVector(*gradientX);
        }
        else
            return -1;
    }
    
    else if (strcmp(variable,"betaFORM") == 0) {
        theInfo.theType = DoubleType;
        theInfo.setVector(beta);
    }
    
    else {
        opserr << "FORMStorage:: unknown variable " << variable << 
            " in getVariable()" << endln;
    }
    
    return 0;
}
