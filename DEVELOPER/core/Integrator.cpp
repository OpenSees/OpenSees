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
// $Source: /usr/local/cvs/OpenSees/SRC/analysis/integrator/Integrator.cpp,v $
                                                                        
                                                                        
// File: ~/analysis/integrator/Integrator.C
// 
// Written: fmk 
// Created: 11/96
// Revision: A
//
// Description: This file contains the implementation of Integrator.
//
// What: "@(#) Integrator.C, revA"

#include <Integrator.h>
#include <Vector.h>

 Integrator::Integrator(int clasTag)
:MovableObject(clasTag)
{
SensitivityKey=0;
}




//////////////////////////////////



Integrator::~Integrator()
{

}

int
Integrator::domainChanged()
{
    return 0;
}


int
Integrator::formSensitivityRHS(int gradNum)
{
    return 0;
}

int
Integrator::formIndependentSensitivityRHS()
{
    return 0;
}

int
Integrator::saveSensitivity(const Vector& v, int gradNum, int numGrads)
{
    return 0;
}

int
Integrator::commitSensitivity(int gradNum, int numGrads)
{
    return 0;
}
///////////////////////Abbas///////////////////////////////////////////

 int Integrator:: formEleTangentSensitivity(FE_Element 
      *theEle, int gradNumber)
{
return 0;
}

double 
Integrator::getLambdaSensitivity(int gradNumber)
{
return 0.0;


}


int
Integrator::computeSensitivities()
{

//should not be called
return 0;

}
int 
Integrator::sensitivityDomainChanged()
{
// I do not think I need it
return 0;
}


bool 
Integrator::shouldComputeAtEachStep(void)
{
  return (analysisTypeTag == 1);
}
bool
Integrator::computeSensitivityAtEachIteration()
{
return false ;
}

// the following function is to activate the sensitivity key, and it will be called in the commands.cpp
bool
Integrator::activateSensitivityKey()
{
SensitivityKey=true;

return SensitivityKey;

}


 ////////////////////////Abbas/////////////////////////////////////
