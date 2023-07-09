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
                                                                                                                                                
                                                                        
// Written: Kevin Mackie
//

#include <ReliabilityStorage.h>
#include <Information.h>

ReliabilityStorage::ReliabilityStorage()
{

}


ReliabilityStorage::~ReliabilityStorage()
{

}


static char unknownClassType[] = {"UnknownReliabilityStorage"};

const char *
ReliabilityStorage::getClassType(void) const
{
  return unknownClassType;
}


int 
ReliabilityStorage::setVariable(const char *variable, Information &theInfo)
{
  return -1;
}


int 
ReliabilityStorage::getVariable(const char *variable, Information &theInfo)
{
  return -1;
}
