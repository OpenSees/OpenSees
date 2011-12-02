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
// $Date: 2008-02-29 19:43:53 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/analysis/telm/PerformanceFunctionCoeff.cpp,v $
                                                                        
#include <PerformanceFunctionCoeff.h>
// for FEM_Object Broker to use
PerformanceFunctionCoeff::PerformanceFunctionCoeff
(int pTag, int pNodeID, int pDirection, double pCoeff)
:DomainComponent(pTag,0)
{
  Tag=pTag;
  NodeID=pNodeID;
  Direction=pDirection;
  Coefficient = pCoeff;
}    
PerformanceFunctionCoeff::~PerformanceFunctionCoeff()
{
}
int PerformanceFunctionCoeff::getTag(void) const
{
    return Tag;
}
int PerformanceFunctionCoeff::getNodeID(void) const
{
    return NodeID;
}
int PerformanceFunctionCoeff::getDirection(void) const
{
    return Direction;
}
double PerformanceFunctionCoeff::getCoefficient(void) const
{
    return Coefficient;
}
int PerformanceFunctionCoeff::sendSelf(int commitTag, Channel &theChannel){ return 0;}
int PerformanceFunctionCoeff::recvSelf(int commitTag, Channel &theChannel, 
			 FEM_ObjectBroker &theBroker){ return 0; }
void PerformanceFunctionCoeff::Print(OPS_Stream &s, int flag){}
