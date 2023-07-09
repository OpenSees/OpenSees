/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 2001, The Regents of the University of California    **
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
** Reliability module developed by:                                   **
**   Terje Haukaas (haukaas@ce.berkeley.edu)                          **
**   Armen Der Kiureghian (adk@ce.berkeley.edu)                       **
**                                                                    **
** ****************************************************************** */
                                                                        
// $Revision: 1.17 $
// $Date: 2010-09-13 21:33:51 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/domain/components/LimitStateFunction.cpp,v $

//
// Written by: 
// Kevin Mackie (kmackie@mail.ucf.edu)
// Michael Scott (mhscott@engr.orst.edu)
//

#include <LimitStateFunction.h>
#include <FunctionEvaluator.h>
#include <classTags.h>

#include <map>
#include <string>
#include <stdlib.h>

LimitStateFunction::LimitStateFunction(int passedTag, const char *passedExpression)
:PerformanceFunction(passedTag, LIMIT_STATE_FUNCTION)
{
	
	int exprLen = strlen(passedExpression);
	theExpression = new char[exprLen+1];
	strcpy(theExpression,passedExpression);
}

LimitStateFunction::~LimitStateFunction()
{
	if (theExpression != 0)
		delete [] theExpression;

}


void
LimitStateFunction::Print(OPS_Stream &s, int flag)  
{
	s << "Limit State Function #" << this->getTag() << endln;
	s << "Expression: " << this->getExpression() << endln;
	s << endln;
}


const char *
LimitStateFunction::getExpression()
{
	return theExpression;
}


int
LimitStateFunction::addGradientExpression(const char *expression, int rvTag)
{
  map<int, string>::iterator theExpr;

  this->removeGradientExpression(rvTag);

  // Check if the expression is already in map
  theExpr = mapOfGradientExpressions.find(rvTag);
  
  // Not there, so add
  if (theExpr == mapOfGradientExpressions.end()) {
    mapOfGradientExpressions.insert(map<int,string>::value_type(rvTag,expression));
    
    // Check if successful
    theExpr = mapOfGradientExpressions.find(rvTag);
    if (theExpr == mapOfGradientExpressions.end()) {
      opserr << "LimitStateFunction::addGradientExpression -- map STL failed to add object with tag: "
	     << rvTag << endln;
      return -1;
    }
  }
  // Already there, give error
  else {
    opserr << "LimitStateFunction::addGradientExpression -- object with tag "
	   << rvTag << " already exists" << endln;    
    return -1;
  }
  
  return 0;
}


int
LimitStateFunction::removeGradientExpression(int rvTag)
{
  map<int, string>::iterator theExpr;

  // Check if the expression is already in map
  theExpr = mapOfGradientExpressions.find(rvTag);

  // If not there, do nothing
  if (theExpr == mapOfGradientExpressions.end())
    return 0;
  // Already there, so remove
  else {
    int ok = mapOfGradientExpressions.erase(rvTag);
    if (ok != 1) {
      opserr << "LimitStateFunction::removeGradientExpression -- map STL failed to remove object with tag: "
	     << rvTag << endln;
      return -1;
    }
  }

  return 0;
}


const char* 
LimitStateFunction::getGradientExpression(int rvTag) 
{
  map<int, string>::iterator theExpr;

  theExpr = mapOfGradientExpressions.find(rvTag);
  if (theExpr == mapOfGradientExpressions.end())
    return 0;
  else
    return (*theExpr).second.c_str();
}

