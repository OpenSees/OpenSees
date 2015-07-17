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
                                                                        
// $Revision: 1.0 $
// $Date: 2010/1/8 21:00:00 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/limitState/WrapperLimitCurve.cpp,v $

// Written: MRL                                                                         

#include <WrapperLimitCurve.h>
#include <string.h>

extern modelState theModelState;

WrapperLimitCurve::WrapperLimitCurve(const char *name, limCrvObject *theLimCrv_)
  :LimitCurve(theLimCrv_->tag,LIMCRV_TAG_WrapperLimitCurve),
  funcName(0),
  theLimCrv(theLimCrv_),
  springForce(0.0), Kdeg(0.0), Fres(0.0), DR(0.0)
{
  /*opserr << "WrapperLimitCurve::WrapperLimitCurve() " << theLimCrv->tag << endln; */

  funcName = new char[strlen(name)+1];
  if (funcName != 0)
    strcpy(funcName, name);

  int isw = ISW_FORM_TANG_AND_RESID;
  int error = 0;
  theLimCrv->limCrvFunctPtr(theLimCrv, &theModelState, &springForce, &Kdeg, &Fres, &isw, &error);

  //  theLimCrv->limCrvObjectPtr = this;
}
WrapperLimitCurve::WrapperLimitCurve()
	:LimitCurve(0, LIMCRV_TAG_WrapperLimitCurve),
	funcName(0),
	theLimCrv(0),
	springForce(0.0), Kdeg(0.0), Fres(0.0), DR(0.0)
{
  
}

// destructor
WrapperLimitCurve::~WrapperLimitCurve()
{
	/*opserr << "WrapperLimitCurve::~WrapperLimitCurve()\n";*/

	if (funcName != 0)
		delete [] funcName;

	if (theLimCrv->theParam != 0)
		delete [] theLimCrv->theParam;

	if (theLimCrv->cState != 0)
		delete [] theLimCrv->cState;

	if (theLimCrv->tState != 0)
		delete [] theLimCrv->tState;

	delete theLimCrv;
}

int 
WrapperLimitCurve::checkElementState(double theSpringForce)
{
	int isw = ISW_FORM_TANG_AND_RESID;
	int error = 0;
	springForce = theSpringForce;
	theLimCrv->limCrvFunctPtr(theLimCrv, &theModelState, &springForce, &Kdeg, &Fres, &isw, &error);

	return error;
}

double 
WrapperLimitCurve::getDegSlope(void)
{
	return Kdeg;
}
double 
WrapperLimitCurve::getResForce(void)
{
	return Fres;
}
double 
WrapperLimitCurve::getUnbalanceForce(void)
{
	return 0.0;
}
double 
WrapperLimitCurve::findLimit(double DR)
{
	return 0.0;
}

int 
WrapperLimitCurve::revertToStart (void)
{
  int isw = ISW_REVERT_TO_START;
  int error = 0;
  theLimCrv->limCrvFunctPtr(theLimCrv, &theModelState, &springForce, &Kdeg , &Fres, &isw, &error);
  return error;
}

LimitCurve *
WrapperLimitCurve::getCopy (void) 
{
    limCrvObject *theLimCrvObject = new limCrvObject;
    theLimCrvObject->tag = theLimCrv->tag;
    theLimCrvObject->nParam = theLimCrv->nParam;
    theLimCrvObject->nState = theLimCrv->nState;

    OPS_AllocateLimitCurve(theLimCrvObject);

    for (int i=0; i<theLimCrv->nParam; i++)
      theLimCrvObject->theParam[i] = theLimCrv->theParam[i];

    for (int i=0; i<theLimCrv->nState; i++) {
      theLimCrvObject->cState[i] = theLimCrv->cState[i];
      theLimCrvObject->tState[i] = theLimCrv->tState[i];
    }

    theLimCrvObject->limCrvFunctPtr = theLimCrv->limCrvFunctPtr;

    WrapperLimitCurve *theResult = new WrapperLimitCurve(funcName, theLimCrvObject);
    return theResult;
}

int 
WrapperLimitCurve::sendSelf(int commitTag, Channel &theChannel)
{
  return -1;
}
int 
WrapperLimitCurve::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  return -1;
}
void 
WrapperLimitCurve::Print(OPS_Stream &s, int flag)
{
  s << "WrapperLimitCurve - wrapping function  limCrvTag: " << theLimCrv->tag << endln;
}



