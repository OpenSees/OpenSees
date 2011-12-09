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
                                                                        
// $Revision: 1.3 $
// $Date: 2001-07-31 18:26:58 $
// $Source: /usr/local/cvs/OpenSees/SRC/recorder/response/FiberResponse.cpp,v $
                                                                        
// Written: MHS 
// Created: Oct 2000
//
// Description: This file contains the FiberResponse class implementation

#include <FiberResponse.h>
#include <Fiber.h>

FiberResponse::FiberResponse(Fiber *fib, int id):
Response(), theFiber(fib), responseID(id)
{

}

FiberResponse::FiberResponse(Fiber *fib, int id, int val):
Response(val), theFiber(fib), responseID(id)
{

}

FiberResponse::FiberResponse(Fiber *fib, int id, double val):
Response(val), theFiber(fib), responseID(id)
{

}

FiberResponse::FiberResponse(Fiber *fib, int id, const ID &val):
Response(val), theFiber(fib), responseID(id)
{

}

FiberResponse::FiberResponse(Fiber *fib, int id, const Vector &val):
Response(val), theFiber(fib), responseID(id)
{

}

FiberResponse::FiberResponse(Fiber *fib, int id, const Matrix &val):
Response(val), theFiber(fib), responseID(id)
{

}


FiberResponse::~FiberResponse()
{

}

int
FiberResponse::getResponse(void)
{
	return theFiber->getResponse(responseID, myInfo);
}
