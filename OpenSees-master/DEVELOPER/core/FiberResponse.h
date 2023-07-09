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
// $Date: 2001-07-31 18:26:59 $
// $Source: /usr/local/cvs/OpenSees/SRC/recorder/response/FiberResponse.h,v $
                                                                        
// Written: MHS 
// Created: Oct 2000
//
// Description: This file contains the FiberResponse class interface

#ifndef FiberResponse_h
#define FiberResponse_h

#include <Response.h>
#include <Information.h>

class Fiber;

class ID;
class Vector;
class Matrix;

class FiberResponse : public Response
{
public:
	FiberResponse(Fiber *fib, int id);
	FiberResponse(Fiber *fib, int id, int val);
	FiberResponse(Fiber *fib, int id, double val);
	FiberResponse(Fiber *fib, int id, const ID &val);
	FiberResponse(Fiber *fib, int id, const Vector &val);
	FiberResponse(Fiber *fib, int id, const Matrix &val);
	~FiberResponse();

	int getResponse(void);

private:
	Fiber *theFiber;
	int responseID;
};

#endif
