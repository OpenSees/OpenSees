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
// $Date: 2009-04-17 23:02:41 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/frictionBearing/frictionModel/FrictionResponse.h,v $

// Written: Andreas Schellenberg (andreas.schellenberg@gmx.net)
// Created: 02/06
// Revision: A
//
// Description: This file contains the FrictionResponse class interface
//
// What: "@(#) FrictionResponse.h, revA"

#ifndef FrictionResponse_h
#define FrictionResponse_h

#include <Response.h>
#include <Information.h>

class FrictionModel;

class ID;
class Vector;
class Matrix;
class Tensor;

class FrictionResponse : public Response
{
public:
    FrictionResponse(FrictionModel *frn, int id);
    FrictionResponse(FrictionModel *frn, int id, int val);
    FrictionResponse(FrictionModel *frn, int id, double val);
    FrictionResponse(FrictionModel *frn, int id, const ID &val);
    FrictionResponse(FrictionModel *frn, int id, const Vector &val);
    FrictionResponse(FrictionModel *frn, int id, const Matrix &val);
    FrictionResponse(FrictionModel *frn, int id, const Tensor &val);
    ~FrictionResponse();
    
    int getResponse(void);
    
private:
    FrictionModel *theFriction;
    int responseID;
};

#endif
