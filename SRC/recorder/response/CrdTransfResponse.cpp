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
                                                                        
// $Revision: 1.4 $
// $Date: 2007-03-02 00:12:50 $
// $Source: /usr/local/cvs/OpenSees/SRC/recorder/response/CrdTransfResponse.cpp,v $
                                                                        
// Written: MHS 
// Created: Oct 2000
//
// Description: This file contains the CrdTransfResponse class implementation

#include <CrdTransfResponse.h>
#include <CrdTransf.h>

CrdTransfResponse::CrdTransfResponse(CrdTransf *mat, int id):
  Response(), theCrdTransf(mat), responseID(id)
{

}

CrdTransfResponse::CrdTransfResponse(CrdTransf *mat, int id, int val):
  Response(val), theCrdTransf(mat), responseID(id)
{

}

CrdTransfResponse::CrdTransfResponse(CrdTransf *mat, int id, double val):
  Response(val), theCrdTransf(mat), responseID(id)
{

}

CrdTransfResponse::CrdTransfResponse(CrdTransf *mat, int id, const ID &val):
  Response(val), theCrdTransf(mat), responseID(id)
{

}

CrdTransfResponse::CrdTransfResponse(CrdTransf *mat, int id, const Vector &val):
  Response(val), theCrdTransf(mat), responseID(id)
{

}

CrdTransfResponse::CrdTransfResponse(CrdTransf *mat, int id, const Matrix &val):
  Response(val), theCrdTransf(mat), responseID(id)
{

}

CrdTransfResponse::~CrdTransfResponse()
{

}

int
CrdTransfResponse::getResponse(void)
{
  if (theCrdTransf != 0)
    return theCrdTransf->getResponse(responseID, myInfo);
  else
    return 0;
}

