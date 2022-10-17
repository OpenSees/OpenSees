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
                                                                        
// $Revision: 1.5 $
// $Date: 2007-03-02 00:12:50 $
// $Source: /usr/local/cvs/OpenSees/SRC/recorder/response/CrdTransfResponse.h,v $
                                                                        
// Written: MHS 
// Created: Oct 2000
//
// Description: This file contains the CrdTransfResponse class interface

#ifndef CrdTransfResponse_h
#define CrdTransfResponse_h

#include <Response.h>
#include <Information.h>

class CrdTransf;
class CrdTransfState;

class ID;
class Vector;
class Matrix;

class CrdTransfResponse : public Response
{
 public:
  CrdTransfResponse(CrdTransf *mat, int id);
  CrdTransfResponse(CrdTransf *mat, int id, int val);
  CrdTransfResponse(CrdTransf *mat, int id, double val);
  CrdTransfResponse(CrdTransf *mat, int id, const ID &val);
  CrdTransfResponse(CrdTransf *mat, int id, const Vector &val);
  CrdTransfResponse(CrdTransf *mat, int id, const Matrix &val);

  ~CrdTransfResponse();
  
  int getResponse(void);

private:
  CrdTransf *theCrdTransf;
  int responseID;
};

#endif
