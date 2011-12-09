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
                                                                        
// $Revision: 1.8 $
// $Date: 2009-12-17 23:50:36 $
// $Source: /usr/local/cvs/OpenSees/SRC/recorder/response/Response.h,v $
                                                                        
// Written: MHS 
// Created: Oct 2000
//
// Description: This file contains the Response class interface

#ifndef Response_h
#define Response_h

class ID;
class Vector;
class Matrix;

#include <Information.h>

class Response
{
 public:
  Response(void);
  Response(int val);
  Response(double val);
  Response(const ID &val);
  Response(const Vector &val);
  Response(const Matrix &val);
  Response(const Vector &val1, const ID &val2);

  virtual ~Response();
  
  virtual int getResponse(void) = 0;
  virtual int getResponseSensitivity(int gradNumber) {return 0;}
  virtual Information &getInformation(void);

  virtual void Print(OPS_Stream &s, int flag = 0);
  virtual void Print(ofstream &s, int flag = 0);

 protected:
  Information myInfo;

 private:

};


#endif
