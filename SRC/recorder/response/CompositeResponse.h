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
// $Date: 2010-05-13 00:15:36 $
// $Source: /usr/local/cvs/OpenSees/SRC/recorder/response/CompositeResponse.h,v $
                                                                        
// Written: fmk
// Created: 05/10
//
// Description: This file contains the CompositeResponse class interface.
// A CompositeResponse is a container holding a number of response objects.

#ifndef CompositeResponse_h
#define CompositeResponse_h

class ID;
class Vector;
class Matrix;

#include <Response.h>

class CompositeResponse: public Response
{
 public:
  CompositeResponse(void);
  virtual ~CompositeResponse();

  int addResponse(Response *);  
  int getResponse(void);

 protected:

 private:
  Response **theResponses;
  int numResponses;
};


#endif
