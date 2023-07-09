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
// $Date: 2004-09-01 03:54:28 $
// $Source: /usr/local/cvs/OpenSees/SRC/damage/DamageResponse.h,v $
                                                                        
// Written: Arash Altoontash, Gregory Deierlein
// Created: Aug 2002
//
// Description: This file contains the DamageResponse class interface

#ifndef DamageResponse_h
#define DamageResponse_h

#include <Information.h>
#include <Response.h>

class ID;
class Vector;
class Matrix;
class Tensor;
class DamageModel;

class DamageResponse : public Response
{
public:
  DamageResponse(DamageModel *dmg, int id);
  DamageResponse(DamageModel *dmg, int id, int val);
  DamageResponse(DamageModel *dmg, int id, double val);
  DamageResponse(DamageModel *dmg, int id, const ID &val);
  DamageResponse(DamageModel *dmg, int id, const Vector &val);
  DamageResponse(DamageModel *dmg, int id, const Matrix &val);
  ~DamageResponse();
  
  int getResponse(void);
  
 private:
  DamageModel *theDamage;
  int responseID;
};

#endif
