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
** ****************************************************************** */
                                                                        
// $Revision: 1.2 $
// $Date: 2007-06-06 19:36:19 $
// $Source: /usr/local/cvs/OpenSees/SRC/domain/component/MaterialStageParameter.h,v $

#ifndef MaterialStageParameter_h
#define MaterialStageParameter_h

#include <Parameter.h>
class Domain;

class MaterialStageParameter : public Parameter
{
 public:
  MaterialStageParameter(int tag, int materialTag);
  MaterialStageParameter();
  virtual ~MaterialStageParameter();

  virtual void Print(OPS_Stream &s, int flag =0);
  
  virtual int sendSelf(int commitTag, Channel &theChannel);  
  virtual int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);

  virtual void setDomain(Domain *theDomain);

 protected:
  
 private:
  Information theMatStageInfo;

  int theMaterialTag;
  int theParameterID;

  Domain *theDomain;
};

#endif
