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
                                                                        
// $Revision: 1.1 $
// $Date: 2007-06-06 19:34:31 $
// $Source: /usr/local/cvs/OpenSees/SRC/domain/component/MatParameter.h,v $


// written: fmk

#ifndef MatParameter_h
#define MatParameter_h

#include <Parameter.h>
class Domain;

class MatParameter : public Parameter
{
 public:
  MatParameter(int tag, int materialTag, const char *materialParameterName);
  MatParameter();
  virtual ~MatParameter();

  virtual void Print(OPS_Stream &s, int flag =0);
  
  virtual int sendSelf(int commitTag, Channel &theChannel);  
  virtual int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);

  virtual void setDomain(Domain *theDomain);

 protected:
  
 private:
  Information theMatInfo;
  char *theParameterName;

  int theMaterialTag;
  int theParameterID;

  Domain *theDomain;
};

#endif
