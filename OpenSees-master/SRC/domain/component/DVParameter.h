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
                                                                        
// $Revision: 1.7 $
// $Date: 2008-08-26 15:38:37 $
// $Source: /usr/local/cvs/OpenSees/SRC/domain/component/Parameter.h,v $

#ifndef DVParameter_h
#define DVParameter_h

#include <Parameter.h>
#include <string.h>

class DesignVariable;
class MovableObject;
class Channel;
class FEM_ObjectBroker;
class Domain;

class DVParameter : public Parameter
{
 public:
  DVParameter(int tag, DesignVariable *theDV, Parameter *theParam = 0);
  virtual ~DVParameter();
  
  virtual void Print(OPS_Stream &s, int flag =0);
  
  virtual int update(int newValue); 
  virtual int update(double newValue); 
  virtual int activate(bool active);
  virtual double getValue(void);
  virtual void setValue(double newValue);

  virtual bool isImplicit(void);
  virtual double getSensitivity(int index);
  virtual double getPerturbation(void);
  virtual const char *getType(void) {return "DesignVariable";}
  virtual int getPointerTag(void);

  virtual void setDomain(Domain *theDomain);
  virtual int sendSelf(int commitTag, Channel &theChannel);  
  virtual int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);

 protected:
  
 private:
  DesignVariable *myDV;
  Parameter *myParam;
  double currentValue;
};

#endif
