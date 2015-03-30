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

#ifndef ElementParameter_h
#define ElementParameter_h

#include <Parameter.h>
#include <ID.h>

class ElementParameter : public Parameter
{
 public:
  ElementParameter(int tag, 
		   int eleTag,
		   const char **argv, 
		   int argc);
  ElementParameter();
  ~ElementParameter();

  int update(int newValue); 
  int update(double newValue); 
  
  int addComponent(int, const char **argv, int argc);  

  void setDomain(Domain *theDomain);

  virtual int sendSelf(int commitTag, Channel &theChannel);  
  virtual int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);

 protected:
  
 private:
  Domain *theDomain;
  ID eleTags; // for each ele a tuple: eleTags, args index;
  char **argv;
  int argc;
  int argvSize;

  int numChannels;
  Channel **theChannels;
};

#endif
