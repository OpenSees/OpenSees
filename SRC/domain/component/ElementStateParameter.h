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
                                                                        
#ifndef ElementStateParameter_h
#define ElementStateParameter_h

#include <Parameter.h>
#include <string.h>

class Node;
class MovableObject;
class Channel;
class FEM_ObjectBroker;
class Domain;

class ElementStateParameter : public Parameter
{
 public:
  ElementStateParameter(double value, 
			const char **argv, 
			int argc, 
			int flag, 
			ID *theEle = 0);

  ElementStateParameter();
  ~ElementStateParameter();
  
  void Print(OPS_Stream &s, int flag =0);

  void setDomain(Domain *theDomain);
  int sendSelf(int commitTag, Channel &theChannel);  
  int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);

 protected:
  
 private:
  double currentValue;

  ID *theEleIDs;
  int flag; // 0 all ele, 1 range ele, 2 list ele

  char **argv; // argc to go to the elements
  int argc;

  int fromFree;
};

#endif
