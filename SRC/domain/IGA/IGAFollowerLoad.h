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

// Written: Felipe Elgueta and Jose Abell
//          01.2021
//
// Description: This file contains the class definition for IGAFollowerLoad, a load class for 
//               applying point forces inside load patterns for IGA analysis.

#ifndef IGAFollowerLoad_h
#define IGAFollowerLoad_h

#include <ElementalLoad.h>

class IGAFollowerLoad : public ElementalLoad
{
  public:
    IGAFollowerLoad(int tag, double xi, double eta, double f1, double f2, double f3, int patchTag);
    IGAFollowerLoad();    
    ~IGAFollowerLoad();

    const Vector &getData(int &type, double loadFactor);

    int sendSelf(int commitTag, Channel &theChannel);  
    int recvSelf(int commitTag, Channel &theChannel,  FEM_ObjectBroker &theBroker);
    void Print(OPS_Stream &s, int flag =0);       

  protected:
	
  private:
    static Vector data;
    double xi;
    double eta;
    double f1;
    double f2;
    double f3;
};

#endif
