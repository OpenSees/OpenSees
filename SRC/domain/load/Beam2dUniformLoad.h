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
// $Date: 2007-10-17 22:11:35 $
// $Source: /usr/local/cvs/OpenSees/SRC/domain/load/Beam2dUniformLoad.h,v $
                                                                        
#ifndef Beam2dUniformLoad_h
#define Beam2dUniformLoad_h

// Written: fmk 

// Purpose: This file contains the class definition for Beam2dUniformLoad.

#include <ElementalLoad.h>

class Beam2dUniformLoad : public ElementalLoad
{
   public:
    Beam2dUniformLoad(int tag, double wTrans, double wAxial, int eleTag);
    Beam2dUniformLoad();    
    ~Beam2dUniformLoad();

    const Vector &getData(int &type, double loadFactor);

    int sendSelf(int commitTag, Channel &theChannel);  
    int recvSelf(int commitTag, Channel &theChannel,  FEM_ObjectBroker &theBroker);
    void Print(OPS_Stream &s, int flag =0);       

    int setParameter(const char **argv, int argc, Parameter &param);
    int updateParameter(int parameterID, Information &info);
    int activateParameter(int paramID);

    const Vector &getSensitivityData(int gradNumber);

  protected:
	
  private:
    double wTrans;
    double wAxial;
    static Vector data;

    int parameterID;
};

#endif

