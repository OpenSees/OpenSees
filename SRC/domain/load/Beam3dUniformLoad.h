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
                                                                        
// $Revision: 1.3 $
// $Date: 2007-10-17 22:11:35 $
// $Source: /usr/local/cvs/OpenSees/SRC/domain/load/Beam3dUniformLoad.h,v $
                                                                        
#ifndef Beam3dUniformLoad_h
#define Beam3dUniformLoad_h

// Written: fmk 

// Purpose: This file contains the class definition for Beam3dUniformLoad.

#include <ElementalLoad.h>

class Beam3dUniformLoad : public ElementalLoad
{
  public:
    Beam3dUniformLoad(int tag, double wy, double wz, double wx, int eleTag);
		      
    Beam3dUniformLoad();    
    ~Beam3dUniformLoad();

    const Vector &getData(int &type, double loadFactor);

    int sendSelf(int commitTag, Channel &theChannel);  
    int recvSelf(int commitTag, Channel &theChannel,  FEM_ObjectBroker &theBroker);
    void Print(OPS_Stream &s, int flag =0);       

  protected:
	
  private:
    double wy;  // Transverse
    double wz;  // Transverse
    double wx;  // Axial
    static Vector data;
};

#endif

