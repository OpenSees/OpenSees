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
// $Source: /usr/local/cvs/OpenSees/SRC/domain/load/Beam3dPointLoad.h,v $
                                                                        
#ifndef Beam3dPointLoad_h
#define Beam3dPointLoad_h

// Written: fmk 

// Purpose: This file contains the class definition for Beam3dPointLoad.

#include <ElementalLoad.h>

class Beam3dPointLoad : public ElementalLoad
{
  public:
    Beam3dPointLoad(int tag, double Py, double Pz, double x, int eleTag, double Pa = 0.0);
    Beam3dPointLoad();    
    ~Beam3dPointLoad();

    const Vector &getData(int &type, double loadFactor);

    int sendSelf(int commitTag, Channel &theChannel);  
    int recvSelf(int commitTag, Channel &theChannel,  FEM_ObjectBroker &theBroker);
    void Print(OPS_Stream &s, int flag =0);       

  protected:
	
  private:
    double Py;    // magnitude of the transverse load
    double Pz;    // magnitude of the transverse load
    double Px;    // magnitude of the axial load
    double x;     // relative distance (x/L) along length from end 1 of element
    static Vector data;
};

#endif

