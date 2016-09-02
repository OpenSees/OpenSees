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
                                                                        
// $Revision: 1.7 $
// $Date: 2006-09-05 20:46:04 $
// $Source: /usr/local/cvs/OpenSees/SRC/domain/node/NodalLoad.h,v $
                                                                        
                                                                        
#ifndef NodalLoad_h
#define NodalLoad_h

// Written: fmk 
// Created: 11/96
// Revision: A
//
// Purpose: This file contains the class interface for NodalLoad.
// NodalLoad is a class for applying nodal loads to the model.

#include <Load.h>
#include <Node.h>
#include <Vector.h>

class NodalLoad : public Load
{
  public:
    NodalLoad(int classTag);
    NodalLoad(int tag, int node, int classTag);
    NodalLoad(int tag, int node, const Vector &load, bool isLoadConstant = false);
    ~NodalLoad();

    virtual void setDomain(Domain *newDomain);
    virtual int getNodeTag(void) const;
    virtual void applyLoad(double loadFactor);
    virtual void applyLoadSensitivity(double loadFactor);
    
    virtual int sendSelf(int commitTag, Channel &theChannel);
    virtual int recvSelf(int commitTag, Channel &theChannel, 
			 FEM_ObjectBroker &theBroker);
    
    virtual void Print(OPS_Stream &s, int flag =0);   
    
    // AddingSensitivity:BEGIN //////////////////////////////////////////
    int setParameter(const char **argv, int argc, Parameter &param);
    int            updateParameter(int parameterID, Information &info);
    int            activateParameter(int parameterID);
    const Vector & getExternalForceSensitivity(int gradNumber);
    // AddingSensitivity:END ///////////////////////////////////////////
  protected:

  private:
    int  myNode;        // tag indicating associated Node objects tag
    Node *myNodePtr;    // pointer to Node object on which load acts
    Vector *load;       // the reference load - pointer to new copy or 0
    bool  konstant;     // true if load is load factor independent
    // AddingSensitivity:BEGIN /////////////////////////////////////
    int parameterID;
    static Vector gradientVector;
    // AddingSensitivity:END ///////////////////////////////////////
};

#endif

