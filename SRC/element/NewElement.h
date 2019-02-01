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
// $Date: 2007-02-02 01:30:47 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/NewElement.h,v $
                                                                        
#ifndef NewElement_h
#define NewElement_h

// Written: fmk 
// Created: 08/01
//
// Description: This file contains the class definition for NewElement. 
//
// What: "@(#) NewElement.h, revA"

#include <Element.h>
#include <Matrix.h>
#include <Vector.h>

class Channel;
class UniaxialMaterial;

#define ELE_TAG_NewElement 100001

class NewElement : public Element
{
  public:
    NewElement(int tag);
    NewElement();    
    ~NewElement();

	const char *getClassType(void) const {return "NewElement";};

    // public methods to obtain information about dof & connectivity    
    int getNumExternalNodes(void) const;
    const ID &getExternalNodes(void);
    Node **getNodePtrs(void);
    int getNumDOF(void);	
    void setDomain(Domain *theDomain);

    // public methods to set the state of the element    
    int commitState(void);
    int revertToLastCommit(void);        
    int revertToStart(void);        
    int update(void);
    
    // public methods to obtain stiffness, mass, damping and 
    // residual information    
    const Matrix &getTangentStiff(void);
    const Matrix &getInitialStiff(void);    

    void zeroLoad(void);	
    int addLoad(const Vector &addP);
    int addInertiaLoadToUnbalance(const Vector &accel);
    const Vector &getResistingForce(void);
    const Vector &getResistingForceIncInertia(void);            

    // public methods for element output
    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, 
		 FEM_ObjectBroker &theBroker);
    int displaySelf(Renderer &, int mode, float fact, const char **displayModes=0, int numModes=0);
    void Print(OPS_Stream &s, int flag =0);    

    Response *setResponse(const char **argv, int argc, OPS_Stream &);
    int getResponse(int responseID, Information &eleInformation);

    int setParameter (const char **argv, int argc, Parameter &param);
    int updateParameter (int parameterID, Information &info);

  protected:
    
  private:
    ID  connectedExternalNodes;   // contains the tags of the end nodes
    Matrix theMatrix;             // matrix to return stiff, damp & mass
    Vector theVector;             // vector to return the residual
};

#endif




