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
                                                                        
// $Revision: 1.2 $
// $Date: 2010-04-02 23:38:05 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/joint/ElasticTubularJoint.h,v $
                                                                        
#ifndef ElasticTubularJoint_h
#define ElasticTubularJoint_h

// Written: Kia & Alanjari
//
// Description: This file contains the interface for the Tubular Joint Element class.
// It defines the class interface and the class attributes. 
// What: "@(#) ElasticTubularJoint.h, revA"

#include <Element.h>
#include <Matrix.h>
#include <Vector.h>

class Node;
class Channel;
class UniaxialMaterial;

#define ELE_TAG_ElasticTubularJoint 2519

class ElasticTubularJoint : public Element
{
  public:
    // constructors
    ElasticTubularJoint(int tag,int iNode , int jNode,
			double Brace_Diameter,
			double Brace_Angle,
			double e , 
			double Chord_Diameter , 
			double Chord_Thickness , 
			double Chord_Angle  );
    
    ElasticTubularJoint();    
    
    // destructor
    ~ElasticTubularJoint();

    
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

    // public methods to obtain stiffness, mass, damping and residual information    
    const Matrix &getTangentStiff(void);
    const Matrix &getInitialStiff(void);
    const Vector &getResistingForce(void);
    
    // public methods for output    
    int sendSelf(int commitTag, Channel &theChannel) ;
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker) ;
    void Print(OPS_Stream &s, int flag =0); 
    int displaySelf(Renderer &, int mode, float fact, const char **displayModes=0, int numModes=0);
 protected:
    
 private:
    double l;
    double cs , sn;
    double E;
    double braceD, braceangle;
    double chordD, chordT, chordangle;
    double InitLJFv;
    double InitLJFipb;
    double TangLJFv;
    double TangLJFipb;
    //UniaxialMaterial *theMaterial;

    Matrix k;    // matrix to return stiff
    Vector p;             // vector to return the resisting force
    Vector displacement;

    Node *theNodes[2];         //  contains the address of end nodes
    ID  connectedExternalNodes;	// contains the tags of the end nodes
};
#endif

