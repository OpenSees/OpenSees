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
                                                                        
// $Revision: 1.11 $
// $Date: 2010/02/04 01:12:57 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/truss/CorotTruss2.h,v $

#ifndef CorotTruss2_h
#define CorotTruss2_h

// Written Y.Lu and M.Panagiotou 2013
//  minor mod to CorotTruss written by  MHS  2001

#include <Element.h>
#include <Matrix.h>
#include <Vector.h>
#include <ConcretewBeta.h>

class Node;
class Channel;
class UniaxialMaterial;

class CorotTruss2 : public Element
{
  public:
    CorotTruss2(int tag, int dim,
	       int Nd1, int Nd2, int oNd1, int oNd2, 
	       UniaxialMaterial &theMaterial,
	       double A, double rho=0.0);
    
    CorotTruss2();    
    ~CorotTruss2();

    const char *getClassType(void) const {return "CorotTruss2";};

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
    const Matrix &getMass(void);    

    void zeroLoad(void);	
    int addLoad(ElementalLoad *theLoad, double loadFactor);
    int addInertiaLoadToUnbalance(const Vector &accel);

    const Vector &getResistingForce(void);
    const Vector &getResistingForceIncInertia(void);            

    // public methods for element output
    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);
    int displaySelf(Renderer &, int mode, float fact, const char **displayModes=0, int numModes=0);
    void Print(OPS_Stream &s, int flag =0);    

    Response *setResponse(const char **argv, int argc, OPS_Stream &s);
    int getResponse(int responseID, Information &eleInformation);

  protected:
    
  private:
	double computeCurrentNormalStrain(void);
   
    // private attributes - a copy for each object of the class
    UniaxialMaterial *theMaterial;  // pointer to a material
    ConcretewBeta *theBetaMaterial; //consider making abstract class for materials that are sensitive to normal strain.
    ID  connectedExternalNodes;     // contains the tags of the end nodes
	ID  connectedExternalOtherNodes;     // contains the tags of the end nodes
    int numDOF;	                    // number of dof for CorotTruss2
    int numDIM;                     // number of dimensions

    double Lo;	        // initial length of truss
    double Ln;		    // current length of truss
    double d21[3];	    // current displacement offsets in basic system
    double v21[3];      // current velocity offsets in basic system
    double A; 	        // area of CorotTruss2
    double rho; 	    // mass density per unit length

	double otherLength;
	double otherLength_new;
	double theta; // for getting the normal strain
	double od21[3]; //displacement offsets for auxiliary nodes.

    Node *theNodes[2];
	Node *theOtherNodes[2];  // for computing et.

    Matrix R;	// Rotation matrix

    Matrix *theMatrix;

    static Matrix M2;
    static Matrix M4;
    static Matrix M6;
    static Matrix M12;
    
    Vector *theVector;

    static Vector V2;
    static Vector V4;
    static Vector V6;
    static Vector V12;
};

#endif




