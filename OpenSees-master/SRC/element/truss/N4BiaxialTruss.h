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


#ifndef N4BiaxialTruss_h
#define N4BiaxialTruss_h

#include <Element.h>
#include <Matrix.h>
#include <ConcretewBeta.h>

class Node;
class Channel;
class UniaxialMaterial;

class N4BiaxialTruss : public Element
{
  public:
    N4BiaxialTruss(int tag, 
	  int dimension,
	  int Nd1, int Nd2, 
	  int GNd1, int GNd2, 
	  UniaxialMaterial &theMaterial,
	  double A, 
	  double rho=0.0, 
	  int doRayleighDamping = 0);
    
    N4BiaxialTruss();    
    ~N4BiaxialTruss();

    const char *getClassType(void) const {return "N4BiaxialTruss";};

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
    const Matrix &getKi(void);
    const Matrix &getTangentStiff(void);
    const Matrix &getInitialStiff(void);
    const Matrix &getDamp(void);    
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
    int computeCurrentStrainBiaxial(void);
	int computeCurrentStrainRate(void);
    
    // private attributes - a copy for each object of the class
    UniaxialMaterial *theMaterial_1;  // pointer to a material for the 1 direction
	ConcretewBeta *theBetaMaterial_1; // pointer for biaxial concrete material (if available - 0 otherwise)
    UniaxialMaterial *theMaterial_2;  // pointer to a material for the 1 direction
	ConcretewBeta *theBetaMaterial_2; // pointer for biaxial concrete material (if available - 0 otherwise)
    ID  connectedExternalNodes;     // contains the tags of the end nodes
    int dimension;                  // truss in 2 or 3d domain
    int numDOF;	                    // number of dof for truss

    Vector *theLoad;     // pointer to the load vector P
    Matrix *theMatrix; // pointer to objects matrix (a class wide Matrix)
    Vector *theVector; // pointer to objects vector (a clas wide Vector)
    Vector *theVector2; // pointer to objects half size vector (a clas wide Vector)

    double L;	    // length of truss based on undeformed configuration
    double L2;	    // length of truss based on undeformed configuration
    double A; 	    // area of truss
    double rho; 	// rho: mass density per unit length
    int doRayleighDamping; 

    double vectorX[3]; // direction cosines
    double vectorY[3]; // direction cosines
	double lengthX; // length of rectangle in relative X coord
	double lengthY; // length of rectangle in relative Y coord
	
    Node *theNodes[4];
    double cosX[3]; // direction cosines
    double cosX2[3]; // direction cosines for second direction
	
	// Output variables:
	double strain_1;
	double strain_2;
	double normalStrain_1;
	double normalStrain_2;
	double strainRate_1;
	double strainRate_2;
	
	// Computed variables
	double oneOverL;
	double LxoverL;
	double LyoverL;
	double oneOver2Lx;
	double oneOver2Ly;

    // static data - single copy for all objects of the class	
    static Matrix trussM2;   // class wide matrix for 2*2
    static Matrix trussM8;   // class wide matrix for 8*8
    static Matrix trussM12;   // class wide matrix for 12*12
    static Matrix trussM24;  // class wide matrix for 24*24
    static Vector trussV2;   // class wide Vector for size 2
    static Vector trussV4;   // class wide Vector for size 4
    static Vector trussV6;   // class wide Vector for size 6
    static Vector trussV8;   // class wide Vector for size 8
    static Vector trussV12;   // class wide Vector for size 12
    static Vector trussV24;  // class wide Vector for size 24
};

#endif




