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
                                                                        
#ifndef ZeroLengthImpact3D_h
#define ZeroLengthImpact3D_h

// Developed by: Prof. Arash E. Zaghi & Majid Cashany, University of Connecticut (UConn) Copyright 2012
//   Based on ZereLengthContact element by Gang Wang.
//
// What: "@(#) ZeroLengthImpact3D.h, revA"

#include <Element.h>
#include <Matrix.h>
#include <Vector.h>

// Tolerance for zero length of element
#define	LENTOL 1.0e-6

class Node;
class Channel;
class UniaxialMaterial;
class Response;

class ZeroLengthImpact3D : public Element
{
  public:
    // constructors
    ZeroLengthImpact3D(int tag, 
	    int Nd1, int Nd2, 
	    int direction, 
		double initGapInput, double fRatio, double Kt, 
		double Kn, double Kn2Input, double Delta_yInput, 
		double c);

    ZeroLengthImpact3D();    
    
    // destructor
    ~ZeroLengthImpact3D();

    
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
    //int update(void);

    // public methods to obtain stiffness
    const Matrix &getTangentStiff(void);
    const Matrix &getInitialStiff(void);
    const Matrix &getDamp(void);
    const Matrix &getMass(void);
    
    void zeroLoad(void);
    int addLoad(ElementalLoad *theLoad, double loadFactor);
    int addInertiaLoadToUnbalance(const Vector &accel);

    // public method to obtain resisting force
    const Vector &getResistingForce(void);
    const Vector &getResistingForceIncInertia(void);

    // public methods for output    
    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);
    int displaySelf(Renderer &, int mode, float fact, const char **displayModes=0, int numModes=0);
    void Print(OPS_Stream &s, int flag =0);    

    Response *setResponse(const char **argv, int argc, OPS_Stream &s);
    int getResponse(int responseID, Information &eleInformation);

  protected:
    
  private:
    // private member functions - only available to objects of the class
    int    directionID;
    ID     connectedExternalNodes;         // contains the tags of the end nodes
    static const int  numberNodes  ;
    Node *nodePointers[2];   // node pointer
    // contact forces
    double pressure;    // contact pressure
    double t1;          // friction in local dir1
    double t2;          // friction in local dir2
    // parameters
    double gap;         // gap of time n+1 step
    double gap_n;       // gap of time n step
    double Kn;			// normal penalty
    double Kt;			// tangential penalty
    double fs;			// friction ratio
    double cohesion;    // cohension
    // contact point and contact plane stuffs
    Vector stickPt;  // (keci_1, keci_2)
    Vector xi;       // trial stick point in local coord
    Vector origin;   // (x0,y0) for circular impact
    // Normal and Tangental Vectors for Elemental Nodes, (6*1)
    Vector N;
    Vector T1;
    Vector T2;
    int ContactFlag;                    // 0: not contact; 1: stick; 2: slide
    int numDOF;	                        // number of dof for ZeroLength
    // detect the contact and set flag
    int contactDetect();
    //form residual and tangent
    void formResidAndTangent(int tang_flag ) ;
    Matrix *Ki; 	    	// pointer to objects matrix (a class Matrix)
    Vector *load;         	// pointer to objects vector (a class Vector)
    //static variables for 3D contact
    static Matrix stiff;   // for stiff matrix
    static Vector resid;   // for force residual vector
    static Matrix zeroMatrix;

	double initGap; //initial gap
	double Kn1; //normal stiffness before yielding
	double Kn2; //normal stiffness after yielding
	double Delta_y; //yielding displacement

	void KnANDpressure(); // impact material implementation
	double pressT;
	double pressC;
	double gapC;
	double gapT;
	double gapD;
	double tangentT;

};
#endif

