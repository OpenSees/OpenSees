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

// Written by: Manish Kumar (mkumar2@buffalo.edu)
// Credits: This element extends the formulation of elastomericBearing element written by Andreas Schellenberg 
// Created: 02/29/2012
// Revision: A

#ifndef ElastomericX_h
#define ElastomericX_h

#include <Element.h>
#include <Matrix.h>
#include <Vector.h>

class Channel;
class UniaxialMaterial;
class Response;

class ElastomericX : public Element
{
public:
    // Constructor
    ElastomericX(int eleTag, int Nd1, int Nd2, double fy, double alpha, double Gr, double Kbulk,
		double D1, double D2, double ts, double tr, int n, const Vector y, const Vector x=0, 
		double kc=10, double PhiM=0.75, double ac=1.0, double sDratio=0.5, double m=0.0, 
		double cd=0.0, double tc=0.0, int tag1=0, int tag2=0, int tag3=0, int tag4=0);
        
    ElastomericX();
       
     // Destructor
    ~ElastomericX();
       
    // Method to get class type
    const char *getClassType() const {return "ElastomericX";};
   
    // Public methods to obtain information about dof & connectivity    
    int getNumExternalNodes() const;
    const ID &getExternalNodes();
    Node **getNodePtrs();
    int getNumDOF();
    void setDomain(Domain *theDomain);
       
    // Public methods to set the state of the element    
    int commitState();
    int revertToLastCommit();        
    int revertToStart();        
    int update();

    // Public methods to obtain stiffness, mass, damping and residual information    
    const Matrix &getTangentStiff();
    const Matrix &getInitialStiff();
	const Matrix &getDamp();
    const Matrix &getMass();
       
    void zeroLoad();
    int addLoad(ElementalLoad *theLoad, double loadFactor);
    int addInertiaLoadToUnbalance(const Vector &accel);
   
    const Vector &getResistingForce();
    const Vector &getResistingForceIncInertia();
   
    // Public methods for element output
    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);
    int displaySelf(Renderer &theViewer, int displayMode, float fact, const char **modes, int numMode);
    void Print(OPS_Stream &s, int flag = 0);    
       
    // Public methods for element recorder
    Response *setResponse(const char **argv, int argc, OPS_Stream &s);
    int getResponse(int responseID, Information &eleInfo);
   
protected:

private:
    // Private methods
    void setUp();
    double sgn(double x);
   
    // Private attributes - a copy for each object of the class
    ID connectedExternalNodes;					// contains the tags of the end nodes
    Node *theNodes[2];							// array of nodes
        UniaxialMaterial *theMaterials[3];		// array of uniaxial materials
   
    // PARAMETERS
	// Horizontal direction
    double k0;									// Initial stiffness of hysteretic component(due to lead)
	double qYield;								// Yield force of hysteretic component/chacteristic strength of bearing
    double ke;									// Stiffness of elastic component(due to rubber)
	double cd;									// Viscous damping parameter

	// Vertical direction
	double S;									// Shape facor
	double Ec;									// Compression modulus
	double Kv0;									// Stiffness at zero horizontal displacement
	double Kv;									// Elastic stiffness in compression and tension
	double kc;									// Quality index of elastomer
	double PhiM;								// Maximum reduction in the cavitation strength of elastomer
	double ac;									// Strength degradation parameter
	double Fcr;									// Critical buckling load at zero lateral deformation
	double ucr;									// Critical buckling deformation at zero lateral deformation
	double Fc;									// Initial cavitation strength
	double uc;									// Initial cavitation deformation
	double tCurrent, tCommit;

	double Kt, Kr;								// Torsional and rotational stiffness of bearing


	// Others
	double G;									// Shear modulus of rubber
	Vector x;									// local x direction
    Vector y;									// local y direction
	int tag1, tag2, tag3, tag4, tag5;			// Vector of tags to include or exclude a particular behavior
    double shearDistI;							// shear distance from node I as fraction of length
    double mass;								// mass of element
	double Tr;									// height of rubber in the bearing
	double D1, D2;								// Inner and outer diameter of the bearing
    double L;									// element length
	double h;									// height of rubber + shims
	double rg;									// radius of gyration of bearing
	double A;									// Bonded rubber area of bearing
	double Ar;									// reduced bonded area due to shear displacement
	double n, ts;								// number of layers and shim thickness

    // State variables
	double Fcrn;								// Current critical buckling force at particular shear deformation
	double ucrn;								// Current critical buckling displacement at particular shear deformation
	double Fcrmin;								// Min value of critical buckling load during whole loading	
	double Fcn, ucn;									// Current cavitation strength and deformation
	double Fmax, umax;								// Maximum force and deformation ever experienced by the elastomer
	//int tag;									// Tag to identify if bearing has failed in buckling

    Vector ub;									// Displacements in basic system
	Vector ubdot;								// Velocities in basic system
    Vector z;									// Hysteretic dimensionless quantity
	Matrix dzdu;								// Tangent of hysteretic evolution parameters
    Vector qb;									// Forces in basic system
    Matrix kb;									// Stiffness matrix in basic system
    Vector ul;									// Displacements in local system
    Matrix Tgl;									// Transformation matrix from global to local system
    Matrix Tlb;									// Transformation matrix from local to basic system
   
    // Committed history variables
	Vector ubC;									// Displacements in basic system
    Vector zC;									// Plastic displacements in basic system
    // Initial stiffness matrix in basic system
    Matrix kbInit;
   
    static Matrix theMatrix;
    static Vector theVector;
    static Vector theLoad;
};

#endif
