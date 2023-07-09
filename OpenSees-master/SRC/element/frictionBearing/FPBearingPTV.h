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

// Written: Manish Kumar (mkumar3@buffalo.edu)
// Created: 2014/12/22
//
// Description: This file contains the implementation of the
// FPBearingPTV class. This class is similar to the existing
// singleFPBearing class. The key difference lies in 
// consideration of dependence of coefficient of friction on
// instantaneous values of sliding velocity, axial pressure on 
// the bearing, and temperature at the sliding surface.

#ifndef FPBearingPTV_h
#define FPBearingPTV_h
#include <Element.h>
#include <Matrix.h>
#include <Vector.h>

class Channel;
class UniaxialMaterial;
class Response;
class Channel;

class FPBearingPTV : public Element
{
public:
    // constructor
    FPBearingPTV(int tag, int Nd1, int Nd2, double MuReference,
		int IsPDependent, double refP, int IsTDependent, double Diffusivity_Steel,
		double Conductivity_Steel, int IsVDependent, double rate_v_mu, double Reff,
		double r_Contact, double kInit, 
		UniaxialMaterial &theMatA, UniaxialMaterial &theMatB,
	    UniaxialMaterial &theMatC, UniaxialMaterial &theMatD, 
        const Vector x = 0, const Vector y = 0,
        double shearDistI = 0.0,
        int addRayleigh = 0, double mass = 0.0,
        int maxIter = 25, double tol = 1E-12, int unit = 1);
    FPBearingPTV();

    
    // destructor
    ~FPBearingPTV();
    
    // method to get class type
    const char *getClassType() const {return "FPBearingPTV";};
    
    // public methods to obtain information about dof & connectivity
    int getNumExternalNodes() const;
    const ID &getExternalNodes();
    Node **getNodePtrs();
    int getNumDOF();
    void setDomain(Domain *theDomain);
    
    // public methods to set the state of the element
    int commitState();
    int revertToLastCommit();
    int revertToStart();
    int update();
    
    // public methods to obtain stiffness, mass, damping and residual information
    const Matrix &getTangentStiff();
    const Matrix &getInitialStiff();
    const Matrix &getDamp();
    const Matrix &getMass();
    
    void zeroLoad();
    int addLoad(ElementalLoad *theLoad, double loadFactor);
    int addInertiaLoadToUnbalance(const Vector &accel);
    
    const Vector &getResistingForce();
    const Vector &getResistingForceIncInertia();
    
    // public methods for element output
    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);
    int displaySelf(Renderer &theViewer, int displayMode, float fact, const char **modes, int numMode);
    void Print(OPS_Stream &s, int flag = 0);
    
    Response *setResponse(const char **argv, int argc, OPS_Stream &s);
    int getResponse(int responseID, Information &eleInformation);

	int setParameter(const char **argv, int argc, Parameter &param);
    //int updateParameter(int parameterID, Information &info);

	/*//AddingSensitivity:BEGIN////////////////////////////////
	int activateParameter(int parameterID);
	double getStressSenitivity(int gradIndex, bool conditional);
	double getInitialTangentSensitivity(int gradIndex);
	int commitSensitivity(double strainGradient, int gradIndex, int numGrads);
	//AddingSensitivity:END//////////////////////////////////*/
    
protected:

private:
    // private methods
    void setUp();
    double sgn(double x);
    
    // private attributes - a copy for each object of the class
    ID connectedExternalNodes;          // contains the tags of the end nodes
    Node *theNodes[2];                  // array of nodes   

	//Materials   
	UniaxialMaterial *theMaterials[4];  // uniaxial materials
    
    // parameters
    double muRef;		//Reference coefficient of friction
	int kpFactor;		//If friction is pressure dependent 1, 0 otherwise
	double refPressure; //Reference axial pressure
	int kTFactor;		//If friction is temperature dependent 1, 0 otherwise
	double diffuse;		//Thermal diffusivity of steel
	double conduct;		//Thermal conductivity of steel
	int kvFactor;		//If friction is velocity dependent 1, 0 otherwise
	double rateParam;	//rate parameter for velocity dependence
	int unit;			//Units of force,displacement etc. 1: N,m,s,C; 2: kN,m,s,C; 3: N,mm,s,C; 4: kN,mm,s,C; 5: lb,in,s,C; 6: kip,in,s,C; 7: lb,ft,s,C; 8: kip,ft,s,C
	double k0;          // initial stiffness of hysteretic component
    Vector x;           // local x direction
    Vector y;           // local y direction
    double shearDistI;  // shear distance from node I as fraction of length
    int addRayleigh;    // flag to add Rayleigh damping
    double mass;        // mass of element
    int maxIter;        // maximum number of iterations
    double tol;         // tolerance for convergence criterion
    double L;           // element length
	double Reffective; //Effective radius of curvature
	double rContact;	//Radius of the area of contact
    
    // state variables
    Vector ub;          // displacements in basic system
    Vector ubPlastic;   // plastic displacements in basic system
    Vector qb;          // forces in basic system
    Matrix kb;          // stiffness matrix in basic system
    Vector ul;          // displacements in local system
    Matrix Tgl;         // transformation matrix from global to local system
    Matrix Tlb;         // transformation matrix from local to basic system
	Vector DomainTime; //Current time in the domain	
	Vector DomainTimeTemp; //To be used to save data while resizing DomainTime
	Vector DomainHeatFlux; //History of heat flux
	Vector DomainHeatFluxTemp; //To be used to save data while resizing DomainHeatFlux	
	Vector kpFTemp;
	Vector kTFTemp;
	Vector kvFTemp;
	Vector TemperatureCenter;	
	Vector MuFactors;
	Vector MuAdjusted;
	Vector HeatFluxCenter;
    
    // committed history variables
    Vector ubPlasticC;  // plastic displacements in basic system
    
    // initial stiffness matrix in basic system
    Matrix kbInit;
	Matrix DomainDisp; //Displacement history	
    
    static Matrix theMatrix;
    static Vector theVector;
    static Vector theLoad;

	// new parameters
	double mu;	
	int iCountTime;

	//AddingSensitivity:Begin////////////////////////////////
	int parameterID;
	//AddingSensitivity:End//////////////////////////////////
};

#endif
