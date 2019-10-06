/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** ****************************************************************** */

#ifndef ElastomericBearingBoucWenMod3d_h
#define ElastomericBearingBoucWenMod3d_h

// Written: Ping Tan, Guangzhou Univ. China
// Created: 2015
// Revision: A
//
// Description: This element extends the formulation of elastomericBearing element
// written by Andreas Schellenberg. Parameters are introduced to consider compresssive
// stress dependency and temperature dependency of horizontal stiffness. Material properties 
// are defined and calculated within this element to simplify users' inputs.

#include <Element.h>
#include <Matrix.h>
#include <Vector.h>
class Channel;
class UniaxialMaterial;
class Response;

class ElastomericBearingBoucWenMod3d : public Element
{
public:
    // constructor
  // using geometry parameters input instead of material definition
  // ai and bi are introduced to consider compressive stress dependency 
  // and temperature dependency.
  ElastomericBearingBoucWenMod3d(int tag, int Nd1, int Nd2, 
				 double kInit, double fy, 
				 double Gr, double Kbulk, double D1, double D2,
				 double ts, double tr, int n,double alpha1,
				 double alpha2 = 0.0, double mu = 2.0,double eta = 1.0,
				 double beta = 0.5, double gamma = 0.5, 
				 double a1 = 0.0, double a2 = 1.0, double T = 23.0,
				 double b1 = 1.0, double b2 = 0.0, double b3 = 0.0, double b4 = 0.0,
				 const Vector y = 0, const Vector x = 0,
				 double shearDistI = 0.5,
				 int addRayleigh = 0, double mass = 0.0,
				 int maxIter = 25, double tol = 1E-12);
  ElastomericBearingBoucWenMod3d();
  
  // destructor
  ~ElastomericBearingBoucWenMod3d();
  
  // method to get class type
  const char *getClassType() const {return "ElastomericBearingBoucWenMod3d";};
  
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
  int displaySelf(Renderer &theViewer, int displayMode, float fact);
  void Print(OPS_Stream &s, int flag = 0);
  
  // public methods for element recorder
  Response *setResponse(const char **argv, int argc, OPS_Stream &s);
  int getResponse(int responseID, Information &eleInfo);
  
protected:
  
 private:
  // private methods
  void setUp();
  double sgn(double x);
  
  // private attributes - a copy for each object of the class
  ID connectedExternalNodes;          // contains the tags of the end nodes
  Node *theNodes[2];                  // array of nodes
  
  
  // parameters
  double G;			// Shear modulus of rubber
  double k0;          // initial stiffness of hysteretic component
  double qYield;      // Yield force of hysteretic component/chacteristic strength of bearing
  double kInit;       // initial stiffness 
  double fy;          // yield strength
  double alpha1;      //post yield stiffness ratio of linear hardening component 
  double alpha2;      //post yield stiffness ratio of non-linear hardening component
  double k2;          // stiffness of linear elastic component
  double k3;          // stiffness of nonlinear elastic component
  double a1, a2;      // parameters of compressive stress dependency 
  double T;           // temperature
  double b1, b2, b3, b4;  //parameters of temperature dependency
  double mu;          // exponent of nonlinear elastic component
  double eta;         // yielding exponent (sharpness of hysteresis loop corners)
  double beta;        // hysteretic shape parameter
  double gamma;       // hysteretic shape parameter
  double A;           // restoring force amplitude parameter
  Vector x;           // local x direction
  Vector y;           // local y direction
  double shearDistI;  // shear distance from node I as fraction of length
  int addRayleigh;    // flag to add Rayleigh damping
  double mass;        // mass of element
  int maxIter;        // maximum number of iterations
  double tol;         // tolerance for convergence criterion
  double L;           // element length
  bool onP0;          // flag to indicate if the element is on P0
  double D1, D2;		// Inner and outer diameter of the bearing
  double Tr;			// height of rubber in the bearing
  double S;			// Shape facor
  double n, ts;		// number of layers and shim thickness
  double h;			// height of rubber + shims
  double Ec;			// Compression modulus
  double Kv0;			// Stiffness at zero horizontal displacement
  double Kt, Kr;		// Torsional and rotational stiffness of bearing
  
  
  // state variables
  Vector ub;          // displacements in basic system
  Vector z;           // hysteretic evolution parameters
  Matrix dzdu;        // tangent of hysteretic evolution parameters
  Vector qb;          // forces in basic system
  Matrix kb;          // stiffness matrix in basic system
  Vector ul;          // displacements in local system
  Matrix Tgl;         // transformation matrix from global to local system
  Matrix Tlb;         // transformation matrix from local to basic system
  
  // committed history variables
  Vector ubC;         // displacements in basic system
  Vector zC;          // hysteretic evolution parameters
  
  // initial stiffness matrix in basic system
  Matrix kbInit;
  
  static Matrix theMatrix;  // a class wide Matrix
  static Vector theVector;  // a class wide Vector
  Vector theLoad;
};

#endif
