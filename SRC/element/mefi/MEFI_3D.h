// Code written/implemented by:	Carlos Lopez Olea (carlos.lopez.o@ug.uchile.cl)
//
// User documentation available at: https://github.com/carloslopezolea/MEFI
//
// Created: 01/2026
//
// Description: The three-dimensional Membrane Fiber element (MEFI_3D) is a four-node element with six degrees of freedom (DOFs) per node: three translational DOFs and three rotational DOFs.
// The in-plane response is based on the MEFI_2D formulation, whereas the out-of-plane response follows the Kirchhoff plate formulation with four integration points per element. Both behaviors are formulated independently, 
// providing an uncoupled representation of membrane and bending actions in reinforced concrete walls.
//
// Reference:
// 1.- Lopez, C. N., Rojas, F., & Massone, L. M. (2022). Membrane fiber element for reinforced concrete walls - the benefits of macro and micro modeling approaches. Engineering Structures, 254, 113819.
// 2.- Suquillo, B., Rojas, F., López, C. et al. MEFI-3D: A membrane fiber element for non-planar reinforced concrete structural walls. Bull Earthquake Eng 24, 211–238 (2026).
//
// Source: /usr/local/cvs/OpenSees/SRC/element/MEFI/MEFI_3D.cpp
//
// Rev: 1.0          

#ifndef MEFI_3D_h
#define MEFI_3D_h

#ifndef _bool_h
#include "bool.h"
#endif

#include <Element.h>
#include <Matrix.h>
#include <Vector.h>
#include <ID.h>

class Node;
class SectionForceDeformation;
class Response;

class MEFI_3D : public Element
{
  public:
      MEFI_3D(int tag,                              // element tag
          int nd1, int nd2, int nd3, int nd4,       // node tags
          int numFibers,                            // number of fiber elements
          SectionForceDeformation** Sections,		// array of material tags
          double* Width,						// array of macro-fiber widths;
          double thickMod,
          double nu);

          
      MEFI_3D();
      ~MEFI_3D();

    const char *getClassType(void) const {return "MEFI_3D";};

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

    void zeroLoad();
    int addLoad(ElementalLoad *theLoad, double loadFactor);
    int addInertiaLoadToUnbalance(const Vector &accel);

    const Vector &getResistingForce(void);
    const Vector &getResistingForceIncInertia(void);            

    // public methods for element output
    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker  &theBroker);
    int displaySelf(Renderer &, int mode, float fact, const char **displayModes=0, int numModes=0);

    void Print(OPS_Stream &s, int flag =0);

    Response *setResponse(const char **argv, int argc, OPS_Stream &s);
    int getResponse(int responseID, Information &eleInformation);

  protected:
    
  private:

    // private static attributes
    static Matrix K;		        // Element stiffness, damping, and mass Matrix(24,24) in global CS
    static Matrix KL;		        // Element stiffness, damping, and mass Matrix(24,24) in local CS
    static Matrix K12;		        // Element stiffness, damping, and mass Matrix(12,12) in local CS
    static Matrix K24;		        // Element stiffness, damping, and mass Matrix(12,12) in local CS
    static Vector P;		        // Element resisting force, and displacement vector(24) in global CS
    static Vector PL;		        // Element resisting force, and displacement vector(24) in local CS
    static Vector P12;		        // Element resisting force, and displacement vector(12) in local CS
    static double detJac;	        // Element resisting force vector
    static Matrix BSD;		        // Element resisting force vector
    static Matrix DB;		        // Element resisting force vector

    static int iMemb[12];		        // Vector that relates membrane and element local DOFs 
    static int iBend[12];		        // Vector that relates bending and element local DOFs 

    static double qdtLocationsB[4][2];
    static double qdtWeightsB[4];
    

    // private attributes - a copy for each object of the class
    Node* theNodes[4];                      // external node pointers 
    SectionForceDeformation** theSection;   // pointer to section objects    
    ID connectedExternalNodes;              // Tags of nodes

    Vector Q;		        // Applied nodal loads
    
    // private member functions - only objects of this class can call these
    void membraneFieldInterpolation(double xi, double eta);
    void plateFieldInterpolation(double xi, double eta);
    const Matrix& transpose(const Matrix& M);
    void setTransformationMatrix(void);
    void setPlateTangent();

    // parameters for bending behavior
    double Eave;						
    double Tave;						
    double Nu;
    double tMod;						
    Matrix plateTangentNDM;


    // Nodal global coordinates
    Vector nd1Crds;
    Vector nd2Crds;
    Vector nd3Crds;
    Vector nd4Crds;

    // Nodal local coordinates
    Vector nd1CrdsL;
    Vector nd2CrdsL;
    Vector nd3CrdsL;
    Vector nd4CrdsL;

    // calculated element parameters
    int nip;					    // number of fibers
    Matrix qdtLocationsM;
    Vector qdtWeightsM;
    double lw;

    Matrix* Ki;
    Matrix T;							
    Matrix T6;							
    Matrix Tt;							
    Matrix KP;						
    Matrix BM;							
    Vector detJacM;
    Matrix BP;							
    Vector detJacP;
    Vector dispB;					    

};

#endif

