// Code written/implemented by:	Carlos Lopez Olea (carlos.lopez.o@ug.uchile.cl)
//
// User documentation available at: https://github.com/carloslopezolea/MEFI
//
// Created: 06/2022
//
// Description: The Membrane Fiber (MEFI) element, is described by four nodes, each containing three degrees of freedom (DOFs), two translations, and one in-plane rotation 
// (drilling) DOF, which incorporates a blended interpolation function for the displacements over the element (Fig. 1b-c). The element formulation accommodates the quadrature 
// points and weights of the classical finite element formulation of membrane elements to resemble strips (fibers), similarly to macroscopic elements.
//
// Reference:
// 1.- López, C. N., Rojas, F., & Massone, L. M. (2022). Membrane fiber element for reinforced concrete walls – the benefits of macro and micro modeling approaches. Engineering Structures, 254, 113819.
//
// Source: /usr/local/cvs/OpenSees/SRC/element/mefi/MEFI.h
//
// Rev: 1.0         

#ifndef MEFI_h
#define MEFI_h

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

class MEFI : public Element
{
  public:
      MEFI(int tag,                                 // element tag
          int nd1, int nd2, int nd3, int nd4,       // node tags
          int numFibers,                            // number of fiber elements
          SectionForceDeformation** Sections,		// array of material tags
          double* Width);							// array of macro-fiber widths;
          
      MEFI();
      ~MEFI();

    const char *getClassType(void) const {return "FourNodeQuad";};

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
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker 
		  &theBroker);

    int displaySelf(Renderer &, int mode, float fact, const char **displayModes=0, int numModes=0);
    void Print(OPS_Stream &s, int flag =0);

    Response *setResponse(const char **argv, int argc, 
			  OPS_Stream &s);

    int getResponse(int responseID, Information &eleInformation);

  protected:
    
  private:
    //private attributes - a copy for each object of the class
    Node* theNodes[4];                       // external node pointers 
    SectionForceDeformation ** theSection;   // pointer to the ND material objects   
    ID connectedExternalNodes;               // Tags of nodes
 
    //static attributes
    static Matrix K;		// Element stiffness, damping, and mass Matrix
    static Vector P;		// Element resisting force vector
    static double detJac;	// Element resisting force vector
    static Matrix BSD;		// Element resisting force vector
    Vector Q;		        // Applied nodal loads

    //private member functions - only objects of this class can call these
    void membraneFieldInterpolation(double xi, double eta);

    //nodal coordinates
    Vector nd1Crds;
    Vector nd2Crds;
    Vector nd3Crds;
    Vector nd4Crds;

    //calculated element parameters
    int nip;					   
    Matrix qdtLocations;
    Vector qdtWeights;
    double lw;

    const Matrix& transpose(const Matrix& M);

    Matrix* Ki;

};

#endif

