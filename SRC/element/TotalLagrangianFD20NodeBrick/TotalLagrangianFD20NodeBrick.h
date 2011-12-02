//===============================================================================
//# COPYRIGHT (C): Woody's license (by BJ):
//                 ``This    source  code is Copyrighted in
//                 U.S.,  for  an  indefinite  period,  and anybody
//                 caught  using it without our permission, will be
//                 mighty good friends of ourn, cause we don't give
//                 a  darn.  Hack it. Compile it. Debug it. Run it.
//                 Yodel  it.  Enjoy it. We wrote it, that's all we
//                 wanted to do.''
//
//# PROJECT:           Object Oriented Finite Element Program
//# PURPOSE:           Finite Deformation Hyper-Elastic classes
//# CLASS:
//#
//# VERSION:           0.6_(1803398874989) (golden section)
//# LANGUAGE:          C++
//# TARGET OS:         all...
//# DESIGN:            Zhao Cheng, Boris Jeremic (jeremic@ucdavis.edu)
//# PROGRAMMER(S):     Zhao Cheng, Boris Jeremic
//#
//#
//# DATE:              Sept2003
//# UPDATE HISTORY:    28May2004, Zhao Optimized the Stiffness Tensor
//#
//#
//===============================================================================
#ifndef TOTALLAGRANGIANFD20BRICK_H
#define TOTALLAGRANGIANFD20BRICK_H

#ifndef _bool_h
#include "bool.h"
#endif

#include <ID.h>
#include <Renderer.h>
#include <Domain.h>
#include <Element.h>
#include <Node.h>

#include <string.h>

#include <OPS_Globals.h>

#include <basics.h>
#include <nDarray.h>
#include <Vector.h>
#include <Matrix.h>
#include <BJtensor.h>
#include <stresst.h>
#include <straint.h>

#include <NDMaterial.h>

#include <Information.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <ElementResponse.h>
#include <ElementalLoad.h>


class TotalLagrangianFD20NodeBrick: public Element
{
  public:
    TotalLagrangianFD20NodeBrick(int tag,
    int node_numb_1,  int node_numb_2,  int node_numb_3,  int node_numb_4,
    int node_numb_5,  int node_numb_6,  int node_numb_7,  int node_numb_8,
    int node_numb_9,  int node_numb_10, int node_numb_11, int node_numb_12,
    int node_numb_13, int node_numb_14, int node_numb_15, int node_numb_16,
    int node_numb_17, int node_numb_18, int node_numb_19, int node_numb_20,
    NDMaterial &m, double b1=0.0, double b2=0.0, double b3=0.0);

    TotalLagrangianFD20NodeBrick ();
    ~TotalLagrangianFD20NodeBrick();
    
    const char *getClassType(void) const {return "TotalLagrangianFD20NodeBrick";};
    
    int getNumExternalNodes () const;
    const ID &getExternalNodes ();
    Node **getNodePtrs();

    int getNumDOF ();
    void setDomain(Domain *theDomain);

    int commitState ();
    int revertToLastCommit ();
    int revertToStart ();
    int update();

    const Matrix &getTangentStiff ();
    const Matrix &getInitialStiff();
    const Matrix &getMass ();

    void zeroLoad ();
    int addLoad(ElementalLoad *theLoad, double loadFactor);
    int addInertiaLoadToUnbalance(const Vector &accel);

    const Vector &getResistingForce ();
    const Vector &getResistingForceIncInertia ();

    // public methods for element output
    int sendSelf (int commitTag, Channel &theChannel);
    int recvSelf (int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);
    int displaySelf (Renderer &theViewer, int displayMode, float fact);
    void Print(OPS_Stream &s, int flag =0);

    Response *setResponse (const char **argv, int argc, Information &eleInformation, OPS_Stream &s);
    int getResponse (int responseID, Information &eleInformation);

//    int setParameter(const char **argv, int argc, Information &info);
//    int updateParameter(int parameterID, Information &info);


  protected:

  private:

    NDMaterial **theMaterial; // Pointer to the NDMaterial objects
    ID  connectedExternalNodes; // Tags of TotalLagrangianFD20Brick nodes
    Node *theNodes[20];

    static Matrix K;    // Element stiffness Matrix
//    static Matrix C;    // Element damping matrix
    static Matrix M;    // Element mass matrix
    static Vector P;    // Element resisting force vector
    static const double pts[3];   // Stores quadrature points
    static const double wts[3];   // Stores quadrature weights
    Vector *Q;     // Applied nodal loads
    Vector bf;    // Body forces

    double rho;    // Mass per unit volume

    double det_of_Jacobian;
    
    Matrix *Ki;

  private:
    
    static const int  NumIntegrationPts;
    static const int  NumTotalGaussPts;
    static const int  NumNodes;
    static const int  NumDof;
    static const int  NumElemDof;
    
    tensor shapeFunction(double , double , double );
    tensor shapeFunctionDerivative(double , double , double );

    tensor Jacobian_3D(double , double , double);
    tensor Jacobian_3Dinv(double , double , double);
    tensor dh_Global(double , double , double);
    tensor getNodesCrds(void);
    tensor getNodesDisp(void);

    tensor getStiffnessTensor(void);
    tensor getRtensor(void);
    tensor getBodyForce(void);
    tensor getSurfaceForce(void);
    tensor getForces(void);

};


#endif

