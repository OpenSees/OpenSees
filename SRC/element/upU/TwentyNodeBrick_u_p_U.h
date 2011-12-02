///////////////////////////////////////////////////////////////////////////////
//
// COPYRIGHT (C):     :-))
// PROJECT:           Object Oriented Finite Element Program
// FILE:             TwentyNodeBrick_u_p_U.h
// CLASS:            TwentyNodeBrick_u_p_U
// MEMBER FUNCTIONS:
//
// MEMBER VARIABLES
//
// PURPOSE:           Finite Element Class for coupled system
//  "Coupled system": Solid and fluid coexist.
//                    u-- Solid displacement
//                    p-- Pore pressure
//                    U-- Absolute fluid displacement
// RETURN:
// VERSION:
// LANGUAGE:          C++.ver >= 3.0
// TARGET OS:         DOS || UNIX || . . .
// DESIGNER:          Boris Jeremic, Xiaoyan Wu
// PROGRAMMER:        Boris Jeremic, Xiaoyan Wu
// DATE:              Aug. 2001
// UPDATE HISTORY:    Modified from EightNodeBrick_u_p_U.h.   Aug. 2001           
//                    01/16/2002    Xiaoyan
//                    Add the permeability tensor and ks, kf  to the constructor  Xiaoyan 
//
//                    Clean-up and re-write by Zhao Cheng, 10/20/2004
//
///////////////////////////////////////////////////////////////////////////////

#ifndef TWENTYNODEBRICK_U_P_U_H
#define TWENTYNODEBRICK_U_P_U_H

#ifndef _bool_h
#include "bool.h"
#endif

#include <Information.h>
#include <ElementResponse.h>
#include <ElementalLoad.h>
#include <Domain.h>
#include <Node.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Renderer.h>
#include <OPS_Globals.h>

#include <Element.h>
#include <NDMaterial.h>

#include <Vector.h>
#include <Matrix.h>
#include <BJtensor.h>
#include <stresst.h>
#include <straint.h>

class TwentyNodeBrick_u_p_U: public Element
{
public:
   TwentyNodeBrick_u_p_U(int element_number,
                   int node_numb_1, int node_numb_2, int node_numb_3, int node_numb_4,
                   int node_numb_5, int node_numb_6, int node_numb_7, int node_numb_8,
                   int node_numb_9, int node_numb_10, int node_numb_11, int node_numb_12,
                   int node_numb_13, int node_numb_14, int node_numb_15, int node_numb_16,
                   int node_numb_17, int node_numb_18, int node_numb_19, int node_numb_20,
		   NDMaterial *Globalmmodel, double b1, double b2, double b3,
		   double nn, double alf, double rs,double rf, 
		   double permb_x, double permb_y, double permb_z, 
		   double kks, double kkf, double pp); 
   TwentyNodeBrick_u_p_U ();
   ~TwentyNodeBrick_u_p_U();
    
    // public methods to obtain inforrmation about dof & connectivity    
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
    int displaySelf(Renderer &theViewer, int displayMode, float fact);    
    void Print(OPS_Stream &s, int flag =0);    

    Response *setResponse(const char **argv, int argc, Information &eleInfo);
    int getResponse(int responseID, Information &eleInformation);

    //int setParameter (const char **argv, int argc, Information &info);
    //int updateParameter (int parameterID, Information &info);

private:
    tensor shapeFunction(double, double, double);
    tensor shapeFunctionDerivative(double, double, double);
    tensor getGaussPts(void);
    tensor getNodesCrds(void);
    tensor getNodesDisp(void);
    tensor Jacobian_3D(tensor dh);
    tensor Jacobian_3Dinv(tensor dh);
    tensor dh_Global(tensor dh);

    tensor getStiffnessTensorKep();
    tensor getStiffnessTensorG12();
    tensor getStiffnessTensorP();
    tensor getMassTensorMsf();
    tensor getDampTensorC123();
    const Matrix &getStiff(int Ki_flag);
    double getPorePressure(double, double, double);
    Vector getInForceS();
    Vector getInForceP();
    Vector getExForceS();
    Vector getExForceF();
        
private:
    ID  connectedExternalNodes;    // tags of nodes
    Node *theNodes[20] ;           // pointers to eight nodes
    NDMaterial **theMaterial;      // pointers to the ND material objects

    static Matrix K;               // Stiffness
    static Matrix C;		   // Damping
    static Matrix M;		   // Mass
    static Vector P;		   // Residual

    static const int  Num_IntegrationPts;
    static const int  Num_TotalGaussPts;
    static const int  Num_Nodes;
    static const int  Num_Dim;
    static const int  Num_Dof;
    static const int  Num_ElemDof;
    static const double pts[3];    // Stores quadrature points
    static const double wts[3];    // Stores quadrature weights
    static tensor perm; 	   // Permeability = k/(rho_f*g)

    Vector bf;  	           // Body forces 
    double poro;                   // Porosity
    double alpha;		   // Coefficient for soil (approximate equal 1)
    double rho_s;      		   // Solid density
    double rho_f;      		   // Fluid density
    double ks;                     // Bulk modulus of solid
    double kf;                     // Bulk modulus of fluid
    double pressure;   		   // Normal surface traction (pressure) over entire element  //?
    
    Vector *Q;
    Vector *eleQ;		           
    Matrix *Ki;                    
};


#endif

