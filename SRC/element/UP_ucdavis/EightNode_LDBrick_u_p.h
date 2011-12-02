///////////////////////////////////////////////////////////////////////////////
//
// COPYRIGHT (C):     :-))
// PROJECT:           Object Oriented Finite Element Program
// FILE:              EightNodeBrick_u_p.cpp
// CLASS:             EightNodeBrick_u_p
// VERSION:
// LANGUAGE:          C++
// TARGET OS:         DOS || UNIX || . . .
// DESIGNER:          Zhao Cheng, Boris Jeremic
// PROGRAMMER:        Zhao Cheng, Boris Jeremic
// DATE:              Aug. 2006
// UPDATE HISTORY:    
//
///////////////////////////////////////////////////////////////////////////////

#ifndef EIGHTNODE_LDBRICK_U_P_H
#define EIGHTNODE_LDBRICK_U_P_H

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

class EightNode_LDBrick_u_p: public Element
{
public:
   EightNode_LDBrick_u_p(int element_number,
                   int node_numb_1, int node_numb_2, int node_numb_3, int node_numb_4,
                   int node_numb_5, int node_numb_6, int node_numb_7, int node_numb_8,
                   NDMaterial *Globalmmodel, double b1, double b2, double b3, 
                   double nf, double rs, double rf, 
                   double permb_x, double permb_y, double permb_z, 
                   double kkf); 
   EightNode_LDBrick_u_p();
   ~EightNode_LDBrick_u_p();
    
    const char *getClassType(void) const {return "EightNode_LDBrick_u_p";};   
    
    int getNumExternalNodes(void) const;
    const ID &getExternalNodes(void);
    Node **getNodePtrs(void);
    int getNumDOF(void);	
    void setDomain(Domain *theDomain);
   
    int commitState(void);
    int revertToLastCommit(void);        
    int revertToStart(void);        
    int update(void);
       
    const Matrix &getTangentStiff(void);
    const Matrix &getInitialStiff(void);    
    const Matrix &getDamp(void);
    const Matrix &getMass(void);

    void zeroLoad(void);	
    int addLoad(ElementalLoad *theLoad, double loadFactor);
    int addInertiaLoadToUnbalance(const Vector &accel);
    const Vector &getResistingForce(void);
    const Vector &getResistingForceIncInertia(void);            

    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);
    int displaySelf(Renderer &theViewer, int displayMode, float fact);    
    void Print(OPS_Stream &s, int flag =0);    

    Response *setResponse(const char **argv, int argc, Information &eleInfo, OPS_Stream &s);
    int getResponse(int responseID, Information &eleInformation);

    //int setParameter (const char **argv, int argc, Information &info);
    //int updateParameter (int parameterID, Information &info);

private:
    tensor shapeFunction(double, double, double);
    tensor shapeFunctionDerivative(double, double, double);
    tensor getNodesCrds();
    tensor getNodesDisp();
    tensor Jacobian_3D(double, double, double);
    tensor Jacobian_3Dinv(double, double, double);
    tensor dh_Global(double, double, double);

    tensor getGaussPts(void);
    double getPorePressure(double, double, double);

    tensor getStiffnessTensorK1();
    tensor getStiffnessTensorK2();
    tensor getMassTensorM1();
    tensor getMassTensorM2();
    tensor getDampingTensorC1();
    tensor getDampingTensorC2();
    tensor getMatStiffness();
    const Matrix& getStiff(int Ki_flag);
    const Matrix& getStiffnessK0();
    const Vector& getInternalForce();
    const Vector& getForceU();
    const Vector& getForceP();
    tensor LagrangianPerm(const tensor& Finv, const tensor& permea, double Jin);
        
private:
    ID  connectedExternalNodes;
    Node *theNodes[8]; 
    NDMaterial **theMaterial;

    static Matrix MCK;
    static Vector P;

    static const int  Num_IntegrationPts;
    static const int  Num_TotalGaussPts;
    static const int  Num_Nodes;
    static const int  Num_Dim;
    static const int  Num_Dof;
    static const int  Num_ElemDof;
    static const double pts[2];
    static const double wts[2]; 
    static tensor perm;            // Permeability = k/(rho_f*g)

    Vector bf;
    double nf;
    double rho_s;
    double rho_f;
    double kf;

    Vector *Q;		           
    Matrix *Ki;                    
};


#endif

