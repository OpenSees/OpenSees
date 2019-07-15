// Written by: Quan Gu_1,  Yongdou Liu_1, Wei Guo_23, Weiquan Li_1, Zhiwu Yu_23, Lizhong Jiang_23 and Hanyun Liu_2
//(1.School of Architecture and Civil Engineering, Xiamen University, 361005, China;
// 2.School of Civil Engineering, Central South University, 410075, China;
// 3.National Engineering Laboratory for High-speed Railway Construction, 410075, China)

// Reference: Quan Gu, Yongdou Liu, et al. A Practical Wheel-Rail Interaction Element for Modelling Vehicle-Track Systems[J]. 
// International Journal of Structural Stability and Dynamics
//
// Created: 09/2017
// revised : on 10/28/2018 by Yongdou Liu
// Copyright by the writers. 
                    
#ifndef WheelRail_h
#define WheelRail_h

// Description: This file contains the class definition for Wheel Rail Dynamic Interaction. 
// What: "@(#) WheelRail.h, revA"

#include <Element.h>
#include <ElasticBeam2d.h>
#include <Matrix.h>
#include <Vector.h>
class Channel;
class UniaxialMaterial;

#define ELE_TAG_WheelRail 10035443001

class WheelRail : public Element
{
  public:
    WheelRail(int pTag, double pDeltT, double pVel, double pInitLocation, int pNd1, 
	      double pRWheel,double pI,double pE,double pA,CrdTransf *ptheCoordTransf,
	      int pnLoad, Vector *pNodeList,
	      Vector * pDeltaYList=0,Vector * pDeltaYLocationList=0);
    ~WheelRail();

    const char *getClassType(void) const {return "WheelRail";};

    // public methods to obtain inforrmation about dof & connectivity    
    int getNumExternalNodes(void) const;
    const ID &getExternalNodes(void);
    Node **getNodePtrs(void);
    int getNumDOF(void);	
    double getResidualOfDeltaU(double pFhz,double uWheel);
    void setDomain(Domain *theDomain);
    void getDeltaY();

    // public methods to set the state of the element    
    int commitState(void);
    int revertToLastCommit(void);        
    int revertToStart(void);        
    int update(void);
    
    // public methods to obtain stiffness, mass, damping and 
    const Matrix &getTangentStiff(void);
    const Matrix &getInitialStiff(void);    
    //---------new Algorithm--------------
    void NewtonBisection(Vector limits,double uWheel);
    double FalsePostionAlgorithm(Vector limits,double uWheel);
    void getActiveDof();
    void getShapeFuns();

    void zeroLoad(void);	
    int addLoad(const Vector &addP);
    int addInertiaLoadToUnbalance(const Vector &accel);
    const Vector &getResistingForce(void);
    const Vector &getResistingForceIncInertia(void);            
    // public methods for element output
    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, 
		 FEM_ObjectBroker &theBroker);
    int displaySelf(Renderer &theViewer, int displayMode, float fact);    
    void Print(OPS_Stream &s, int flag =0);    
    
    Response *setResponse(const char **argv, int argc, OPS_Stream &);
    int getResponse(int responseID, Information &eleInformation);
    
    int setParameter (const char **argv, int argc, Parameter &param);
    int updateParameter (int parameterID, Information &info);
    
 protected:
    
 private:
//--------------------members in the construtor--------------------------

    ID  connectedExternalNodes;    // must reset using above four values   // contains the tags of the end nodes
    ID activeDof;
    
    Node **theNodes;
    Vector *P;
    Matrix *theTangent;
    
    CrdTransf *theCoordTransf;
    Domain *theDomain;
    
    Vector shapFun1,shapFun2;
    
    double deltT,vel,initLocation;
    double I,E,A;
    double currentLocation;  // from left node of the last commited beam element.
    double deltaU,Fhz,uF,theDeltaY;
    double theEleLength,a,b;
    double G, rollingRadiusWheel;
    
    int nLoad, loadStep, activeBeamIndex;
    int wheelNodeNum;
    int numRailNodeList;
    int theNumOfDeltaYList;
    
    Vector rearRailNode,frontRailNode,railDisp;
    
    Vector * theNodeList;
    Vector * theDeltaYList;
    Vector * theDeltaYLocationList;
    
    static Vector contactData,localActiveForce,activeData;


};

#endif




