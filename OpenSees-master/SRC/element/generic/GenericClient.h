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

#ifndef GenericClient_h
#define GenericClient_h

// Written: Andreas Schellenberg (andreas.schellenberg@gmail.com)
// Created: 11/06
// Revision: A
//
// Description: This file contains the class definition for GenericClient.
// GenericClient is a generic element defined by any number of nodes and 
// the degrees of freedom at those nodes. The element communicates with 
// OpenFresco through a tcp/ip connection.

#include <Element.h>
#include <Matrix.h>

#define RemoteTest_open              1
#define RemoteTest_setup             2
#define RemoteTest_setTrialResponse  3
#define RemoteTest_execute           4
#define RemoteTest_commitState       5
#define RemoteTest_getDaqResponse    6
#define RemoteTest_getDisp           7
#define RemoteTest_getVel            8
#define RemoteTest_getAccel          9
#define RemoteTest_getForce         10
#define RemoteTest_getTime          11
#define RemoteTest_getInitialStiff  12
#define RemoteTest_getTangentStiff  13
#define RemoteTest_getDamp          14
#define RemoteTest_getMass          15
#define RemoteTest_DIE              99

class Channel;


class GenericClient : public Element
{
public:
    // constructors
    GenericClient(int tag, ID nodes, ID *dof,
          int port, char *machineInetAddr = 0,
          int ssl = 0, int udp = 0, int dataSize = 256,
          int addRayleigh = 1);
    GenericClient();
    
    // destructor
    ~GenericClient();
    
    // method to get class type
    const char *getClassType() const {return "GenericClient";};
    
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
    
    // public methods to obtain other response in global system
    //const Vector &getTime();
    
    // public methods to obtain other response in basic system
    //const Vector &getBasicDisp();
    //const Vector &getBasicVel();
    //const Vector &getBasicAccel();
    
    // public methods for element output
    int sendSelf(int commitTag, Channel &sChannel);
    int recvSelf(int commitTag, Channel &rChannel, FEM_ObjectBroker &theBroker);
    int displaySelf(Renderer &theViewer, int displayMode, float fact, const char **modes, int numMode);
    void Print(OPS_Stream &s, int flag = 0);
    
    // public methods for element recorder
    Response *setResponse(const char **argv, int argc, OPS_Stream &s);
    int getResponse(int responseID, Information &eleInfo);
    
protected:
    
private:
    // private attributes - a copy for each object of the class
    ID connectedExternalNodes;  // contains the tags of the end nodes
    ID *theDOF;                 // array with the dof of the end nodes
    ID basicDOF;                // contains the basic dof
    
    int numExternalNodes;       // number of external nodes
    int numDOF;                 // number of total DOF
    int numBasicDOF;            // number of used DOF
    
    int port;                   // ipPort
    char *machineInetAddr;      // ipAddress
    int ssl;                    // secure socket layer flag
    int udp;                    // udp socket flag
    int dataSize;               // data size of send/recv vectors
    int addRayleigh;            // flag to add Rayleigh damping
    
    Matrix theMatrix;           // objects matrix
    Vector theVector;           // objects vector
    Vector theLoad;             // load vector
    Matrix theInitStiff;        // initial stiffness matrix
    Matrix theMass;             // mass matrix
    
    Channel *theChannel;        // channel
    double *sData;              // send data array
    Vector *sendData;           // send vector
    double *rData;              // receive data array
    Vector *recvData;           // receive vector
    
    Vector *db;                 // trial displacements in basic system
    Vector *vb;                 // trial velocities in basic system
    Vector *ab;                 // trial accelerations in basic system
    Vector *t;                  // trial time
    
    Vector *qDaq;               // daq forces in basic system
    Matrix *rMatrix;            // receive matrix
    
    Vector dbCtrl;              // ctrl displacements in basic system
    Vector vbCtrl;              // ctrl velocities in basic system
    Vector abCtrl;              // ctrl accelerations in basic system
    
    bool initStiffFlag;
    bool massFlag;
    
    Node **theNodes;
    
    int setupConnection();
};

#endif
