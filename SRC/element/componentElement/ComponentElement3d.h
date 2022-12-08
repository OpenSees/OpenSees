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
                                                                        
// Written: fmk 09/15, MHS 11/22 for 3d
// Revised:
//
// Purpose: This file contains the class definition for ComponentElement3d.
// ComponentElement3d is an elastic element with hinges at the ends. It uses
// static condensation to deal with the internal dof

#ifndef ComponentElement3d_h
#define ComponentElement3d_h

#include <Element.h>
#include <Node.h>
#include <Matrix.h>
#include <Vector.h>

class Channel;
class Information;
class CrdTransf;
class Response;
class Renderer;
class UniaxialMaterial;

class ComponentElement3d : public Element
{
  public:
    ComponentElement3d();        
    ComponentElement3d(int tag, double A, double E, double Iz,
		       double Iy, double G, double J,
		       int Nd1, int Nd2, CrdTransf &theTransf, 
		       UniaxialMaterial *end1z, UniaxialMaterial *end2z,
		       UniaxialMaterial *end1y, UniaxialMaterial *end2y,		       
		       double rho = 0.0, int cMass = 0);
    ComponentElement3d(int tag, double A, double E, double Iz,
		       double Iy, double G, double J,
		       int Nd1, int Nd2, CrdTransf &theTransf, 
		       double kzI, double kzJ, double kyI, double kyJ,
		       double rho = 0.0, int cMass = 0);  
    ~ComponentElement3d();

    const char *getClassType(void) const {return "ComponentElement3d";};
    
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
    const Matrix &getMass(void);    

    void zeroLoad(void);	
    int addLoad(ElementalLoad *theLoad, double loadFactor);
    int addInertiaLoadToUnbalance(const Vector &accel);

    const Vector &getResistingForce(void);
    const Vector &getResistingForceIncInertia(void);            
    
    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);
    
    void Print(OPS_Stream &s, int flag = 0);    
    int displaySelf(Renderer &theViewer, int displayMode, float fact, const char **modes = 0, int numModes = 0);

    Response *setResponse (const char **argv, int argc, OPS_Stream &s);
    int getResponse (int responseID, Information &info);
 
    int setParameter (const char **argv, int argc, Parameter &param);
    int updateParameter (int parameterID, Information &info);

  private:
  double A,E,Iz,Iy,G,J;     // area, elastic modulus, moment of inertia
    double rho;       // mass per unit length
    int cMass;        // consistent mass flag

    Vector Q;
    
    Vector q;
    double q0[5];  // Fixed end forces in basic system
    double p0[5];  // Reactions in basic system
    
    Node *theNodes[2];
    ID  connectedExternalNodes;    
    CrdTransf *theCoordTransf;

    UniaxialMaterial *end1zHinge;
    UniaxialMaterial *end2zHinge;
    UniaxialMaterial *end1yHinge;
    UniaxialMaterial *end2yHinge;  
  //double cD1z, cD2z, tD1z, tD2z; // committed and trial interior dof displacements
  //double cD1y, cD2y, tD1y, tD2y;
    Matrix kzTrial;
  //Vector Rz;
    Vector uzTrial;
    Vector uzCommit;
    Matrix kyTrial;
  //Vector Ry;
    Vector uyTrial;
    Vector uyCommit;  
  //Vector rTrial;
  //Vector rCommit;
    Matrix kb;

    static Vector P;
    static Matrix K;
    bool init;

  double L;
  double EAoverL;		// EA/L
  double EIzoverL2;		// 2EIz/L
  double EIzoverL4;		// 4EIz/L
  double EIyoverL2;		// 2EIy/L
  double EIyoverL4;		// 4EIy/L
  double GJoverL;               // GJ/L
};

#endif
