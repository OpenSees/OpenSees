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
                                                                        
// $Revision: 1.11 $
// $Date: 2006-08-04 19:08:07 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/elasticBeamColumn/ElasticBeam2d.h,v $
                                                                        
                                                                        
// Written: fmk 11/95
// Revised:
//
// Purpose: This file contains the class definition for ElasticBeam2d.
// ElasticBeam2d is a plane frame member.

#ifndef ElasticBeam2d_h
#define ElasticBeam2d_h

#include <Element.h>
#include <Node.h>
#include <Matrix.h>
#include <Vector.h>

class Channel;
class Information;
class CrdTransf2d;
class Response;
class Renderer;

class ElasticBeam2d : public Element
{
  public:
    ElasticBeam2d();        
    ElasticBeam2d(int tag, double A, double E, double I, 
		  int Nd1, int Nd2, CrdTransf2d &theTransf, 
		  double alpha = 0.0, double d = 0.0, double rho = 0.0);
    ~ElasticBeam2d();

    const char *getClassType(void) const {return "ElasticBeam2d";};

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
    int displaySelf(Renderer &theViewer, int displayMode, float fact);

    Response *setResponse (const char **argv, int argc, Information &info, OPS_Stream &s);
    int getResponse (int responseID, Information &info);
 
    int setParameter (const char **argv, int argc, Information &info);
    int updateParameter (int parameterID, Information &info);

  private:
    double A,E,I;
    double alpha, d;
    
    double rho;
    
    static Matrix K;
    static Vector P;
    Vector Q;
    
    static Matrix kb;
    Vector q;
    double q0[3];  // Fixed end forces in basic system
    double p0[3];  // Reactions in basic system
    
    Node *theNodes[2];
    
    ID  connectedExternalNodes;    

    CrdTransf2d *theCoordTransf;
};

#endif


