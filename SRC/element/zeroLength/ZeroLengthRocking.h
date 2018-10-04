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
                                                                        
// $Revision: 1.14 $
// $Date: 2010-02-04 01:17:46 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/ZeroLengthRocking/ZeroLengthRocking.h,v $
                                                                        
                                                                        
#ifndef ZeroLengthRocking_h
#define ZeroLengthRocking_h

// File: ~/element/ZeroLengthRocking/ZeroLengthRocking.h
// 
// Written: KRM
// Created: 12/2012
// Revision: A
//
// Description: This file contains the class definition for ZeroLengthRocking.
// A ZeroLengthRocking element is defined by two nodes with the same coordinate. It is based 
// on the rocking element in Barthes, C. Design of Earthquake Resistant Bridges Using Rocking 
// Columns, PhD Dissertation, UC Berkeley 2012.
//
// What: "@(#) ZeroLengthRocking.h, revA"

#include <Element.h>
#include <Matrix.h>
#include <Vector.h>

// Tolerance for zero length of element
#define	LENTOL 1.0e-6


class Node;
class Channel;
class Response;

class ZeroLengthRocking : public Element
{
  public:
    
  // Constructor
  ZeroLengthRocking(int tag, 			      
        int dimension,
        int Nd1, int Nd2, 
        const Vector& x,
        const Vector& yprime,
        double kr, double R,
        double theta, double kappa,
        double xi, double dispTol, double velTol);

    ZeroLengthRocking();    
    ~ZeroLengthRocking();

    const char *getClassType(void) const {return "ZeroLengthRocking";};

    // public methods to obtain information about dof & connectivity    
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
    int displaySelf(Renderer &, int mode, float fact, const char **displayModes=0, int numModes=0);
    void Print(OPS_Stream &s, int flag =0);    

    Response *setResponse(const char **argv, int argc, OPS_Stream &s);
    int getResponse(int responseID, Information &eleInformation);
    
    int setParameter (const char **argv, int argc, Parameter &param);
    int updateParameter (int parameterID, Information &info);

  protected:
    
  private:
    // private methods
    void   setUp ( int Nd1, int Nd2, const Vector& x, const Vector& y);

    // private attributes - a copy for each object of the class
    ID  connectedExternalNodes;         // contains the tags of the end nodes
    int dimension;                      // = 1, 2, or 3 dimensions
    int numDOF;	                        // number of dof for ZeroLengthRocking
    Matrix transformation;		// transformation matrix for orientation
	
    Node *theNodes[2];

    Matrix *theMatrix; 	    	// pointer to objects matrix (a class Matrix)
    Vector *theVector;      	// pointer to objects vector (a class Vector)

    Matrix *Llocal;
    Vector *constraint;
    Vector *vb;
    
    // element data
    double ktheta;      // rotational stiffness (set to zero for free rocking)
    double Rrock;       // rocking radius
    double Trock;       // initial orientation of rocking surface (0 for horizontal)
    double kappa;       // penalty factor
    double xi;          // bound parameter for smoothing
    double dispTol;     // displacemnt tol for deactivating rocking
    double velTol;      // velocity tol for deactivating rocking
    
    int Rocking;        // state variable
    int RockingCounter; // state variable
    double Moment;      // state variable
    double d31plusT;    // state variable

    // static data - single copy for all objects of the class	
    static Matrix ZeroLengthRockingM6;   // class wide matrix for 6*6
    static Matrix ZeroLengthRockingM12;  // class wide matrix for 12*12
    static Vector ZeroLengthRockingV6;   // class wide Vector for size 6
    static Vector ZeroLengthRockingV12;  // class wide Vector for size 12

};

#endif




