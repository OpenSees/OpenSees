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
                                                                        
// $Revision: 1.2 $
// $Date: 2008/12/10 00:05:21 $
// $Source: /usr/local/cvs/OpenSees/PACKAGES/NewElement/cpp/Truss2D.h,v $

// Written by: Gonzalo Torrisi, Universidad Nacional de Cuyo

#ifndef MasonPan12_h
#define MasonPan12_h

#include <Element.h>
#include <Matrix.h>
#include <Vector.h>

class UniaxialMaterial;

class MasonPan12 : public Element
{
  public:
    // constructors
    MasonPan12(int tag, 
              int Nd1, int Nd2, int Nd3, int Nd4, int Nd5, int Nd6,
              int Nd7, int Nd8, int Nd9, int Nd10, int Nd11, int Nd12,
             UniaxialMaterial &theMaterial, UniaxialMaterial &theMaterial2,
             double thick,  double wr,double w1);



   MasonPan12();    
    
    // destructor
    ~MasonPan12();

    
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

    // public methods to obtain stiffness
    const Matrix &getTangentStiff(void);
    const Matrix &getInitialStiff(void);

    // public method to obtain resisting force
    const Vector &getResistingForce(void);

    // public methods for output    
    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker); 
     int displaySelf(Renderer &theViewer, int displayMode, float fact, const char **modes, int numMode);
    void Print(OPS_Stream &s, int flag =0);    

    Response *setResponse(const char **argv, int argc, OPS_Stream &s);
    int getResponse(int responseID, Information &eleInformation);

  protected:
    
  private:
    // private member functions - only available to objects of the class
    double computeCurrentStrain(int mat) const;
    
    // private attributes - a copy for each object of the class
    UniaxialMaterial **theMaterial;       // pointer to a material
	UniaxialMaterial *theMaterial2;
    ID  externalNodes;          	 // contains the id's of end nodes
    Matrix trans;       // hold the transformation matrix
	Vector rig1;
	Vector rig2;
	Vector rig3;
	double TH;
	double W1;
	double WR;

	double Cdeltares;
		double Tdeltares;
		    int numDOF;	                    // number of dof for truss

	Node *theNodes[12];  // node pointers

    // static data - single copy for all objects of the class
    static Matrix PanelK;   // class wide matrix for returning stiffness
    static Vector PanelR;   // class wide vector for returning residual
};
#endif

