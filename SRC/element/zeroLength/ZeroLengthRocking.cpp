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
                                                                        
// $Revision: 1.27 $
// $Date: 2010-02-04 01:17:46 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/ZeroLengthRocking/ZeroLengthRocking.cpp,v $

// Written: KRM
// Created: 12/2012
// Revision: A
//
// Description: This file contains the implementation for the ZeroLengthRocking class.

#include "ZeroLengthRocking.h"
#include <Information.h>
#include <Parameter.h>
#include <Domain.h>
#include <Node.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Renderer.h>
#include <ElementResponse.h>
#include <Matrix.h>
#include <Vector.h>

#include <math.h>
#include <float.h>
#include <stdlib.h>
#include <string.h>
#include <elementAPI.h>

// initialise the class wide variables
Matrix ZeroLengthRocking::ZeroLengthRockingM6(6,6);
Matrix ZeroLengthRocking::ZeroLengthRockingM12(12,12);
Vector ZeroLengthRocking::ZeroLengthRockingV6(6);
Vector ZeroLengthRocking::ZeroLengthRockingV12(12);

void * OPS_ADD_RUNTIME_VPV(OPS_ZeroLengthRocking)
{
    
    int ndm = OPS_GetNDM(); // the spatial dimension of the problem
    
    //
    // first scan the command line to obtain eleID, iNode, jNode, and the orientation 
    // of ele xPrime and yPrime not along the global x and y axis
    //
    
    // a quick check on number of args
    if (OPS_GetNumRemainingInputArgs() < 7) {
        opserr << "WARNING too few arguments " <<
	    "want - element ZeroLengthRocking eleTag? iNode? jNode? " <<
	    "kr? radius? theta0? kappa? <-orient x1? x2? x3? y1? y2? y3?>\n";
        
        return 0;
    }

    // eleTag, iNode, jNode
    int idata[3];
    int numdata = 3;
    if (OPS_GetIntInput(&numdata, idata) < 0) {
	opserr << "WARNING invalied int inputs " <<
	    "- element ZeroLengthRocking eleTag? iNode? jNode? " <<
	    "kr? radius? theta0? kappa? <-orient x1? x2? x3? y1? y2? y3?>\n";
	return 0;
    }
    int eleTag = idata[0];
    int iNode = idata[1];
    int jNode = idata[2];
    
    // look for rocking required inputs
    double ddata[4];
    numdata = 4;
    if (OPS_GetDoubleInput(&numdata, ddata) < 0) {
	opserr << "WARNING invalied double inputs " <<
	    "- element ZeroLengthRocking eleTag? iNode? jNode? " <<
	    "kr? radius? theta0? kappa? <-orient x1? x2? x3? y1? y2? y3?>\n";
	return 0;
    }
    double kr = ddata[0];
    double R = ddata[1];
    double theta = ddata[2];
    double kap = ddata[3];
    

    // create the vectors for the element orientation
    Vector x(3); x(0) = 1.0; x(1) = 0.0; x(2) = 0.0;
    Vector y(3); y(0) = 0.0; y(1) = 1.0; y(2) = 0.0;
    double xi = 1.0e-8;
    double dTol = 1.0e-7;
    double vTol = 1.0e-7;
    
    while (OPS_GetNumRemainingInputArgs() > 0) {
	const char* flag = OPS_GetString();
	numdata = 1;
        if (strcmp(flag,"-orient") == 0) {
            if (OPS_GetNumRemainingInputArgs() < 6) {
                opserr << "WARNING not enough parameters after -orient flag for ele " << eleTag <<
		    "- element ZeroLengthRocking eleTag? iNode? jNode? " <<
		    "kr? radius? theta0? kappa? <-orient x1? x2? x3? y1? y2? y3?>\n";	      
                return 0;
                
            } else {
                double value;
                // read the x values
                for (int i=0; i<3; i++)  {
		    if (OPS_GetDoubleInput(&numdata, &value) < 0) {
                        opserr << "WARNING invalid -orient value for ele  " << eleTag <<
			    "- element ZeroLength eleTag? iNode? jNode? " <<
			    "kr? radius? theta0? kappa? <-orient x1? x2? x3? y1? y2? y3?>\n";
                        return 0;
                    } else {
                        x(i) = value;
                    }
                }
                // read the y values
                for (int j=0; j<3; j++)  {
		    if (OPS_GetDoubleInput(&numdata, &value) < 0) {
                        opserr << "WARNING invalid -orient value for ele  " <<
			    eleTag << 
			    "- element ZeroLength eleTag? iNode? jNode? " <<
			    "kr? radius? theta0? kappa? <-orient x1? x2? x3? y1? y2? y3?>\n";
                        return 0;
                    } else {
                        y(j) = value;
                    }
                }
            }
            
        } else if (strcmp(flag,"-xi") == 0) {
            if (OPS_GetNumRemainingInputArgs() < 1) {
                opserr << "WARNING not enough parameters after -xi flag for ele " << eleTag << endln;
                return 0;
            } else {
		if (OPS_GetDoubleInput(&numdata, &xi) < 0) {
                    opserr << "WARNING invalid -xi value for ele  " << eleTag << endln;
                    return 0;
                }
            }
            
        } else if (strcmp(flag,"-dTol") == 0) {
            if (OPS_GetNumRemainingInputArgs() < 1) {
                opserr << "WARNING not enough parameters after -dTol flag for ele " << eleTag << endln;
                return 0;
            } else {
		if (OPS_GetDoubleInput(&numdata, &dTol) < 0) {
                    opserr << "WARNING invalid -dTol value for ele  " << eleTag << endln;
                    return 0;
		}
            }

        } else if (strcmp(flag,"-vTol") == 0) {
            if (OPS_GetNumRemainingInputArgs() < 1) {
                opserr << "WARNING not enough parameters after -vTol flag for ele " << eleTag << endln;
                return 0;
            } else {
		if (OPS_GetDoubleInput(&numdata, &vTol) < 0) {
                    opserr << "WARNING invalid -vTol value for ele  " << eleTag << endln;
                    return 0;
                }
            }
            
        }
    }
    
    //
    // now we create the element and add it to the domain
    //
    return new ZeroLengthRocking(eleTag, ndm, iNode, jNode, x, y, 
				 kr, R, theta, kap, xi, dTol, vTol);
}

//  Constructor:
//  responsible for allocating the necessary space needed by each object
//  and storing the tags of the ZeroLengthRocking end nodes.
ZeroLengthRocking::ZeroLengthRocking(int tag,
		       int dim,
		       int Nd1, int Nd2, 
		       const Vector &x, const Vector &yp,
		       double kr, double R, double theta, double kap,
               double eb, double dtol, double vtol)
 :Element(tag,ELE_TAG_ZeroLengthRocking),     
  connectedExternalNodes(2),
  dimension(dim), numDOF(0), transformation(3,3), 
  theMatrix(0), theVector(0),  
  ktheta(kr), Rrock(R), Trock(theta), kappa(kap),
  xi(eb), dispTol(dtol), velTol(vtol)
{

    // establish the connected nodes and set up the transformation matrix for orientation
    this->setUp( Nd1, Nd2, x, yp);

    // do some quick parameter checking
    if (ktheta < 0) {
        opserr << "ZeroLengthRocking:: cannot have negative rocking stiffness, setting to zero" 
            << endln;
        ktheta = 0;
    }
    
    // initialization
    Rocking = 0;
    RockingCounter = 0;
    Moment = 0;
    d31plusT = 0;
}


//   Constructor:
//   invoked by a FEM_ObjectBroker - blank object that recvSelf needs to be invoked upon
ZeroLengthRocking::ZeroLengthRocking(void)
  :Element(0,ELE_TAG_ZeroLengthRocking),     
  connectedExternalNodes(2),
  dimension(0), numDOF(0), transformation(3,3),
  theMatrix(0), theVector(0),
  ktheta(0), Rrock(0), Trock(0), kappa(0), 
  xi(0), dispTol(0), velTol(0)
{
    // ensure the connectedExternalNode ID is of correct size 
    if (connectedExternalNodes.Size() != 2)
      opserr << "FATAL ZeroLengthRocking::ZeroLengthRocking - failed to create an ID of correct size\n";

}


//  Destructor: 
//  delete must be invoked on any objects created by the object
ZeroLengthRocking::~ZeroLengthRocking()
{
    if (Llocal != 0)
        delete Llocal;
    if (constraint != 0)
        delete constraint;
    if (vb != 0)
        delete vb;
}


int
ZeroLengthRocking::getNumExternalNodes(void) const
{
    return 2;
}


const ID &
ZeroLengthRocking::getExternalNodes(void) 
{
    return connectedExternalNodes;
}


Node **
ZeroLengthRocking::getNodePtrs(void) 
{
  return theNodes;
}


int
ZeroLengthRocking::getNumDOF(void) 
{
    return numDOF;
}


// method: setDomain()
//    to set a link to the enclosing Domain and to set the node pointers.
//    also determines the number of dof associated
//    with the ZeroLengthRocking element, we set matrix and vector pointers,
//    allocate space for t matrix and define it as the basic deformation-
//    displacement transformation matrix.
void
ZeroLengthRocking::setDomain(Domain *theDomain)
{
    // check Domain is not null - invoked when object removed from a domain
    if (theDomain == 0) {
        theNodes[0] = 0;
        theNodes[1] = 0;
        return;
    }

    // set default values for error conditions
    numDOF = 3;
    theMatrix = &ZeroLengthRockingM6;
    theVector = &ZeroLengthRockingV6;
    
    // first set the node pointers
    int Nd1 = connectedExternalNodes(0);
    int Nd2 = connectedExternalNodes(1);
    theNodes[0] = theDomain->getNode(Nd1);
    theNodes[1] = theDomain->getNode(Nd2);	

    // if can't find both - send a warning message
    if ( theNodes[0] == 0 || theNodes[1] == 0 ) {
      if (theNodes[0] == 0) 
        opserr << "WARNING ZeroLengthRocking::setDomain() - Nd1: " << Nd1 << " does not exist in ";
      else
        opserr << "WARNING ZeroLengthRocking::setDomain() - Nd2: " << Nd2 << " does not exist in ";

      opserr << "model for ZeroLengthRocking ele: " << this->getTag() << endln;

      return;
    }

    // now determine the number of dof and the dimension    
    int dofNd1 = theNodes[0]->getNumberDOF();
    int dofNd2 = theNodes[1]->getNumberDOF();	

    // if differing dof at the ends - print a warning message
    if ( dofNd1 != dofNd2 ) {
      opserr << "WARNING ZeroLengthRocking::setDomain(): nodes " << Nd1 << " and " << Nd2 <<
            "have differing dof at ends for ZeroLengthRocking " << this->getTag() << endln;
      return;
    }	

    // Check that length is zero within tolerance
    const Vector &end1Crd = theNodes[0]->getCrds();
    const Vector &end2Crd = theNodes[1]->getCrds();	
    Vector diff = end1Crd - end2Crd;
    double L  = diff.Norm();
    double v1 = end1Crd.Norm();
    double v2 = end2Crd.Norm();
    double vm;
    
    vm = (v1<v2) ? v2 : v1;

    if (L > LENTOL*vm)
        opserr << "WARNING ZeroLengthRocking::setDomain(): Element " << this->getTag() << " has L= " << L 
               << ", which is greater than the tolerance\n";
        
    // call the base class method
    this->DomainComponent::setDomain(theDomain);
    
    // set the number of dof for element and set matrix and vector pointer
    if (dimension == 2 && dofNd1 == 3) {
        numDOF = 6;	
        theMatrix = &ZeroLengthRockingM6;
        theVector = &ZeroLengthRockingV6;
        
        Llocal = new Matrix(2,6);
        constraint = new Vector(2);
        vb = new Vector(1);
    }
    else if (dimension == 3 && dofNd1 == 6) {
        numDOF = 12;	    
        theMatrix = &ZeroLengthRockingM12;
        theVector = &ZeroLengthRockingV12;
        
        Llocal = new Matrix(4,12);
        constraint = new Vector(4);
        vb = new Vector(3);
    }
    else {
        opserr << "WARNING ZeroLengthRocking::setDomain cannot handle " << dimension << 
                  "dofs at nodes in " << dofNd1 << " d problem\n"; 
        return;
    }
   
}   	 


int
ZeroLengthRocking::commitState()
{
    // get trial displacements and take difference
    const Vector& disp1 = theNodes[0]->getDisp();
    const Vector& disp2 = theNodes[1]->getDisp();
    Vector  diff  = disp2-disp1;
    const Vector& vel1  = theNodes[0]->getVel();
    const Vector& vel2  = theNodes[1]->getVel();
    Vector  diffv = vel2-vel1;

    // increment rocking counter
    RockingCounter++;
    
    if ( Rocking == 0 ) {
        // check for activation
        if ( (Moment > 0) && (RockingCounter > 10) ) {
            opserr << "Rocking activated in element " << this->getTag() << " with counter at "
                   << RockingCounter << endln;
            RockingCounter = 0;
            Rocking = 1;
        }
    } else {
        // check for deactivation
        if ( (fabs(diff(2)) <= dispTol) && (fabs(diffv(2)) <= velTol) && (RockingCounter >= 50) ) {
            opserr << "Rocking deactivated in element " << this->getTag() << " with counter at "
                   << RockingCounter << endln;
            RockingCounter = 0;
            Rocking = 0;
        }
    }
    
    int code = 0;
    
    // call element commitState to do any base class stuff
    if ((code = this->Element::commitState()) != 0) {
        opserr << "ZeroLengthRocking::commitState () - failed in base class";
    }    

    return code;
}


int
ZeroLengthRocking::revertToLastCommit()
{
    int code = 0;
    // ZeroLengthRocking is memoryless
    return code;
}


int
ZeroLengthRocking::revertToStart()
{   
    int code = 0;
    
    Rocking = 0;
    RockingCounter = 0;
    Moment = 0;
    d31plusT = 0;
    
    return code;
}


int
ZeroLengthRocking::update(void)
{
    // get trial displacements and take difference
    const Vector& disp1 = theNodes[0]->getTrialDisp();
    const Vector& disp2 = theNodes[1]->getTrialDisp();
    Vector  diff  = disp2-disp1;
    const Vector& vel1  = theNodes[0]->getTrialVel();
    const Vector& vel2  = theNodes[1]->getTrialVel();
    Vector  diffv = vel2-vel1;
    
    // FROM HERE ON specific to 2D problem, need to fail properly if user attempts 3D until 
    // someone implements it......
    
    // need to store basic deformations for later calls
    (*vb)(0) = diff(2);
    
    // now compute and store rocking and constraint matrices
    double RockingSign = ( diff(2) > 0 ? 1.0:-1.0 );
    if ( diff(2) == 0.0 )
        RockingSign = 0;
    
    d31plusT = disp1(2)+Trock;
    
    // need to size and store Llocal and constraint in setDomain
    (*Llocal)(0,0) = -cos(d31plusT);
    (*Llocal)(0,1) = -sin(d31plusT);
    (*Llocal)(0,2) = -diff(0)*sin(d31plusT) + diff(1)*cos(d31plusT) - RockingSign*Rrock*sin(diff(2));
    (*Llocal)(0,3) = cos(d31plusT);
    (*Llocal)(0,4) = sin(d31plusT);
    (*Llocal)(0,5) = RockingSign*Rrock*sin(diff(2));
    
    (*Llocal)(1,0) = sin(d31plusT);
    (*Llocal)(1,1) = -cos(d31plusT);
    (*Llocal)(1,2) = -diff(0)*cos(d31plusT) - diff(1)*sin(d31plusT) + RockingSign*Rrock*cos(diff(2));
    (*Llocal)(1,3) = -sin(d31plusT);
    (*Llocal)(1,4) = cos(d31plusT);
    (*Llocal)(1,5) = -RockingSign*Rrock*cos(diff(2));
    
    (*constraint)(0) = diff(0)*cos(d31plusT) + diff(1)*sin(d31plusT) + RockingSign*Rrock*(1-cos(diff(2)));
    (*constraint)(1) = -diff(0)*sin(d31plusT) + diff(1)*cos(d31plusT) - RockingSign*Rrock*sin(diff(2));
    
    if ( fabs(diff(2)) >= xi ) {
        // do nothing, already defined above
        
    } else {
        double eps2 = xi*xi;
        double alpha = -1./(8.*eps2)*sin(xi)-1./(8.*xi*eps2)*cos(xi);
        double beta = -1./2.*sin(xi) - 6.*alpha*eps2;
        double gamma = sin(xi) - alpha*eps2*eps2 - beta*eps2;
        
        (*Llocal)(1,2) = -diff(0)*cos(d31plusT) - diff(1)*sin(d31plusT) + Rrock*(4.*alpha*pow(diff(2),3) + 2.*beta*diff(2));
        (*Llocal)(1,5) = -Rrock*(4.*alpha*pow(diff(2),3) + 2.*beta*diff(2));
        
        (*constraint)(1) = -diff(0)*sin(d31plusT) + diff(1)*cos(d31plusT) - Rrock*(alpha*pow(diff(2),4) + beta*pow(diff(2),2) + gamma);
    }
    
    return 0;
}


const Matrix &
ZeroLengthRocking::getTangentStiff(void)
{
    // stiff is a reference to the matrix holding the stiffness matrix
    Matrix& stiff = *theMatrix;
    
    // transform basic to global
    stiff.addMatrixTransposeProduct(0.0,*Llocal,*Llocal,kappa);

    // mod rotation terms
    stiff(2,2) += ktheta;
    stiff(2,5) -= ktheta;
    stiff(5,2) -= ktheta;
    stiff(5,5) += ktheta;
    
    if (Rocking == 0) {
        stiff(2,2) += kappa;
        stiff(2,5) -= kappa;
        stiff(5,2) -= kappa;
        stiff(5,5) += kappa;
    }
    
    return stiff;
}


const Matrix &
ZeroLengthRocking::getInitialStiff(void)
{
    // for now do the same thing as tangent stiffness, but in the future 
    // compute the exact Llocal with all the sin and cos terms that go to zero
    
    // stiff is a reference to the matrix holding the stiffness matrix
    Matrix& stiff = *theMatrix;
    
    // transform basic to global
    stiff.addMatrixTransposeProduct(0.0,*Llocal,*Llocal,kappa);
    
    // mod rotation terms
    stiff(2,2) += ktheta;
    stiff(2,5) -= ktheta;
    stiff(5,2) -= ktheta;
    stiff(5,5) += ktheta;
    
    if (Rocking == 0) {
        stiff(2,2) += kappa;
        stiff(2,5) -= kappa;
        stiff(5,2) -= kappa;
        stiff(5,5) += kappa;
    }
    
    return stiff;
}
    

const Matrix &
ZeroLengthRocking::getDamp(void)
{
    // NYI
    // damp is a reference to the matrix holding the damping matrix
    Matrix& damp = *theMatrix;

    // zero damping matrix
    damp.Zero();

    return damp;
}


const Matrix &
ZeroLengthRocking::getMass(void)
{
    // no mass 
    theMatrix->Zero();    
    return *theMatrix; 
}


void 
ZeroLengthRocking::zeroLoad(void)
{
    // does nothing now
}

int 
ZeroLengthRocking::addLoad(ElementalLoad *theLoad, double loadFactor)
{
    opserr << "ZeroLengthRocking::addLoad - load type unknown for element with tag: " << this->getTag() << endln;

    return -1;
}

int 
ZeroLengthRocking::addInertiaLoadToUnbalance(const Vector &accel)
{
    // does nothing as element has no mass yet!
    return 0;
}


const Vector &
ZeroLengthRocking::getResistingForce()
{
    // force is a reference to the vector holding the resisting force
    Vector& force = *theVector;

    // basic to global
    force.addMatrixTransposeVector(0.0,*Llocal,*constraint,kappa);
    
    // mod rotation terms
    force(2) -= ktheta * (*vb)(0);
    force(5) += ktheta * (*vb)(0);
    
    if (Rocking == 0) {
        force(2) -= kappa * (*vb)(0);
        force(5) += kappa * (*vb)(0);
    }
    
    // update moment state
    Moment = fabs(force(5)-force(2)) - Rrock*sin(d31plusT)*(force(3)-force(0)) + Rrock*cos(d31plusT)*(force(4)-force(1));
    
    return force;
}


const Vector &
ZeroLengthRocking::getResistingForceIncInertia()
{	
    // we don't have any inertia
    return this->getResistingForce();

}


int
ZeroLengthRocking::sendSelf(int commitTag, Channel &theChannel)
{
	int res = 0;

	// note: we don't check for dataTag == 0 for Element
	// objects as that is taken care of in a commit by the Domain
	// object - don't want to have to do the check if sending data
	int dataTag = this->getDbTag();

	// ZeroLengthRocking packs its data into an ID and sends this to theChannel
	// along with its dbTag and the commitTag passed in the arguments
    static ID idData(7);

	idData(0) = this->getTag();
	idData(1) = dimension;
	idData(2) = numDOF;
	idData(3) = connectedExternalNodes(0);
	idData(4) = connectedExternalNodes(1);
    idData(5) = Rocking;
    idData(6) = RockingCounter;

	res += theChannel.sendID(dataTag, commitTag, idData);
	if (res < 0) {
	  opserr << "ZeroLengthRocking::sendSelf -- failed to send ID data\n";
	  return res;
	}
    
    // now send vector of double data
    static Vector dData(9);
    
    dData(0) = ktheta;
    dData(1) = Rrock;
    dData(2) = Trock;
    dData(3) = kappa;
    dData(4) = xi;
    dData(5) = dispTol;
    dData(6) = velTol;
    dData(7) = Moment;
    dData(8) = d31plusT;
    
    res += theChannel.sendVector(dataTag, commitTag, dData);
	if (res < 0) {
        opserr << "ZeroLengthRocking::sendSelf -- failed to send Vector data\n";
        return res;
	}
    
    // note we are not sending transformation because it is NYI

	return res;
}


int
ZeroLengthRocking::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
    int res = 0;

    int dataTag = this->getDbTag();

    // ZeroLengthRocking creates a ID, receives the ID and then sets the 
    // internal data with the data in the ID
    static ID idData(7);
    res += theChannel.recvID(dataTag, commitTag, idData);
    if (res < 0) {
        opserr << "ZeroLengthRocking::recvSelf -- failed to receive ID data\n";
        return res;
    }

    this->setTag(idData(0));
    dimension = idData(1);
    numDOF = idData(2);
    connectedExternalNodes(0) = idData(3);
    connectedExternalNodes(1) = idData(4);
    Rocking = idData(5);
    RockingCounter = idData(6);
    
    // now receive Vector
    static Vector dData(9);
    res += theChannel.recvVector(dataTag, commitTag, dData);
    if (res < 0) {
        opserr << "ZeroLengthRocking::recvSelf -- failed to receive Vector data\n";
        return res;
    }
    
    ktheta = dData(0);
    Rrock = dData(1);
    Trock = dData(2);
    kappa = dData(3);
    xi = dData(4);
    dispTol = dData(5);
    velTol = dData(6);
    Moment = dData(7);
    d31plusT = dData(8);

    return res;
}


int
ZeroLengthRocking::displaySelf(Renderer &theViewer, int displayMode, float fact, const char **modes, int numMode)
{
    // ensure setDomain() worked
    if (theNodes[0] == 0 || theNodes[1] == 0 )
       return 0;

    // get the end point display coords    
    static Vector v1(3);
    static Vector v2(3);
    theNodes[0]->getDisplayCrds(v1, fact, displayMode);
    theNodes[1]->getDisplayCrds(v2, fact, displayMode);

    // draw the line (don't display strain or force)
    return theViewer.drawLine(v1, v2, 0.0, 0.0, this->getTag());
}


void
ZeroLengthRocking::Print(OPS_Stream &s, int flag)
{
    if (flag == 0) { // print everything
        s << "Element: " << this->getTag(); 
        s << " type: ZeroLengthRocking  iNode: " << connectedExternalNodes(0);
        s << " jNode: " << connectedExternalNodes(1) << endln;
        s << " Moment: " << Moment << " and rocking state: " << Rocking << endln;
    } else if (flag == 1) {
        s << this->getTag() << "  " << vb << "  ";
    }
}


Response*
ZeroLengthRocking::setResponse(const char **argv, int argc, OPS_Stream &output)
{
    Response *theResponse = 0;

    output.tag("ElementOutput");
    output.attr("eleType","ZeroLengthRocking");
    output.attr("eleTag",this->getTag());
    output.attr("node1",connectedExternalNodes[0]);
    output.attr("node2",connectedExternalNodes[1]);

    char outputData[10];

    if ((strcmp(argv[0],"force") == 0) || (strcmp(argv[0],"forces") == 0)) {

            int numDOFperNode = numDOF/2;
            for (int i=0; i<numDOFperNode; i++) {
                sprintf(outputData,"P1_%d", i+1);
                output.tag("ResponseType", outputData);
            }
            for (int j=0; j<numDOFperNode; j++) {
                sprintf(outputData,"P2_%d", j+1);
                output.tag("ResponseType", outputData);
            }
            theResponse = new ElementResponse(this, 1, Vector(numDOF));

    } else if (strcmp(argv[0],"defo") == 0 || strcmp(argv[0],"deformations") == 0 ||
	       strcmp(argv[0],"deformation") == 0 || strcmp(argv[0],"basicDeformation") == 0) {

            for (int i=0; i<vb->Size(); i++) {
                sprintf(outputData,"vb%d",i+1);
                output.tag("ResponseType",outputData);
            }
            theResponse = new ElementResponse(this, 2, Vector(vb->Size()));

    }
    
    output.endTag();

    return theResponse;
}


int 
ZeroLengthRocking::getResponse(int responseID, Information &eleInformation)
{
    const Vector& disp1 = theNodes[0]->getTrialDisp();
    const Vector& disp2 = theNodes[1]->getTrialDisp();
    const Vector  diff  = disp2-disp1;

    switch (responseID) {
    case -1:
        return -1;

    case 1:
        return eleInformation.setVector(this->getResistingForce());

    case 2:
        return eleInformation.setVector(*vb);

    default:
        return -1;
    }
}


// Private methods


// Establish the external nodes and set up the transformation matrix for orientation
void
ZeroLengthRocking::setUp( int Nd1, int Nd2,
		   const Vector &x,
		   const Vector &yp )
{ 
    // ensure the connectedExternalNode ID is of correct size & set values
    if (connectedExternalNodes.Size() != 2)
      opserr << "FATAL ZeroLengthRocking::setUp - failed to create an ID of correct size\n";
    
    connectedExternalNodes(0) = Nd1;
    connectedExternalNodes(1) = Nd2;

    for (int i=0; i<2; i++)
      theNodes[i] = 0;

    // check that vectors for orientation are correct size
    if ( x.Size() != 3 || yp.Size() != 3 )
	opserr << "FATAL ZeroLengthRocking::setUp - incorrect dimension of orientation vectors\n";

    // establish orientation of element for the transformation matrix
    // z = x cross yp
    Vector z(3);
    z(0) = x(1)*yp(2) - x(2)*yp(1);
    z(1) = x(2)*yp(0) - x(0)*yp(2);
    z(2) = x(0)*yp(1) - x(1)*yp(0);

    // y = z cross x
    Vector y(3);
    y(0) = z(1)*x(2) - z(2)*x(1);
    y(1) = z(2)*x(0) - z(0)*x(2);
    y(2) = z(0)*x(1) - z(1)*x(0);

    // compute length(norm) of vectors
    double xn = x.Norm();
    double yn = y.Norm();
    double zn = z.Norm();

    // check valid x and y vectors, i.e. not parallel and of zero length
    if (xn == 0 || yn == 0 || zn == 0) {
      opserr << "FATAL ZeroLengthRocking::setUp - invalid vectors to constructor\n";
    }
    
    // create transformation matrix of direction cosines
    for (int i=0; i<3; i++) {
        transformation(0,i) = x(i)/xn;
        transformation(1,i) = y(i)/yn;
        transformation(2,i) = z(i)/zn;
    }

}


int
ZeroLengthRocking::setParameter(const char **argv, int argc, Parameter &param)
{
    if (argc < 1)
        return -1;
    
    // rotational stiffness
    if (strcmp(argv[0],"kr") == 0)
        return param.addObject(1, this);
    
    // kappa
    if (strcmp(argv[0],"kappa") == 0)
        return param.addObject(2, this);
    
    // xi
    if (strcmp(argv[0],"xi") == 0)
        return param.addObject(3, this);
    
    return -1;
}


int
ZeroLengthRocking::updateParameter (int parameterID, Information &info)
{
	switch (parameterID) {
        case -1:
            return -1;
        case 1:
            ktheta = info.theDouble;
            return 0;
        case 2:
            kappa = info.theDouble;
            return 0;
        case 3:
            xi = info.theDouble;
            return 0;
        default:
            return -1;
	}
}
