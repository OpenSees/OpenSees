/* ****************************************************************** **
**    Openers - Open System for Earthquake Engineering Simulation    **
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
                                                                        
// $Revision: 6715 $
// $Date: 2018-05-03 06:39:06 -0700 (Thu, 03 May 2018) $
// $URL: svn://peera.berkeley.edu/usr/local/svn/OpenSees/trunk/SRC/element/elasticBeamColumn/ElasticBeam2d.cpp $
                                                                        
                                                                        
// File: ~/model/ElasticBeam2d.C
//
// Written: fmk 11/95
// Revised:
//
// Purpose: This file contains the class definition for ElasticBeam2d.
// ElasticBeam2d is a 3d beam element. As such it can only
// connect to a node with 6-dof. 

#include <ElasticBeam2d.h>
#include <ElementalLoad.h>

#include <Domain.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>

#include <CrdTransf.h>
#include <SectionForceDeformation.h>
#include <Information.h>
#include <Parameter.h>
#include <ElementResponse.h>
#include <Renderer.h>

#include <math.h>
#include <stdlib.h>
#include <elementAPI.h>
#include <string>
#include <ElementIter.h>

Matrix ElasticBeam2d::K(6,6);
Vector ElasticBeam2d::P(6);
Matrix ElasticBeam2d::kb(3,3);

void* OPS_ElasticBeam2d(const ID &info)
{
    if(OPS_GetNumRemainingInputArgs() < 5) {
	opserr<<"insufficient arguments:eleTag,iNode,jNode,<A,E,Iz>or<sectionTag>,transfTag\n";
	return 0;
    }

    int ndm = OPS_GetNDM();
    int ndf = OPS_GetNDF();
    if(ndm != 2 || ndf != 3) {
	opserr<<"ndm must be 2 and ndf must be 3\n";
	return 0;
    }

    // inputs: 
    int iData[3];
    int numData = 3;
    if(OPS_GetIntInput(&numData,&iData[0]) < 0) {
	opserr<<"WARNING failed to read integers\n";
	return 0;
    }

    bool section = false;
    int sectionTag;
    double data[3];
    if (OPS_GetNumRemainingInputArgs() > 3) {
      // Read A, E, Iz
      numData = 3;
      if(OPS_GetDoubleInput(&numData,&data[0]) < 0) {
	opserr<<"WARNING failed to read doubles\n";
	return 0;
      }
    } else {
      // Read a section tag
      numData = 1;
      if(OPS_GetIntInput(&numData,&sectionTag) < 0) {
	opserr<<"WARNING sectionTag is not integer\n";
	return 0;
      }
      section = true;
    }
    numData = 1;
    int transfTag;
    if(OPS_GetIntInput(&numData,&transfTag) < 0) {
	opserr<<"WARNING transfTag is not integer\n";
	return 0;
    }
    
    // options
    double mass = 0.0, alpha=0.0, depth=0.0;
    int cMass = 0;
    while(OPS_GetNumRemainingInputArgs() > 0) {
	std::string type = OPS_GetString();
	if(type == "-alpha") {
	    if(OPS_GetNumRemainingInputArgs() > 0) {
		if(OPS_GetDoubleInput(&numData,&alpha) < 0) return 0;
	    }
	} else if(type == "-depth") {
	    if(OPS_GetNumRemainingInputArgs() > 0) {
		if(OPS_GetDoubleInput(&numData,&depth) < 0) return 0;
	    }

	} else if(type == "-mass") {
	    if(OPS_GetNumRemainingInputArgs() > 0) {
		if(OPS_GetDoubleInput(&numData,&mass) < 0) return 0;
	    }
	} else if(type == "-cMass") {
	    cMass = 1;
	}
    }

    // check transf
    CrdTransf* theTransf = OPS_getCrdTransf(transfTag);
    if(theTransf == 0) {
	opserr<<"coord transfomration not found\n";
	return 0;
    }

    if (section) {
      SectionForceDeformation *theSection = OPS_getSectionForceDeformation(sectionTag);
      if (theSection == 0) {
	opserr << "section not found\n";
	return 0;
      }
      return new ElasticBeam2d(iData[0],iData[1],iData[2],*theSection,
			       *theTransf,alpha,depth,mass,cMass);
    } else {
      return new ElasticBeam2d(iData[0],data[0],data[1],data[2],iData[1],iData[2],
			       *theTransf,alpha,depth,mass,cMass);
    }
}

int OPS_ElasticBeam2d(Domain& theDomain, const ID& elenodes, ID& eletags)
{
    if(OPS_GetNumRemainingInputArgs() < 4) {
	opserr<<"insufficient arguments:A,E,Iz,transfTag\n";
	return -1;
    }

    // inputs: 
    double data[3];
    int numData = 3;
    if(OPS_GetDoubleInput(&numData,&data[0]) < 0) return -1;

    numData = 1;
    int transfTag;
    if(OPS_GetIntInput(&numData,&transfTag) < 0) return -1;
    
    // options
    double mass = 0.0, alpha=0.0, depth=0.0;
    int cMass = 0;
    while(OPS_GetNumRemainingInputArgs() > 0) {
	std::string type = OPS_GetString();
	if(type == "-alpha") {
	    if(OPS_GetNumRemainingInputArgs() > 0) {
		if(OPS_GetDoubleInput(&numData,&alpha) < 0) return -1;
	    }
	} else if(type == "-depth") {
	    if(OPS_GetNumRemainingInputArgs() > 0) {
		if(OPS_GetDoubleInput(&numData,&depth) < 0) return -1;
	    }

	} else if(type == "-mass") {
	    if(OPS_GetNumRemainingInputArgs() > 0) {
		if(OPS_GetDoubleInput(&numData,&mass) < 0) return -1;
	    }
	} else if(type == "-cMass") {
	    cMass = 1;
	}
    }

    // check transf
    CrdTransf* theTransf = OPS_getCrdTransf(transfTag);
    if(theTransf == 0) {
	opserr<<"coord transfomration not found\n";
	return -1;
    }

    // create elements
    ElementIter& theEles = theDomain.getElements();
    Element* theEle = theEles();
    int currTag = 0;
    if (theEle != 0) {
	currTag = theEle->getTag();
    }
    eletags.resize(elenodes.Size()/2);
    for (int i=0; i<elenodes.Size()/2; i++) {
	theEle = new ElasticBeam2d(--currTag,data[0],data[1],data[2],elenodes(2*i),elenodes(2*i+1),
				   *theTransf,alpha,depth,mass,cMass);
	if (theEle == 0) {
	    opserr<<"WARING: run out of memory for creating element\n";
	    return -1;
	}
	if (theDomain.addElement(theEle) == false) {
	    opserr<<"WARNING: failed to add element to domain\n";
	    delete theEle;
	    return -1;
	}
	eletags(i) = currTag;
    }

    return 0;
}

ElasticBeam2d::ElasticBeam2d()
  :Element(0,ELE_TAG_ElasticBeam2d), 
  A(0.0), E(0.0), I(0.0), alpha(0.0), d(0.0), rho(0.0), cMass(0),
  Q(6), q(3), connectedExternalNodes(2), theCoordTransf(0)
{
  // does nothing
  q0[0] = 0.0;
  q0[1] = 0.0;
  q0[2] = 0.0;

  p0[0] = 0.0;
  p0[1] = 0.0;
  p0[2] = 0.0;

  // set node pointers to NULL
  for (int i=0; i<2; i++)
    theNodes[i] = 0;      
}

ElasticBeam2d::ElasticBeam2d(int tag, double a, double e, double i, 
			     int Nd1, int Nd2, CrdTransf &coordTransf,
			     double Alpha, double depth, double r, int cm)
  :Element(tag,ELE_TAG_ElasticBeam2d), 
  A(a), E(e), I(i), alpha(Alpha), d(depth), rho(r), cMass(cm),
  Q(6), q(3), connectedExternalNodes(2), theCoordTransf(0)
{
  connectedExternalNodes(0) = Nd1;
  connectedExternalNodes(1) = Nd2;
    
  theCoordTransf = coordTransf.getCopy2d();
    
  if (!theCoordTransf) {
    opserr << "ElasticBeam2d::ElasticBeam2d -- failed to get copy of coordinate transformation\n";
    exit(01);
  }

  q0[0] = 0.0;
  q0[1] = 0.0;
  q0[2] = 0.0;

  p0[0] = 0.0;
  p0[1] = 0.0;
  p0[2] = 0.0;

  // set node pointers to NULL
  theNodes[0] = 0;
  theNodes[1] = 0;
}

ElasticBeam2d::ElasticBeam2d(int tag, int Nd1, int Nd2, SectionForceDeformation &section,  
			     CrdTransf &coordTransf, double Alpha, double depth, double r, int cm)
  :Element(tag,ELE_TAG_ElasticBeam2d), alpha(Alpha), d(depth), rho(r), cMass(cm),
  Q(6), q(3), connectedExternalNodes(2), theCoordTransf(0)
{
  E = 1.0;
  rho = r;
  cMass = cm;

  const Matrix &sectTangent = section.getInitialTangent();
  const ID &sectCode = section.getType();
  for (int i=0; i<sectCode.Size(); i++) {
    int code = sectCode(i);
    switch(code) {
    case SECTION_RESPONSE_P:
      A = sectTangent(i,i);
      break;
    case SECTION_RESPONSE_MZ:
      I = sectTangent(i,i);
      break;
    default:
      break;
    }
  }
  
  connectedExternalNodes(0) = Nd1;
  connectedExternalNodes(1) = Nd2;
  
  theCoordTransf = coordTransf.getCopy2d();
  
  if (!theCoordTransf) {
    opserr << "ElasticBeam2d::ElasticBeam2d -- failed to get copy of coordinate transformation\n";
    exit(-1);
  }

  q0[0] = 0.0;
  q0[1] = 0.0;
  q0[2] = 0.0;

  p0[0] = 0.0;
  p0[1] = 0.0;
  p0[2] = 0.0;

  // set node pointers to NULL
  theNodes[0] = 0;
  theNodes[1] = 0;      
}

ElasticBeam2d::~ElasticBeam2d()
{
    if (theCoordTransf)
	delete theCoordTransf;
}

int
ElasticBeam2d::getNumExternalNodes(void) const
{
    return 2;
}

const ID &
ElasticBeam2d::getExternalNodes(void) 
{
    return connectedExternalNodes;
}

Node **
ElasticBeam2d::getNodePtrs(void) 
{
  return theNodes;
}

int
ElasticBeam2d::getNumDOF(void)
{
    return 6;
}

void
ElasticBeam2d::setDomain(Domain *theDomain)
{
  if (theDomain == 0) {
    opserr << "ElasticBeam2d::setDomain -- Domain is null\n";
    exit(-1);
  }
    
    theNodes[0] = theDomain->getNode(connectedExternalNodes(0));
    theNodes[1] = theDomain->getNode(connectedExternalNodes(1));    
    
    if (theNodes[0] == 0) {
      opserr << "ElasticBeam2d::setDomain -- Node 1: " << connectedExternalNodes(0) << " does not exist\n";
      exit(-1);
    }
			      
    if (theNodes[1] == 0) {
      opserr << "ElasticBeam2d::setDomain -- Node 2: " << connectedExternalNodes(1) << " does not exist\n";
      exit(-1);
    }

    int dofNd1 = theNodes[0]->getNumberDOF();
    int dofNd2 = theNodes[1]->getNumberDOF();    
    
    if (dofNd1 != 3) {
      opserr << "ElasticBeam2d::setDomain -- Node 1: " << connectedExternalNodes(0) 
	     << " has incorrect number of DOF\n";
      exit(-1);
    }
    
    if (dofNd2 != 3) {
      opserr << "ElasticBeam2d::setDomain -- Node 2: " << connectedExternalNodes(1) 
	     << " has incorrect number of DOF\n";
      exit(-1);
    }
	
    this->DomainComponent::setDomain(theDomain);
    
    if (theCoordTransf->initialize(theNodes[0], theNodes[1]) != 0) {
	opserr << "ElasticBeam2d::setDomain -- Error initializing coordinate transformation\n";
	exit(-1);
    }
    
    double L = theCoordTransf->getInitialLength();

    if (L == 0.0) {
      opserr << "ElasticBeam2d::setDomain -- Element has zero length\n";
      exit(-1);
    }
}

int
ElasticBeam2d::commitState()
{
  int retVal = 0;
  // call element commitState to do any base class stuff
  if ((retVal = this->Element::commitState()) != 0) {
    opserr << "ElasticBeam2d::commitState () - failed in base class";
  }    
  retVal += theCoordTransf->commitState();
  return retVal;
}

int
ElasticBeam2d::revertToLastCommit()
{
    return theCoordTransf->revertToLastCommit();
}

int
ElasticBeam2d::revertToStart()
{
    return theCoordTransf->revertToStart();
}

int
ElasticBeam2d::update(void)
{
  return theCoordTransf->update();
}

const Matrix &
ElasticBeam2d::getTangentStiff(void)
{
  const Vector &v = theCoordTransf->getBasicTrialDisp();
  
  double L = theCoordTransf->getInitialLength();

  double EoverL   = E/L;
  double EAoverL  = A*EoverL;			// EA/L
  double EIoverL2 = 2.0*I*EoverL;		// 2EI/L
  double EIoverL4 = 2.0*EIoverL2;		// 4EI/L
  
  // determine q = kv + q0
  q(0) = EAoverL*v(0);
  q(1) = EIoverL4*v(1) + EIoverL2*v(2);
  q(2) = EIoverL2*v(1) + EIoverL4*v(2);

  q(0) += q0[0];
  q(1) += q0[1];
  q(2) += q0[2];
  
  kb(0,0) = EAoverL;
  kb(1,1) = kb(2,2) = EIoverL4;
  kb(2,1) = kb(1,2) = EIoverL2;

  return theCoordTransf->getGlobalStiffMatrix(kb, q);
}

const Matrix &
ElasticBeam2d::getInitialStiff(void)
{
  double L = theCoordTransf->getInitialLength();

  double EoverL   = E/L;
  double EAoverL  = A*EoverL;			// EA/L
  double EIoverL2 = 2.0*I*EoverL;		// 2EI/L
  double EIoverL4 = 2.0*EIoverL2;		// 4EI/L
  
  kb(0,0) = EAoverL;
  kb(1,1) = kb(2,2) = EIoverL4;
  kb(2,1) = kb(1,2) = EIoverL2;
  
  return theCoordTransf->getInitialGlobalStiffMatrix(kb);
}

const Matrix &
ElasticBeam2d::getMass(void)
{ 
    K.Zero();
    
    if (rho > 0.0)  {
        // get initial element length
        double L = theCoordTransf->getInitialLength();
        if (cMass == 0)  {
            // lumped mass matrix
            double m = 0.5*rho*L;
            K(0,0) = K(1,1) = K(3,3) = K(4,4) = m;
        } else  {
            // consistent mass matrix
            static Matrix ml(6,6);
            double m = rho*L/420.0;
            ml(0,0) = ml(3,3) = m*140.0;
            ml(0,3) = ml(3,0) = m*70.0;

            ml(1,1) = ml(4,4) = m*156.0;
            ml(1,4) = ml(4,1) = m*54.0;
            ml(2,2) = ml(5,5) = m*4.0*L*L;
            ml(2,5) = ml(5,2) = -m*3.0*L*L;
            ml(1,2) = ml(2,1) = m*22.0*L;
            ml(4,5) = ml(5,4) = -ml(1,2);
            ml(1,5) = ml(5,1) = -m*13.0*L;
            ml(2,4) = ml(4,2) = -ml(1,5);
            
            // transform local mass matrix to global system
            K = theCoordTransf->getGlobalMatrixFromLocal(ml);
        }
    }
    
    return K;
}

void 
ElasticBeam2d::zeroLoad(void)
{
  Q.Zero();

  q0[0] = 0.0;
  q0[1] = 0.0;
  q0[2] = 0.0;

  p0[0] = 0.0;
  p0[1] = 0.0;
  p0[2] = 0.0;

  return;
}

int 
ElasticBeam2d::addLoad(ElementalLoad *theLoad, double loadFactor)
{
  int type;
  const Vector &data = theLoad->getData(type, loadFactor);
  double L = theCoordTransf->getInitialLength();

  if (type == LOAD_TAG_Beam2dUniformLoad) {
    double wt = data(0)*loadFactor;  // Transverse (+ve upward)
    double wa = data(1)*loadFactor;  // Axial (+ve from node I to J)

    double V = 0.5*wt*L;
    double M = V*L/6.0; // wt*L*L/12
    double P = wa*L;

    // Reactions in basic system
    p0[0] -= P;
    p0[1] -= V;
    p0[2] -= V;

    // Fixed end forces in basic system
    q0[0] -= 0.5*P;
    q0[1] -= M;
    q0[2] += M;
  }

  else if (type == LOAD_TAG_Beam2dPointLoad) {
    double P = data(0)*loadFactor;
    double N = data(1)*loadFactor;
    double aOverL = data(2);

    if (aOverL < 0.0 || aOverL > 1.0)
      return 0;

    double a = aOverL*L;
    double b = L-a;

    // Reactions in basic system
    p0[0] -= N;
    double V1 = P*(1.0-aOverL);
    double V2 = P*aOverL;
    p0[1] -= V1;
    p0[2] -= V2;

    double L2 = 1.0/(L*L);
    double a2 = a*a;
    double b2 = b*b;

    // Fixed end forces in basic system
    q0[0] -= N*aOverL;
    double M1 = -a * b2 * P * L2;
    double M2 = a2 * b * P * L2;
    q0[1] += M1;
    q0[2] += M2;
  }
  
  else if (type == LOAD_TAG_Beam2dTempLoad) {
    double Ttop1 = data(0)* loadFactor;
    double Tbot1 = data(1)* loadFactor;
    double Ttop2 = data(2)* loadFactor;
    double Tbot2 = data(3)* loadFactor;
        
    // fixed end forces due to a linear thermal load
    double dT1 = Ttop1-Tbot1;
    double dT = (Ttop2-Tbot2)-(Ttop1-Tbot1);
    double a = alpha/d;  // constant based on temp difference at top and bottom, 
    // coefficient of thermal expansion and beam depth
    double M1 = a*E*I*(-dT1+(4.0/3.0)*dT); //Fixed End Moment end 1
    double M2 = a*E*I*(dT1+(5.0/3.0)*dT); //Fixed End Moment end 2
    double F = alpha*(((Ttop2+Ttop1)/2+(Tbot2+Tbot1)/2)/2)*E*A; // Fixed End Axial Force
    double M1M2divL =(M1+M2)/L; // Fixed End Shear
    
    // Reactions in basic system
    p0[0] += 0;
    p0[1] += M1M2divL;
    p0[2] -= M1M2divL;

    // Fixed end forces in basic system
    q0[0] -= F;
    q0[1] += M1;
    q0[2] += M2;
  }

  else {
    opserr << "ElasticBeam2d::addLoad()  -- load type unknown for element with tag: " << this->getTag() << endln;
    return -1;
  }

  return 0;
}

int
ElasticBeam2d::addInertiaLoadToUnbalance(const Vector &accel)
{
  if (rho == 0.0)
    return 0;

  // get R * accel from the nodes
  const Vector &Raccel1 = theNodes[0]->getRV(accel);
  const Vector &Raccel2 = theNodes[1]->getRV(accel);
	
  if (3 != Raccel1.Size() || 3 != Raccel2.Size()) {
    opserr << "ElasticBeam2d::addInertiaLoadToUnbalance matrix and vector sizes are incompatable\n";
    return -1;
  }
    
  // want to add ( - fact * M R * accel ) to unbalance
  if (cMass == 0)  {
    // take advantage of lumped mass matrix
    double L = theCoordTransf->getInitialLength();
    double m = 0.5*rho*L;

    Q(0) -= m * Raccel1(0);
    Q(1) -= m * Raccel1(1);

    Q(3) -= m * Raccel2(0);
    Q(4) -= m * Raccel2(1);
  } else  {
    // use matrix vector multip. for consistent mass matrix
    static Vector Raccel(6);
    for (int i=0; i<3; i++)  {
      Raccel(i)   = Raccel1(i);
      Raccel(i+3) = Raccel2(i);
    }
    Q.addMatrixVector(1.0, this->getMass(), Raccel, -1.0);
  }
  
  return 0;
}

const Vector &
ElasticBeam2d::getResistingForceIncInertia()
{	
  P = this->getResistingForce();
  
  // subtract external load P = P - Q
  P.addVector(1.0, Q, -1.0);
  
  // add the damping forces if rayleigh damping
  if (alphaM != 0.0 || betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0)
    P.addVector(1.0, this->getRayleighDampingForces(), 1.0);
    
  if (rho == 0.0)
    return P;

  // add inertia forces from element mass
  const Vector &accel1 = theNodes[0]->getTrialAccel();
  const Vector &accel2 = theNodes[1]->getTrialAccel();    
  
  if (cMass == 0)  {
    // take advantage of lumped mass matrix
    double L = theCoordTransf->getInitialLength();
    double m = 0.5*rho*L;

    P(0) += m * accel1(0);
    P(1) += m * accel1(1);

    P(3) += m * accel2(0);
    P(4) += m * accel2(1);
  } else  {
    // use matrix vector multip. for consistent mass matrix
    static Vector accel(6);
    for (int i=0; i<3; i++)  {
      accel(i)   = accel1(i);
      accel(i+3) = accel2(i);
    }
    P.addMatrixVector(1.0, this->getMass(), accel, 1.0);
  }
  
  return P;
}


const Vector &
ElasticBeam2d::getResistingForce()
{
  theCoordTransf->update();
  
  const Vector &v = theCoordTransf->getBasicTrialDisp();
  double L = theCoordTransf->getInitialLength();

  double EoverL   = E/L;
  double EAoverL  = A*EoverL;			// EA/L
  double EIoverL2 = 2.0*I*EoverL;		// 2EI/L
  double EIoverL4 = 2.0*EIoverL2;		// 4EI/L
  
  // determine q = kv + q0
  q(0) = EAoverL*v(0);
  q(1) = EIoverL4*v(1) + EIoverL2*v(2);
  q(2) = EIoverL2*v(1) + EIoverL4*v(2);

  q(0) += q0[0];
  q(1) += q0[1];
  q(2) += q0[2];
  
  // Vector for reactions in basic system
  Vector p0Vec(p0, 3);
  
  P = theCoordTransf->getGlobalResistingForce(q, p0Vec);

  return P;
}

int
ElasticBeam2d::sendSelf(int cTag, Channel &theChannel)
{
  int res = 0;

    static Vector data(16);
    
    data(0) = A;
    data(1) = E; 
    data(2) = I; 
    data(3) = rho;
    data(4) = cMass;
    data(5) = this->getTag();
    data(6) = connectedExternalNodes(0);
    data(7) = connectedExternalNodes(1);
    data(8) = theCoordTransf->getClassTag();
    	
    int dbTag = theCoordTransf->getDbTag();
    
    if (dbTag == 0) {
      dbTag = theChannel.getDbTag();
      if (dbTag != 0)
	theCoordTransf->setDbTag(dbTag);
    }

    data(9) = dbTag;
    data(10) = alpha;
    data(11) = d;

    data(12) = alphaM;
    data(13) = betaK;
    data(14) = betaK0;
    data(15) = betaKc;

    // Send the data vector
    res += theChannel.sendVector(this->getDbTag(), cTag, data);
    if (res < 0) {
      opserr << "ElasticBeam2d::sendSelf -- could not send data Vector\n";
      return res;
    }

    // Ask the CoordTransf to send itself
    res += theCoordTransf->sendSelf(cTag, theChannel);
    if (res < 0) {
      opserr << "ElasticBeam2d::sendSelf -- could not send CoordTransf\n";
      return res;
    }
    
    return res;
}

int
ElasticBeam2d::recvSelf(int cTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
    int res = 0;
	
    static Vector data(16);

    res += theChannel.recvVector(this->getDbTag(), cTag, data);
    if (res < 0) {
      opserr << "ElasticBeam2d::recvSelf -- could not receive data Vector\n";
      return res;
    }

    A = data(0);
    E = data(1); 
    I = data(2); 
    alpha = data(10);
    d = data(11);

    alphaM = data(12);
    betaK  = data(13);
    betaK0 = data(14);
    betaKc = data(15);

    rho = data(3);
    cMass = (int)data(4);
    this->setTag((int)data(5));
    connectedExternalNodes(0) = (int)data(6);
    connectedExternalNodes(1) = (int)data(7);

    // Check if the CoordTransf is null; if so, get a new one
    int crdTag = (int)data(8);
    if (theCoordTransf == 0) {
      theCoordTransf = theBroker.getNewCrdTransf(crdTag);
      if (theCoordTransf == 0) {
	opserr << "ElasticBeam2d::recvSelf -- could not get a CrdTransf2d\n";
	exit(-1);
      }
    }
    
    // Check that the CoordTransf is of the right type; if not, delete
    // the current one and get a new one of the right type
    if (theCoordTransf->getClassTag() != crdTag) {
      delete theCoordTransf;
      theCoordTransf = theBroker.getNewCrdTransf(crdTag);
      if (theCoordTransf == 0) {
	opserr << "ElasticBeam2d::recvSelf -- could not get a CrdTransf2d\n";
	exit(-1);
      }
    }
	
    // Now, receive the CoordTransf
    theCoordTransf->setDbTag((int)data(9));
    res += theCoordTransf->recvSelf(cTag, theChannel, theBroker);
    if (res < 0) {
      opserr << "ElasticBeam2d::recvSelf -- could not receive CoordTransf\n";
      return res;
    }
    
    return res;
}

void
ElasticBeam2d::Print(OPS_Stream &s, int flag)
{
  // to update forces!
  this->getResistingForce();

  if (flag == -1) {
    int eleTag = this->getTag();
    s << "EL_BEAM\t" << eleTag << "\t";
    s << 0 << "\t" << 0 << "\t" << connectedExternalNodes(0) << "\t" << connectedExternalNodes(1) ;
    s << "0\t0.0000000\n";
  }

  if (flag == OPS_PRINT_CURRENTSTATE) {
    this->getResistingForce();
    s << "\nElasticBeam2d: " << this->getTag() << endln;
    s << "\tConnected Nodes: " << connectedExternalNodes ;
    s << "\tCoordTransf: " << theCoordTransf->getTag() << endln;
    s << "\tmass density:  " << rho << ", cMass: " << cMass << endln;
    double P  = q(0);
    double M1 = q(1);
    double M2 = q(2);
    double L = theCoordTransf->getInitialLength();
    double V = (M1+M2)/L;
    s << "\tEnd 1 Forces (P V M): " << -P+p0[0]
      << " " << V+p0[1] << " " << M1 << endln;
    s << "\tEnd 2 Forces (P V M): " << P
      << " " << -V+p0[2] << " " << M2 << endln;
  }

  if (flag == OPS_PRINT_PRINTMODEL_JSON) {
    s << "\t\t\t{";
	s << "\"name\": " << this->getTag() << ", ";
	s << "\"type\": \"ElasticBeam2d\", ";
    s << "\"nodes\": [" << connectedExternalNodes(0) << ", " << connectedExternalNodes(1) << "], ";
	s << "\"E\": " << E << ", ";
	s << "\"A\": "<< A << ", ";
    s << "\"Iz\": "<< I << ", ";
    s << "\"massperlength\": "<< rho << ", ";
    s << "\"crdTransformation\": \"" << theCoordTransf->getTag() << "\"}";
  }
}

int
ElasticBeam2d::displaySelf(Renderer &theViewer, int displayMode, float fact, const char **modes, int numMode)
{
  static Vector v1(3);
  static Vector v2(3);
  static Vector vp(3);

  theNodes[0]->getDisplayCrds(v1, fact);
  theNodes[1]->getDisplayCrds(v2, fact);

  float d1 = 0.0;
  float d2 = 0.0;
  float d3 = 0.0;

  int res = 0;

  if (displayMode > 0 && numMode == 0) {

    res += theViewer.drawLine(v1, v2, d1, d1, this->getTag(), 0);
    
  } else if (displayMode < 0) {
    
    theNodes[0]->getDisplayCrds(v1, 0.);
    theNodes[1]->getDisplayCrds(v2, 0.);
    
    // add eigenvector values
    int mode = displayMode  *  -1;

    const Matrix &eigen1 = theNodes[0]->getEigenvectors();
    const Matrix &eigen2 = theNodes[1]->getEigenvectors();
    if (eigen1.noCols() >= mode) {
      for (int i = 0; i < 2; i++) {
	v1(i) += eigen1(i,mode-1)*fact;
	v2(i) += eigen2(i,mode-1)*fact;    
      }    
    }

    res = theViewer.drawLine (v1, v2, 0.0, 0.0, this->getTag(), 0);
  }

  if (numMode > 0) {
    // calculate q for potential need below
    this->getResistingForce();
    vp = theCoordTransf->getBasicTrialDisp();
  }
  
  for (int i=0; i<numMode; i++) {

    const char *theMode = modes[i];
    if (strcmp(theMode, "axialForce") == 0) {
      d1 = q(0); 
      res +=theViewer.drawLine(v1, v2, d1, d1, this->getTag(), i);
      
    } else if (strcmp(theMode, "endMoments") == 0) {

      d1 = q(1);
      d2 = q(2);
      static Vector delta(3); delta = v2-v1; delta/=20.;
      res += theViewer.drawPoint(v1+delta, d1, this->getTag(), i);
      res += theViewer.drawPoint(v2-delta, d2, this->getTag(), i);

    } else if (strcmp(theMode, "localForces") == 0) {
      d1 = q(0);
      d2 = q(1);
      d3 = q(2);
      static Vector delta(3); delta = v2-v1; delta/=20;
      res += theViewer.drawPoint(v1+delta, d2, this->getTag(), i);
      res += theViewer.drawPoint(v2-delta, d3, this->getTag(), i);
      res +=theViewer.drawLine(v1, v2, d1, d1, this->getTag(), i);

    } else if (strcmp(theMode, "axialDeformation") == 0) {
      d1 = vp(0); 
      res +=theViewer.drawLine(v1, v2, d1, d1, this->getTag(), i);
      
    } else if (strcmp(theMode, "endRotations") == 0) {

      d1 = vp(1);
      d2 = vp(2);
      static Vector delta(3); delta = v2-v1; delta/=20.;
      res += theViewer.drawPoint(v1+delta, d1, this->getTag(), i);
      res += theViewer.drawPoint(v2-delta, d2, this->getTag(), i);

    } else if (strcmp(theMode, "localDeformations") == 0) {
      d1 = vp(0);
      d2 = vp(1);
      d3 = vp(2);
      static Vector delta(3); delta = v2-v1; delta/=20;
      res += theViewer.drawPoint(v1+delta, d2, this->getTag(), i);
      res += theViewer.drawPoint(v2-delta, d3, this->getTag(), i);
      res +=theViewer.drawLine(v1, v2, d1, d1, this->getTag(), i);

    } else if (strcmp(theMode, "plasticDeformations") == 0) {
      d1 = 0.;
      d2 = 0.;
      d3 = 0.;
      static Vector delta(3); delta = v2-v1; delta/=20;
      res += theViewer.drawPoint(v1+delta, d2, this->getTag(), i);
      res += theViewer.drawPoint(v2-delta, d3, this->getTag(), i);
      res +=theViewer.drawLine(v1, v2, d1, d1, this->getTag(), i);
    }

  }    

  return res;
}

Response*
ElasticBeam2d::setResponse(const char **argv, int argc, OPS_Stream &output)
{

  Response *theResponse = 0;

  output.tag("ElementOutput");
  output.attr("eleType","ElasticBeam2d");
  output.attr("eleTag",this->getTag());
  output.attr("node1",connectedExternalNodes[0]);
  output.attr("node2",connectedExternalNodes[1]);

    // global forces
  if (strcmp(argv[0],"force") == 0 || strcmp(argv[0],"forces") == 0 ||
      strcmp(argv[0],"globalForce") == 0 || strcmp(argv[0],"globalForces") == 0) {

    output.tag("ResponseType","Px_1");
    output.tag("ResponseType","Py_1");
    output.tag("ResponseType","Mz_1");
    output.tag("ResponseType","Px_2");
    output.tag("ResponseType","Py_2");
    output.tag("ResponseType","Mz_2");

    theResponse =  new ElementResponse(this, 2, P);
  
  // local forces
  }    else if (strcmp(argv[0],"localForce") == 0 || strcmp(argv[0],"localForces") == 0) {

    output.tag("ResponseType","N_1");
    output.tag("ResponseType","V_1");
    output.tag("ResponseType","M_1");
    output.tag("ResponseType","N_2");
    output.tag("ResponseType","V_2");
    output.tag("ResponseType","M_2");
    
    theResponse = new ElementResponse(this, 3, P);

  // basic forces
  }    else if (strcmp(argv[0],"basicForce") == 0 || strcmp(argv[0],"basicForces") == 0) {

    output.tag("ResponseType","N");
    output.tag("ResponseType","M_1");
    output.tag("ResponseType","M_2");
    
    theResponse = new ElementResponse(this, 4, Vector(3));

    // deformations
  }  else if (strcmp(argv[0],"deformatons") == 0 || 
	      strcmp(argv[0],"basicDeformations") == 0) {
    
    output.tag("ResponseType","eps");
    output.tag("ResponseType","theta1");
    output.tag("ResponseType","theta2");
    theResponse = new ElementResponse(this, 5, Vector(3));
  
  // chord rotation -
  } else if (strcmp(argv[0],"chordRotation") == 0 || strcmp(argv[0],"chordDeformation") == 0 
	     || strcmp(argv[0],"basicDeformation") == 0) {

    output.tag("ResponseType","eps");
    output.tag("ResponseType","theta1");
    output.tag("ResponseType","theta2");

    theResponse =  new ElementResponse(this, 5, Vector(3));
  }
  output.endTag(); // ElementOutput
  
  return theResponse;
}

int
ElasticBeam2d::getResponse (int responseID, Information &eleInfo)
{
  double N, M1, M2, V;
  double L = theCoordTransf->getInitialLength();
  this->getResistingForce();

  switch (responseID) {
  case 1: // stiffness
    return eleInfo.setMatrix(this->getTangentStiff());
    
  case 2: // global forces
    return eleInfo.setVector(this->getResistingForce());
    
  case 3: // local forces
    // Axial
    N = q(0);
    P(3) =  N;
    P(0) = -N+p0[0];
    // Moment
    M1 = q(1);
    M2 = q(2);
    P(2) = M1;
    P(5) = M2;
    // Shear
    V = (M1+M2)/L;
    P(1) =  V+p0[1];
    P(4) = -V+p0[2];
    return eleInfo.setVector(P);
    
  case 4: // basic forces
    return eleInfo.setVector(q);

  case 5:
    return eleInfo.setVector(theCoordTransf->getBasicTrialDisp());

  default:
    return -1;
  }
}

int
ElasticBeam2d::setParameter(const char **argv, int argc, Parameter &param)
{
  if (argc < 1)
    return -1;

  // E of the beam interior
  if (strcmp(argv[0],"E") == 0)
    return param.addObject(1, this);

  // A of the beam interior
  if (strcmp(argv[0],"A") == 0)
    return param.addObject(2, this);
  
  // I of the beam interior
  if (strcmp(argv[0],"I") == 0)
    return param.addObject(3, this);
  
  return -1;
}

int
ElasticBeam2d::updateParameter (int parameterID, Information &info)
{
	switch (parameterID) {
	case -1:
		return -1;
	case 1:
		E = info.theDouble;
		return 0;
	case 2:
		A = info.theDouble;
		return 0;
	case 3:
		I = info.theDouble;
		return 0;
	default:
		return -1;
	}
}

