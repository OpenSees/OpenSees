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
                                                                        
// $Revision: 1.13 $
// $Date: 2003-05-08 20:17:54 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/elasticBeamColumn/ElasticBeam3d.cpp,v $
                                                                        
                                                                        
// File: ~/model/ElasticBeam3d.C
//
// Written: fmk 11/95
// Revised:
//
// Purpose: This file contains the class definition for ElasticBeam3d.
// ElasticBeam3d is a 3d beam element. As such it can only
// connect to a node with 6-dof. 

#include <ElasticBeam3d.h>
#include <Domain.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>

#include <CrdTransf3d.h>
#include <Information.h>
#include <ElementResponse.h>
#include <ElementalLoad.h>
#include <Renderer.h>
#include <SectionForceDeformation.h>
#include <ID.h>
#include <math.h>
#include <stdlib.h>

Matrix ElasticBeam3d::K(12,12);
Vector ElasticBeam3d::P(12);
Matrix ElasticBeam3d::kb(6,6);

ElasticBeam3d::ElasticBeam3d()
  :Element(0,ELE_TAG_ElasticBeam3d), 
  A(0), E(0), G(0), Jx(0), Iy(0), Iz(0), rho(0.0),
  Q(12), q(6), connectedExternalNodes(2), theCoordTransf(0)
{
  // does nothing
  q0[0] = 0.0;
  q0[1] = 0.0;
  q0[2] = 0.0;
  q0[3] = 0.0;
  q0[4] = 0.0;

  p0[0] = 0.0;
  p0[1] = 0.0;
  p0[2] = 0.0;
  p0[3] = 0.0;
  p0[4] = 0.0;
  
  // set node pointers to NULL
  for (int i=0; i<2; i++)
    theNodes[i] = 0;      
}

ElasticBeam3d::ElasticBeam3d(int tag, double a, double e, double g, 
			     double jx, double iy, double iz, int Nd1, int Nd2, 
			     CrdTransf3d &coordTransf, double r, int sectTag)
  :Element(tag,ELE_TAG_ElasticBeam3d), 
   A(a), E(e), G(g), Jx(jx), Iy(iy), Iz(iz), rho(r), sectionTag(sectTag),
  Q(12), q(6), connectedExternalNodes(2), theCoordTransf(0)
{
  connectedExternalNodes(0) = Nd1;
  connectedExternalNodes(1) = Nd2;
  
  theCoordTransf = coordTransf.getCopy();
  
  if (!theCoordTransf) {
    opserr << "ElasticBeam3d::ElasticBeam3d -- failed to get copy of coordinate transformation\n";
    exit(-1);
  }

  q0[0] = 0.0;
  q0[1] = 0.0;
  q0[2] = 0.0;
  q0[3] = 0.0;
  q0[4] = 0.0;

  p0[0] = 0.0;
  p0[1] = 0.0;
  p0[2] = 0.0;
  p0[3] = 0.0;
  p0[4] = 0.0;

  // set node pointers to NULL
  for (int i=0; i<2; i++)
    theNodes[i] = 0;      
}

ElasticBeam3d::ElasticBeam3d(int tag, int Nd1, int Nd2, SectionForceDeformation *section, 			     
			     CrdTransf3d &coordTransf, double r)
  :Element(tag,ELE_TAG_ElasticBeam3d), 
  Q(12), q(6), connectedExternalNodes(2), theCoordTransf(0)
{
  if (section != 0) {
    sectionTag = section->getTag();
    E = 1.0;
    G = 1.0;
    Jx = 0.0;
    rho = r;

    const Matrix &sectTangent = section->getSectionTangent();
    const ID &sectCode = section->getType();
    for (int i=0; i<sectCode.Size(); i++) {
      int code = sectCode(i);
      switch(code) {
      case SECTION_RESPONSE_P:
	A = sectTangent(i,i);
	break;
      case SECTION_RESPONSE_MZ:
	Iz = sectTangent(i,i);
	break;
      case SECTION_RESPONSE_MY:
	Iy = sectTangent(i,i);
	break;
      case SECTION_RESPONSE_T:
	Jx = sectTangent(i,i);
	break;
      default:
	break;
      }
    }
  }    
  
  if (Jx == 0.0) {
    opserr << "ElasticBeam3d::ElasticBeam3d -- no torsion in section -- setting GJ = 1.0e10\n";
    Jx = 1.0e10;
  }

  connectedExternalNodes(0) = Nd1;
  connectedExternalNodes(1) = Nd2;
  
  theCoordTransf = coordTransf.getCopy();
  
  if (!theCoordTransf) {
    opserr << "ElasticBeam3d::ElasticBeam3d -- failed to get copy of coordinate transformation\n";
    exit(-1);
  }

  q0[0] = 0.0;
  q0[1] = 0.0;
  q0[2] = 0.0;
  q0[3] = 0.0;
  q0[4] = 0.0;

  p0[0] = 0.0;
  p0[1] = 0.0;
  p0[2] = 0.0;
  p0[3] = 0.0;
  p0[4] = 0.0;

  // set node pointers to NULL
  for (int i=0; i<2; i++)
    theNodes[i] = 0;      
}

ElasticBeam3d::~ElasticBeam3d()
{
  if (theCoordTransf)
    delete theCoordTransf;
}

int
ElasticBeam3d::getNumExternalNodes(void) const
{
    return 2;
}

const ID &
ElasticBeam3d::getExternalNodes(void) 
{
    return connectedExternalNodes;
}

Node **
ElasticBeam3d::getNodePtrs(void) 
{
  return theNodes;
}

int
ElasticBeam3d::getNumDOF(void)
{
    return 12;
}

void
ElasticBeam3d::setDomain(Domain *theDomain)
{
  if (theDomain == 0) {
    opserr << "ElasticBeam3d::setDomain -- Domain is null\n";
    exit(-1);
  }
    
    theNodes[0] = theDomain->getNode(connectedExternalNodes(0));
    theNodes[1] = theDomain->getNode(connectedExternalNodes(1));    
    

    if (theNodes[0] == 0) {
      opserr << "ElasticBeam3d::setDomain -- Node 1: " << connectedExternalNodes(0) << " does not exist\n";
      exit(-1);
    }
			      
    if (theNodes[1] == 0) {
      opserr << "ElasticBeam3d::setDomain -- Node 2: " << connectedExternalNodes(1) << " does not exist\n";
      exit(-1);
    }

    int dofNd1 = theNodes[0]->getNumberDOF();
    int dofNd2 = theNodes[1]->getNumberDOF();    
    
    if (dofNd1 != 6) {
      opserr << "ElasticBeam3d::setDomain -- Node 1: " << connectedExternalNodes(0) 
	     << " has incorrect number of DOF\n";
      exit(-1);
    }
    
    if (dofNd2 != 6) {
      opserr << "ElasticBeam3d::setDomain -- Node 2: " << connectedExternalNodes(1) 
	     << " has incorrect number of DOF\n";
      exit(-1);
    }
	
    this->DomainComponent::setDomain(theDomain);
    
    if (theCoordTransf->initialize(theNodes[0], theNodes[1]) != 0) {
	opserr << "ElasticBeam3d::setDomain -- Error initializing coordinate transformation\n";
	exit(-1);
    }
    
    double L = theCoordTransf->getInitialLength();

    if (L == 0.0) {
      opserr << "ElasticBeam3d::setDomain -- Element has zero length\n";
      exit(-1);
    }
}

int
ElasticBeam3d::commitState()
{
  int retVal = 0;
  // call element commitState to do any base class stuff
  if ((retVal = this->Element::commitState()) != 0) {
    opserr << "ElasticBeam3d::commitState () - failed in base class";
  }    
  retVal += theCoordTransf->commitState();
  return retVal;
}

int
ElasticBeam3d::revertToLastCommit()
{
    return theCoordTransf->revertToLastCommit();
}

int
ElasticBeam3d::revertToStart()
{
    return theCoordTransf->revertToStart();
}

int
ElasticBeam3d::update(void)
{
  return theCoordTransf->update();
}

const Matrix &
ElasticBeam3d::getTangentStiff(void)
{
  const Vector &v = theCoordTransf->getBasicTrialDisp();
  
  double L = theCoordTransf->getInitialLength();
  double oneOverL = 1.0/L;
  double EoverL   = E*oneOverL;
  double EAoverL  = A*EoverL;			// EA/L
  double EIzoverL2 = 2.0*Iz*EoverL;		// 2EIz/L
  double EIzoverL4 = 2.0*EIzoverL2;		// 4EIz/L
  double EIyoverL2 = 2.0*Iy*EoverL;		// 2EIy/L
  double EIyoverL4 = 2.0*EIyoverL2;		// 4EIy/L
  double GJoverL = G*Jx*oneOverL;         // GJ/L
  
  q(0) = EAoverL*v(0);
  q(1) = EIzoverL4*v(1) + EIzoverL2*v(2);
  q(2) = EIzoverL2*v(1) + EIzoverL4*v(2);
  q(3) = EIyoverL4*v(3) + EIyoverL2*v(4);
  q(4) = EIyoverL2*v(3) + EIyoverL4*v(4);    
  q(5) = GJoverL*v(5);

  q(0) += q0[0];
  q(1) += q0[1];
  q(2) += q0[2];
  q(3) += q0[3];
  q(4) += q0[4];
  
  kb(0,0) = EAoverL;
  kb(1,1) = kb(2,2) = EIzoverL4;
  kb(2,1) = kb(1,2) = EIzoverL2;
  kb(3,3) = kb(4,4) = EIyoverL4;
  kb(4,3) = kb(3,4) = EIyoverL2;
  kb(5,5) = GJoverL;
  
  return theCoordTransf->getGlobalStiffMatrix(kb,q);
}


const Matrix &
ElasticBeam3d::getInitialStiff(void)
{
  //  const Vector &v = theCoordTransf->getBasicTrialDisp();
  
  double L = theCoordTransf->getInitialLength();
  double oneOverL = 1.0/L;
  double EoverL   = E*oneOverL;
  double EAoverL  = A*EoverL;			// EA/L
  double EIzoverL2 = 2.0*Iz*EoverL;		// 2EIz/L
  double EIzoverL4 = 2.0*EIzoverL2;		// 4EIz/L
  double EIyoverL2 = 2.0*Iy*EoverL;		// 2EIy/L
  double EIyoverL4 = 2.0*EIyoverL2;		// 4EIy/L
  double GJoverL = G*Jx*oneOverL;         // GJ/L
  
  kb(0,0) = EAoverL;
  kb(1,1) = kb(2,2) = EIzoverL4;
  kb(2,1) = kb(1,2) = EIzoverL2;
  kb(3,3) = kb(4,4) = EIyoverL4;
  kb(4,3) = kb(3,4) = EIyoverL2;
  kb(5,5) = GJoverL;
  
  return theCoordTransf->getInitialGlobalStiffMatrix(kb);
}

const Matrix &
ElasticBeam3d::getMass(void)
{ 
  K.Zero();

  if (rho > 0.0) {
    double L = theCoordTransf->getInitialLength();
    double m = 0.5*rho*L;
    
    K(0,0) = m;
    K(1,1) = m;
    K(2,2) = m;
    
    K(6,6) = m;
    K(7,7) = m;
    K(8,8) = m;
  }
  
  return K;
}

void 
ElasticBeam3d::zeroLoad(void)
{
  Q.Zero();

  q0[0] = 0.0;
  q0[1] = 0.0;
  q0[2] = 0.0;
  q0[3] = 0.0;
  q0[4] = 0.0;

  p0[0] = 0.0;
  p0[1] = 0.0;
  p0[2] = 0.0;
  p0[3] = 0.0;
  p0[4] = 0.0;

  return;
}

int 
ElasticBeam3d::addLoad(ElementalLoad *theLoad, double loadFactor)
{
  int type;
  const Vector &data = theLoad->getData(type, loadFactor);
  double L = theCoordTransf->getInitialLength();

  if (type == LOAD_TAG_Beam3dUniformLoad) {
    double wy = data(0)*loadFactor;  // Transverse
    double wz = data(1)*loadFactor;  // Transverse
    double wx = data(2)*loadFactor;  // Axial (+ve from node I to J)

    double Vy = 0.5*wy*L;
    double Mz = Vy*L/6.0; // wy*L*L/12
    double Vz = 0.5*wz*L;
    double My = Vz*L/6.0; // wz*L*L/12
    double P = wx*L;

    // Reactions in basic system
    p0[0] -= P;
    p0[1] -= Vy;
    p0[2] -= Vy;
    p0[3] -= Vz;
    p0[4] -= Vz;

    // Fixed end forces in basic system
    q0[0] -= 0.5*P;
    q0[1] -= Mz;
    q0[2] += Mz;
    q0[3] += My;
    q0[4] -= My;
  }
  else if (type == LOAD_TAG_Beam3dPointLoad) {
    double Py = data(0)*loadFactor;
    double Pz = data(1)*loadFactor;
    double N  = data(2)*loadFactor;
    double aOverL = data(3);
    double a = aOverL*L;
    double b = L-a;

    // Reactions in basic system
    p0[0] -= N;
    double V1, V2;
    V1 = Py*(1.0-aOverL);
    V2 = Py*aOverL;
    p0[1] -= V1;
    p0[2] -= V2;
    V1 = Pz*(1.0-aOverL);
    V2 = Pz*aOverL;
    p0[3] -= V1;
    p0[4] -= V2;

    double L2 = 1.0/(L*L);
    double a2 = a*a;
    double b2 = b*b;

    // Fixed end forces in basic system
    q0[0] -= N*aOverL;
    double M1, M2;
    M1 = -a * b2 * Py * L2;
    M2 = a2 * b * Py * L2;
    q0[1] += M1;
    q0[2] += M2;
    M1 = -a * b2 * Pz * L2;
    M2 = a2 * b * Pz * L2;
    q0[3] -= M1;
    q0[4] -= M2;
  }
  else {
    opserr << "ElasticBeam3d::addLoad()  -- load type unknown for element with tag: " << this->getTag() << endln;
    return -1;
  }

  return 0;
}


int
ElasticBeam3d::addInertiaLoadToUnbalance(const Vector &accel)
{
  if (rho == 0.0)
    return 0;

  // Get R * accel from the nodes
  const Vector &Raccel1 = theNodes[0]->getRV(accel);
  const Vector &Raccel2 = theNodes[1]->getRV(accel);
	
  if (6 != Raccel1.Size() || 6 != Raccel2.Size()) {
    opserr << "ElasticBeam3d::addInertiaLoadToUnbalance matrix and vector sizes are incompatable\n";
    return -1;
  }

  // Want to add ( - fact * M R * accel ) to unbalance
  // Take advantage of lumped mass matrix
  double L = theCoordTransf->getInitialLength();
  double m = 0.5*rho*L;
  
  Q(0) -= m * Raccel1(0);
  Q(1) -= m * Raccel1(1);
  Q(2) -= m * Raccel1(2);
    
  Q(6) -= m * Raccel2(0);    
  Q(7) -= m * Raccel2(1);
  Q(8) -= m * Raccel2(2);    
  
  return 0;
}



const Vector &
ElasticBeam3d::getResistingForceIncInertia()
{	
  P = this->getResistingForce();

  // add the damping forces if rayleigh damping
  if (alphaM != 0.0 || betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0)
    P += this->getRayleighDampingForces();
    
  if (rho == 0.0)
    return P;

  else{
    const Vector &accel1 = theNodes[0]->getTrialAccel();
    const Vector &accel2 = theNodes[1]->getTrialAccel();    
    
    double L = theCoordTransf->getInitialLength();
    double m = 0.5*rho*L;
    
    P(0) += m * accel1(0);
    P(1) += m * accel1(1);
    P(2) += m * accel1(2);
    
    P(6) += m * accel2(0);    
    P(7) += m * accel2(1);
    P(8) += m * accel2(2);    
    
    return P;
  }
}


const Vector &
ElasticBeam3d::getResistingForce()
{
  const Vector &v = theCoordTransf->getBasicTrialDisp();
  
  double L = theCoordTransf->getInitialLength();
  double oneOverL = 1.0/L;
  double EoverL   = E*oneOverL;
  double EAoverL  = A*EoverL;			// EA/L
  double EIzoverL2 = 2.0*Iz*EoverL;		// 2EIz/L
  double EIzoverL4 = 2.0*EIzoverL2;		// 4EIz/L
  double EIyoverL2 = 2.0*Iy*EoverL;		// 2EIy/L
  double EIyoverL4 = 2.0*EIyoverL2;		// 4EIy/L
  double GJoverL = G*Jx*oneOverL;         // GJ/L
  
  q(0) = EAoverL*v(0);
  q(1) = EIzoverL4*v(1) + EIzoverL2*v(2);
  q(2) = EIzoverL2*v(1) + EIzoverL4*v(2);
  q(3) = EIyoverL4*v(3) + EIyoverL2*v(4);
  q(4) = EIyoverL2*v(3) + EIyoverL4*v(4);    
  q(5) = GJoverL*v(5);
  
  q(0) += q0[0];
  q(1) += q0[1];
  q(2) += q0[2];
  q(3) += q0[3];
  q(4) += q0[4];
  
  Vector p0Vec(p0, 5);
  
  //  opserr << q;

  P = theCoordTransf->getGlobalResistingForce(q, p0Vec);

  // opserr << P;
  
  // P = P - Q;
  P.addVector(1.0, Q, -1.0);
  
  return P;
}

int
ElasticBeam3d::sendSelf(int cTag, Channel &theChannel)
{
    int res = 0;

    static Vector data(12);
    
    data(0) = A;
    data(1) = E; 
    data(2) = G; 
    data(3) = Jx; 
    data(4) = Iy; 
    data(5) = Iz;     
    data(6) = rho;
    data(7) = this->getTag();
    data(8) = connectedExternalNodes(0);
    data(9) = connectedExternalNodes(1);
    data(10) = theCoordTransf->getClassTag();    	
	
	int dbTag = theCoordTransf->getDbTag();

	if (dbTag == 0) {
		dbTag = theChannel.getDbTag();
		if (dbTag != 0)
			theCoordTransf->setDbTag(dbTag);
	}

	data(11) = dbTag;

	// Send the data vector
    res += theChannel.sendVector(this->getDbTag(), cTag, data);
    if (res < 0) {
      opserr << "ElasticBeam3d::sendSelf -- could not send data Vector\n";
      return res;
    }

    // Ask the CoordTransf to send itself
    res += theCoordTransf->sendSelf(cTag, theChannel);
    if (res < 0) {
      opserr << "ElasticBeam3d::sendSelf -- could not send CoordTransf\n";
      return res;
    }

    return res;
}

int
ElasticBeam3d::recvSelf(int cTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
	int res = 0;
	
	static Vector data(12);

    res += theChannel.recvVector(this->getDbTag(), cTag, data);
    if (res < 0) {
      opserr << "ElasticBeam3d::recvSelf -- could not receive data Vector\n";
      return res;
    }

    A = data(0);
    E = data(1); 
    G = data(2); 
    Jx = data(3); 
    Iy = data(4); 
    Iz = data(5);     
    rho = data(6);
    this->setTag((int)data(7));
    connectedExternalNodes(0) = (int)data(8);
    connectedExternalNodes(1) = (int)data(9);

    // Check if the CoordTransf is null; if so, get a new one
    int crdTag = (int)data(10);
    if (theCoordTransf == 0) {
      theCoordTransf = theBroker.getNewCrdTransf3d(crdTag);
      if (theCoordTransf == 0) {
	opserr << "ElasticBeam3d::recvSelf -- could not get a CrdTransf3d\n";
	exit(-1);
      }
    }

    // Check that the CoordTransf is of the right type; if not, delete
    // the current one and get a new one of the right type
    if (theCoordTransf->getClassTag() != crdTag) {
      delete theCoordTransf;
      theCoordTransf = theBroker.getNewCrdTransf3d(crdTag);
      if (theCoordTransf == 0) {
	opserr << "ElasticBeam3d::recvSelf -- could not get a CrdTransf3d\n";
	exit(-1);
      }
    }

    // Now, receive the CoordTransf
    theCoordTransf->setDbTag((int)data(11));
    res += theCoordTransf->recvSelf(cTag, theChannel, theBroker);
    if (res < 0) {
      opserr << "ElasticBeam3d::recvSelf -- could not receive CoordTransf\n";
      return res;
    }

    // Revert the crdtrasf to its last committed state
    theCoordTransf->revertToLastCommit();

    return res;
}

void
ElasticBeam3d::Print(OPS_Stream &s, int flag)
{
   if (flag == -1) { 
    int eleTag = this->getTag();
    s << "EL_BEAM\t" << eleTag << "\t";
    s << sectionTag << "\t" << sectionTag; 
    s  << "\t" << connectedExternalNodes(0) << "\t" << connectedExternalNodes(1);
    s << "\t0\t0.0000000\n";
   }  else if (flag < -1) {
     int counter = (flag + 1) * -1;
     int eleTag = this->getTag();
     int i;
     const Vector &force = this->getResistingForce();
     s << "FORCE\t" << eleTag << "\t" << counter << "\t0";
     for (i=0; i<3; i++)
	s << "\t" << force(i);
     s << endln;
     s << "FORCE\t" << eleTag << "\t" << counter << "\t1";
     for (i=0; i<3; i++)
       s << "\t" << force(i+6);
     s << endln;
     s << "MOMENT\t" << eleTag << "\t" << counter << "\t0";
     for (i=3; i<6; i++)
       s << "\t" << force(i);
     s << endln;
     s << "MOMENT\t" << eleTag << "\t" << counter << "\t1";
     for (i=3; i<6; i++)
       s << "\t" << force(i+6);
     s << endln;
   }

   else if (flag == 2){

     static Vector xAxis(3);
     static Vector yAxis(3);
     static Vector zAxis(3);
     
     theCoordTransf->getLocalAxes(xAxis, yAxis, zAxis);
                        
     s << "#ElasticBeamColumn3D\n";
     s << "#LocalAxis " << xAxis(0) << " " << xAxis(1) << " " << xAxis(2);
     s << " " << yAxis(0) << " " << yAxis(1) << " " << yAxis(2) << " ";
     s << zAxis(0) << " " << zAxis(1) << " " << zAxis(2) << endln;

     const Vector &node1Crd = theNodes[0]->getCrds();
     const Vector &node2Crd = theNodes[1]->getCrds();	
     const Vector &node1Disp = theNodes[0]->getDisp();
     const Vector &node2Disp = theNodes[1]->getDisp();    
     
     s << "#NODE " << node1Crd(0) << " " << node1Crd(1) << " " << node1Crd(2)
       << " " << node1Disp(0) << " " << node1Disp(1) << " " << node1Disp(2)
       << " " << node1Disp(3) << " " << node1Disp(4) << " " << node1Disp(5) << endln;
     
     s << "#NODE " << node2Crd(0) << " " << node2Crd(1) << " " << node2Crd(2)
       << " " << node2Disp(0) << " " << node2Disp(1) << " " << node2Disp(2)
       << " " << node2Disp(3) << " " << node2Disp(4) << " " << node2Disp(5) << endln;

    double N, Mz1, Mz2, Vy, My1, My2, Vz, T;
    double L = theCoordTransf->getInitialLength();
    double oneOverL = 1.0/L;
    
    N   = q(0);
    Mz1 = q(1);
    Mz2 = q(2);
    Vy  = (Mz1+Mz2)*oneOverL;
    My1 = q(3);
    My2 = q(4);
    Vz  = -(My1+My2)*oneOverL;
    T   = q(5);
    
    s << "#END_FORCES " << -N+p0[0] << ' ' <<  Vy+p0[1] << ' ' << Vz+p0[3] << ' ' 
      << -T << ' ' << My1 << ' ' <<  Mz1 << endln;
    s << "#END_FORCES " <<  N << ' ' << -Vy+p0[2] << ' ' << -Vz+p0[4] << ' ' 
      << T << ' ' << My2 << ' ' << Mz2 << endln;
   }
   else {

    s << "\nElasticBeam3d: " << this->getTag() << endln;
    s << "\tConnected Nodes: " << connectedExternalNodes ;
    s << "\tCoordTransf: " << theCoordTransf->getTag() << endln;
    
    double N, Mz1, Mz2, Vy, My1, My2, Vz, T;
    double L = theCoordTransf->getInitialLength();
    double oneOverL = 1.0/L;
    
    N   = q(0);
    Mz1 = q(1);
    Mz2 = q(2);
    Vy  = (Mz1+Mz2)*oneOverL;
    My1 = q(3);
    My2 = q(4);
    Vz  = -(My1+My2)*oneOverL;
    T   = q(5);
    
    s << "\tEnd 1 Forces (P Mz Vy My Vz T): "
      << -N+p0[0] << ' ' << Mz1 << ' ' <<  Vy+p0[1] << ' ' << My1 << ' ' <<  Vz+p0[3] << ' ' << -T << endln;
    s << "\tEnd 2 Forces (P Mz Vy My Vz T): "
      <<  N << ' ' << Mz2 << ' ' << -Vy+p0[2] << ' ' << My2 << ' ' << -Vz+p0[4] << ' ' <<  T << endln;
  }
}

int
ElasticBeam3d::displaySelf(Renderer &theViewer, int displayMode, float fact)
{
    // first determine the end points of the quad based on
    // the display factor (a measure of the distorted image)
    const Vector &end1Crd = theNodes[0]->getCrds();
    const Vector &end2Crd = theNodes[1]->getCrds();	

    const Vector &end1Disp = theNodes[0]->getDisp();
    const Vector &end2Disp = theNodes[1]->getDisp();

    static Vector v1(3);
    static Vector v2(3);
    
    for (int i = 0; i < 3; i++) {
      v1(i) = end1Crd(i) + end1Disp(i)*fact;
      v2(i) = end2Crd(i) + end2Disp(i)*fact;    
    }
    
    return theViewer.drawLine (v1, v2, 1.0, 1.0);
}

Response*
ElasticBeam3d::setResponse(const char **argv, int argc, Information &info)
{
    // stiffness
    if (strcmp(argv[0],"stiffness") == 0)
		return new ElementResponse(this, 1, K);

    // global forces
    else if (strcmp(argv[0],"force") == 0 || strcmp(argv[0],"forces") == 0 ||
		strcmp(argv[0],"globalForce") == 0 || strcmp(argv[0],"globalForces") == 0)
		return new ElementResponse(this, 2, P);

	// local forces
    else if (strcmp(argv[0],"localForce") == 0 || strcmp(argv[0],"localForces") == 0)
		return new ElementResponse(this, 3, P);

    else
		return 0;
}

int
ElasticBeam3d::getResponse (int responseID, Information &eleInfo)
{
  double N, V, M1, M2, T;
  double L = theCoordTransf->getInitialLength();
  double oneOverL = 1.0/L;

  switch (responseID) {
  case 1: // stiffness
    return eleInfo.setMatrix(this->getTangentStiff());
    
  case 2: // global forces
    return eleInfo.setVector(this->getResistingForce());
    
  case 3: // local forces
    // Axial
    N = q(0);
    P(6) =  N;
    P(0) = -N+p0[0];
    
    // Torsion
    T = q(5);
    P(9) =  T;
    P(3) = -T;
    
    // Moments about z and shears along y
    M1 = q(1);
    M2 = q(2);
    P(5)  = M1;
    P(11) = M2;
    V = (M1+M2)*oneOverL;
    P(1) =  V+p0[1];
    P(7) = -V+p0[2];
    
    // Moments about y and shears along z
    M1 = q(3);
    M2 = q(4);
    P(4)  = M1;
    P(10) = M2;
    V = -(M1+M2)*oneOverL;
    P(2) = -V+p0[3];
    P(8) =  V+p0[4];
    
    return eleInfo.setVector(P);
    
  default:
    return -1;
  }
}
