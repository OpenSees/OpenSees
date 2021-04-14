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
                                                                        
// $Revision: 1.20 $
// $Date: 2008/09/23 22:50:33 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/elasticBeamColumn/ElasticBeamWarping3d.cpp,v $
                                                                        
                                                                        
// File: ~/model/ElasticBeamWarping3d.C
//
// Written: fmk 11/95
// Revised:
//
// Purpose: This file contains the class definition for ElasticBeamWarping3d.
// ElasticBeamWarping3d is a 3d beam element. As such it can only
// connect to a node with 6-dof. 
// Modified by Xi Zhang from University of Sydney, Australia (include warping degrees of freedom). Refer to 
// Formulation and Implementation of Three-dimensional Doubly Symmetric Beam-Column Analyses with Warping Effects in OpenSees
// Research Report R917, School of Civil Engineering, University of Sydney.

#include <ElasticBeamWarping3d.h>
#include <Domain.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>

#include <CrdTransf.h>
#include <Information.h>
#include <ElementResponse.h>
#include <ElementalLoad.h>
#include <Renderer.h>
#include <SectionForceDeformation.h>
#include <ID.h>
#include <math.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <elementAPI.h>

Matrix ElasticBeamWarping3d::K(14,14);
Vector ElasticBeamWarping3d::P(14);
Matrix ElasticBeamWarping3d::kb(9,9);
using std::string;
using namespace std;

void* OPS_ElasticBeamWarping3d(void)
{
    int numArgs = OPS_GetNumRemainingInputArgs();
    if(numArgs < 11 && numArgs != 6) {
	opserr<<"insufficient arguments:eleTag,iNode,jNode,<A,E,G,J,Iy,Iz>or<sectionTag>,transfTag,Cw\n";
	return 0;
    }

    int ndm = OPS_GetNDM();
    int ndf = OPS_GetNDF();
    if(ndm != 3 || ndf != 7) {
	opserr<<"ndm must be 3 and ndf must be 7\n";
	return 0;
    }

    // inputs: 
    int iData[3];
    int numData = 3;
    if(OPS_GetIntInput(&numData,&iData[0]) < 0) return 0;

    SectionForceDeformation* theSection = 0;
    CrdTransf* theTrans = 0;
    double data[6];
    int transfTag, secTag;
    
    if(numArgs == 6) {
	numData = 1;
	if(OPS_GetIntInput(&numData,&secTag) < 0) return 0;
	if(OPS_GetIntInput(&numData,&transfTag) < 0) return 0;

	theSection = OPS_getSectionForceDeformation(secTag);
	if(theSection == 0) {
	    opserr<<"no section is found\n";
	    return 0;
	}
	theTrans = OPS_getCrdTransf(transfTag);
	if(theTrans == 0) {
	    opserr<<"no CrdTransf is found\n";
	    return 0;
	}
    } else {
	numData = 6;
	if(OPS_GetDoubleInput(&numData,&data[0]) < 0) return 0;
	numData = 1;
	if(OPS_GetIntInput(&numData,&transfTag) < 0) return 0;
	theTrans = OPS_getCrdTransf(transfTag);
	if(theTrans == 0) {
	    opserr<<"no CrdTransf is found\n";
	    return 0;
	}
    }

    // Read Cw
    numData = 1;
    double Cw;
    if(OPS_GetDoubleInput(&numData,&Cw) < 0)
      return 0;    
    
    // options
    double mass = 0.0;
    int cMass = 0;
    while(OPS_GetNumRemainingInputArgs() > 0) {
	std::string theType = OPS_GetString();
	if (theType == "-mass") {
	    if(OPS_GetNumRemainingInputArgs() > 0) {
		if(OPS_GetDoubleInput(&numData,&mass) < 0) return 0;
	    }
	}
    }

    if (theSection != 0) {
      return new ElasticBeamWarping3d(iData[0],iData[1],iData[2],theSection,*theTrans,Cw,mass); 
    } else {
      return new ElasticBeamWarping3d(iData[0],data[0],data[1],data[2],data[3],data[4],
				      data[5],iData[1],iData[2],*theTrans, Cw, mass);
    }
}

ElasticBeamWarping3d::ElasticBeamWarping3d()
  :Element(0,ELE_TAG_ElasticBeamWarping3d), 
  A(0), E(0), G(0), Jx(0), Iy(0), Iz(0), Cw(0), rho(0.0),
  Q(14), q(9), connectedExternalNodes(2), theCoordTransf(0)
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

ElasticBeamWarping3d::ElasticBeamWarping3d(int tag, double a, double e, double g, 
			     double jx, double iy, double iz, int Nd1, int Nd2, 
			     CrdTransf &coordTransf,  double cw, double r)
  :Element(tag,ELE_TAG_ElasticBeamWarping3d), 
   A(a), E(e), G(g), Jx(jx), Iy(iy), Iz(iz), Cw(cw), rho(r),
  Q(14), q(9), connectedExternalNodes(2), theCoordTransf(0)
{
  connectedExternalNodes(0) = Nd1;
  connectedExternalNodes(1) = Nd2;
  
  theCoordTransf = coordTransf.getCopy3d();
  
  if (!theCoordTransf) {
    opserr << "ElasticBeamWarping3d::ElasticBeamWarping3d -- failed to get copy of coordinate transformation\n";
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

ElasticBeamWarping3d::ElasticBeamWarping3d(int tag, int Nd1, int Nd2, SectionForceDeformation *section, 			     
					   CrdTransf &coordTransf, double cw, double r)
  :Element(tag,ELE_TAG_ElasticBeamWarping3d), Cw(cw), rho(r),
  Q(14), q(9), connectedExternalNodes(2), theCoordTransf(0)
{
  if (section != 0) {
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
    opserr << "ElasticBeamWarping3d::ElasticBeamWarping3d -- no torsion in section -- setting GJ = 1.0e10\n";
    Jx = 1.0e10;
  }

  connectedExternalNodes(0) = Nd1;
  connectedExternalNodes(1) = Nd2;
  
  theCoordTransf = coordTransf.getCopy3d();
  
  if (!theCoordTransf) {
    opserr << "ElasticBeamWarping3d::ElasticBeamWarping3d -- failed to get copy of coordinate transformation\n";
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

ElasticBeamWarping3d::~ElasticBeamWarping3d()
{
  if (theCoordTransf)
    delete theCoordTransf;
}

int
ElasticBeamWarping3d::getNumExternalNodes(void) const
{
    return 2;
}

const ID &
ElasticBeamWarping3d::getExternalNodes(void) 
{
    return connectedExternalNodes;
}

Node **
ElasticBeamWarping3d::getNodePtrs(void) 
{
  return theNodes;
}

int
ElasticBeamWarping3d::getNumDOF(void)
{
    return 14;
}

void
ElasticBeamWarping3d::setDomain(Domain *theDomain)
{
  if (theDomain == 0) {
    opserr << "ElasticBeamWarping3d::setDomain -- Domain is null\n";
    exit(-1);
  }
    
    theNodes[0] = theDomain->getNode(connectedExternalNodes(0));
    theNodes[1] = theDomain->getNode(connectedExternalNodes(1));    
    

    if (theNodes[0] == 0) {
      opserr << "ElasticBeamWarping3d::setDomain -- Node 1: " << connectedExternalNodes(0) << " does not exist\n";
      exit(-1);
    }
			      
    if (theNodes[1] == 0) {
      opserr << "ElasticBeamWarping3d::setDomain -- Node 2: " << connectedExternalNodes(1) << " does not exist\n";
      exit(-1);
    }

    int dofNd1 = theNodes[0]->getNumberDOF();
    int dofNd2 = theNodes[1]->getNumberDOF();    
    
    //if (dofNd1 != 6) {
     // opserr << "ElasticBeamWarping3d::setDomain -- Node 1: " << connectedExternalNodes(0) 
	    // << " has incorrect number of DOF\n";
     // exit(-1);
   // }
    
    /*if (dofNd2 != 6) {
      opserr << "ElasticBeamWarping3d::setDomain -- Node 2: " << connectedExternalNodes(1) 
	     << " has incorrect number of DOF\n";
      exit(-1);
    }*/
	
    this->DomainComponent::setDomain(theDomain);
    
    if (theCoordTransf->initialize(theNodes[0], theNodes[1]) != 0) {
	opserr << "ElasticBeamWarping3d::setDomain -- Error initializing coordinate transformation\n";
	exit(-1);
    }
    
    double L = theCoordTransf->getInitialLength();

    if (L == 0.0) {
      opserr << "ElasticBeamWarping3d::setDomain -- Element has zero length\n";
      exit(-1);
    }
}

int
ElasticBeamWarping3d::commitState()
{
  int retVal = 0;
  // call element commitState to do any base class stuff
  if ((retVal = this->Element::commitState()) != 0) {
    opserr << "ElasticBeamWarping3d::commitState () - failed in base class";
  }    
  retVal += theCoordTransf->commitState();
  return retVal;
}

int
ElasticBeamWarping3d::revertToLastCommit()
{
    return theCoordTransf->revertToLastCommit();
}

int
ElasticBeamWarping3d::revertToStart()
{
    return theCoordTransf->revertToStart();
}

int
ElasticBeamWarping3d::update(void)
{
  return theCoordTransf->update();
}

const Matrix &
ElasticBeamWarping3d::getTangentStiff(void)
{
  const Vector &v = theCoordTransf->getBasicTrialDisp();
  /*for (int j=0; j<9; j++)
				{
					uuoutput << "v"<<"\t\t"<<j<<"\t\t"<<v(j)<<endln;
					
			}
 uuoutput<<"------------------------------------------------------------"<<endln;
 uuoutput<<"finish v"<<endln;*/
  //double Cw=134.46;
  double L = theCoordTransf->getInitialLength();
  double oneOverL = 1.0/L;
  double EoverL   = E*oneOverL;
  double EAoverL  = A*EoverL;			// EA/L
  double EIzoverL2 = 2.0*Iz*EoverL;		// 2EIz/L
  double EIzoverL4 = 2.0*EIzoverL2;		// 4EIz/L
  double EIyoverL2 = 2.0*Iy*EoverL;		// 2EIy/L
  double EIyoverL4 = 2.0*EIyoverL2;		// 4EIy/L
  double GJoverL = G*Jx*oneOverL;         // GJ/L
  double ECoverL3=E*Cw/L/L/L;
  double ECoverL2=E*Cw/L/L;
  double ECoverL=E*Cw/L;
  double GJover10=G*Jx/10.0;
  double GJL=G*Jx*L;
  
  q(0) = (12.0*ECoverL3+6.0/5.0*GJoverL)*(v(0)-v(4))+(GJover10+6.0*ECoverL2)*(v(3)+v(7));
  q(1) = EIzoverL4*v(1) + EIzoverL2*v(5);
  q(2) = EIyoverL4*v(2) + EIyoverL2*v(6);
  q(3) = (GJover10+6.0*ECoverL2)*(v(0)-v(4))+(4.0*ECoverL+2.0/15.0*GJL)*v(3)+(2.0*ECoverL-1.0/30.0*GJL)*v(7);
  q(4) = (12.*ECoverL3+6./5.*GJoverL)*(v(4)-v(0))-(GJover10+6.*ECoverL2)*(v(3)+v(7));
  q(5) = EIzoverL2*v(1) + EIzoverL4*v(5);
  q(6) = EIyoverL2*v(2) + EIyoverL4*v(6);
  q(7) = (GJover10+6.*ECoverL2)*(v(0)-v(4))+(2.*ECoverL-1./30.*GJL)*v(3)+(4.*ECoverL+2./15.*GJL)*v(7);
  q(8) = EAoverL*v(8);

  q(0) += q0[0];
  q(1) += q0[1];
  q(2) += q0[2];
  q(3) += q0[3];
  q(4) += q0[4];

  kb(0,0) = 12.*ECoverL3+6./5.*GJoverL;
  kb(0,3) = kb(3,0) = kb(0,7) = kb(7,0) = GJover10+6.*ECoverL2;
  kb(0,4) = kb(4,0) = -12.*ECoverL3-6./5.*GJoverL;
  kb(1,1) = kb(5,5) = EIzoverL4;
  kb(1,5) = kb(5,1) = EIzoverL2;
  kb(2,2) = kb(6,6) = EIyoverL4;
  kb(2,6) = kb(6,2) = EIyoverL2;
  kb(3,3) = 4.*ECoverL+2./15.*GJL;
  kb(3,4) = kb(4,3) = -kb(3,0);
  kb(3,7) = kb(7,3) = 2.*ECoverL-1./30.*GJL;
  kb(4,4) = -kb(4,0);
  kb(4,7) = kb(7,4) = -GJover10-6.*ECoverL2;
  kb(7,7) = 4.*ECoverL+2./15.*GJL;
  kb(8,8) = EAoverL;
 
  return theCoordTransf->getGlobalStiffMatrix(kb,q);
}


const Matrix &
ElasticBeamWarping3d::getInitialStiff(void)
{
  //  const Vector &v = theCoordTransf->getBasicTrialDisp();
  double L = theCoordTransf->getInitialLength();
  //double Cw=134.46;
  double oneOverL = 1.0/L;
  double EoverL   = E*oneOverL;
  double EAoverL  = A*EoverL;			// EA/L
  double EIzoverL2 = 2.0*Iz*EoverL;		// 2EIz/L
  double EIzoverL4 = 2.0*EIzoverL2;		// 4EIz/L
  double EIyoverL2 = 2.0*Iy*EoverL;		// 2EIy/L
  double EIyoverL4 = 2.0*EIyoverL2;		// 4EIy/L
  double GJoverL = G*Jx*oneOverL;         // GJ/L
  double ECoverL3=E*Cw/L/L/L;
  double ECoverL2=E*Cw/L/L;
  double ECoverL=E*Cw/L;
  double GJover10=G*Jx/10.;
  double GJL=G*Jx*L;
  

  kb(0,0) = 12.*ECoverL3+6./5.*GJoverL;
  kb(0,3) = kb(3,0) = kb(0,7) = kb(7,0) = GJover10+6.*ECoverL2;
  kb(0,4) = kb(4,0) = -12.*ECoverL3-6./5.*GJoverL;
  kb(1,1) = kb(5,5) = EIzoverL4;
  kb(1,5) = kb(5,1) = EIzoverL2;
  kb(2,2) = kb(6,6) = EIyoverL4;
  kb(2,6) = kb(6,2) = EIyoverL2;
  kb(3,3) = 4.*ECoverL+2./15.*GJL;
  kb(3,4) = kb(4,3) = -kb(3,0);
  kb(3,7) = kb(7,3) = 2.*ECoverL-1./30.*GJL;
  kb(4,4) = -kb(4,0);
  kb(4,7) = kb(7,4) = -GJover10-6.*ECoverL2;
  kb(7,7) = 4.*ECoverL+2./15.*GJL;
  kb(8,8) = EAoverL;

  return theCoordTransf->getInitialGlobalStiffMatrix(kb);
}

const Matrix &
ElasticBeamWarping3d::getMass(void)
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
ElasticBeamWarping3d::zeroLoad(void)
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
ElasticBeamWarping3d::addLoad(ElementalLoad *theLoad, double loadFactor)
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

    if (aOverL < 0.0 || aOverL > 1.0)
      return 0;

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
    opserr << "ElasticBeamWarping3d::addLoad()  -- load type unknown for element with tag: " << this->getTag() << endln;
    return -1;
  }

  return 0;
}


int
ElasticBeamWarping3d::addInertiaLoadToUnbalance(const Vector &accel)
{
  if (rho == 0.0)
    return 0;

  // Get R * accel from the nodes
  const Vector &Raccel1 = theNodes[0]->getRV(accel);
  const Vector &Raccel2 = theNodes[1]->getRV(accel);
	
  if (6 != Raccel1.Size() || 6 != Raccel2.Size()) {
    opserr << "ElasticBeamWarping3d::addInertiaLoadToUnbalance matrix and vector sizes are incompatable\n";
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
ElasticBeamWarping3d::getResistingForceIncInertia()
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
ElasticBeamWarping3d::getResistingForce()
{
  const Vector &v = theCoordTransf->getBasicTrialDisp();
  // opserr << "v" <<v<< endln;
  double L = theCoordTransf->getInitialLength();
  //double Cw=134.46;
  double oneOverL = 1.0/L;
  double EoverL   = E*oneOverL;
  double EAoverL  = A*EoverL;			// EA/L
  double EIzoverL2 = 2.0*Iz*EoverL;		// 2EIz/L
  double EIzoverL4 = 2.0*EIzoverL2;		// 4EIz/L
  double EIyoverL2 = 2.0*Iy*EoverL;		// 2EIy/L
  double EIyoverL4 = 2.0*EIyoverL2;		// 4EIy/L
  double GJoverL = G*Jx*oneOverL;         // GJ/L
   double ECoverL3=E*Cw/L/L/L;
  double ECoverL2=E*Cw/L/L;
  double ECoverL=E*Cw/L;
  double GJover10=G*Jx/10.;
  double GJL=G*Jx*L;

  q(0) = (12.0*ECoverL3+6.0/5.0*GJoverL)*(v(0)-v(4))+(GJover10+6.0*ECoverL2)*(v(3)+v(7));
  q(1) = EIzoverL4*v(1) + EIzoverL2*v(5);
  q(2) = EIyoverL4*v(2) + EIyoverL2*v(6);
  q(3) = (GJover10+6.0*ECoverL2)*(v(0)-v(4))+(4.0*ECoverL+2.0/15.0*GJL)*v(3)+(2.0*ECoverL-1.0/30.0*GJL)*v(7);
  q(4) = (12.*ECoverL3+6./5.*GJoverL)*(v(4)-v(0))-(GJover10+6.*ECoverL2)*(v(3)+v(7));
  q(5) = EIzoverL2*v(1) + EIzoverL4*v(5);
  q(6) = EIyoverL2*v(2) + EIyoverL4*v(6);
  q(7) = (GJover10+6.*ECoverL2)*(v(0)-v(4))+(2.*ECoverL-1./30.*GJL)*v(3)+(4.*ECoverL+2./15.*GJL)*v(7);
  q(8) = EAoverL*v(8);

  
  q(0) += q0[0];
  q(1) += q0[1];
  q(2) += q0[2];
  q(3) += q0[3];
  q(4) += q0[4];

  Vector p0Vec(p0, 5);

  P = theCoordTransf->getGlobalResistingForce(q, p0Vec); 

  // P = P - Q;
  P.addVector(1.0, Q, -1.0);

  return P;
}

int
ElasticBeamWarping3d::sendSelf(int cTag, Channel &theChannel)
{
    int res = 0;

    static Vector data(16);
    
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
    
    data(12) = alphaM;
    data(13) = betaK;
    data(14) = betaK0;
    data(15) = betaKc;
    
    // Send the data vector
    res += theChannel.sendVector(this->getDbTag(), cTag, data);
    if (res < 0) {
      opserr << "ElasticBeamWarping3d::sendSelf -- could not send data Vector\n";
      return res;
    }

    // Ask the CoordTransf to send itself
    res += theCoordTransf->sendSelf(cTag, theChannel);
    if (res < 0) {
      opserr << "ElasticBeamWarping3d::sendSelf -- could not send CoordTransf\n";
      return res;
    }

    return res;
}

int
ElasticBeamWarping3d::recvSelf(int cTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  int res = 0;
  static Vector data(16);

  res += theChannel.recvVector(this->getDbTag(), cTag, data);
  if (res < 0) {
    opserr << "ElasticBeamWarping3d::recvSelf -- could not receive data Vector\n";
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
  
  alphaM = data(12);
  betaK = data(13);
  betaK0 = data(14);
  betaKc = data(15);
  
  // Check if the CoordTransf is null; if so, get a new one
  int crdTag = (int)data(10);
  if (theCoordTransf == 0) {
    theCoordTransf = theBroker.getNewCrdTransf(crdTag);
    if (theCoordTransf == 0) {
      opserr << "ElasticBeamWarping3d::recvSelf -- could not get a CrdTransf3d\n";
      exit(-1);
    }
  }
  
  // Check that the CoordTransf is of the right type; if not, delete
  // the current one and get a new one of the right type
  if (theCoordTransf->getClassTag() != crdTag) {
    delete theCoordTransf;
    theCoordTransf = theBroker.getNewCrdTransf(crdTag);
    if (theCoordTransf == 0) {
      opserr << "ElasticBeamWarping3d::recvSelf -- could not get a CrdTransf3d\n";
      exit(-1);
    }
  }
  
  // Now, receive the CoordTransf
  theCoordTransf->setDbTag((int)data(11));
  res += theCoordTransf->recvSelf(cTag, theChannel, theBroker);
  if (res < 0) {
    opserr << "ElasticBeamWarping3d::recvSelf -- could not receive CoordTransf\n";
    return res;
  }
  
  // Revert the crdtrasf to its last committed state
  theCoordTransf->revertToLastCommit();
  
  return res;
}

void
ElasticBeamWarping3d::Print(OPS_Stream &s, int flag)
{
   if (flag == -1) { 
    int eleTag = this->getTag();
    s << "EL_BEAM\t" << eleTag << "\t";
    s  << "\t" << connectedExternalNodes(0) << "\t" << connectedExternalNodes(1);
    s << "\t0\t0.0000000\n";
   }  else if (flag < -1) {
     int counter = (flag + 1) * -1;
     int eleTag = this->getTag();
     const Vector &force = this->getResistingForce();

    double P, MZ1, MZ2, VY, MY1, MY2, VZ, T;
    double L = theCoordTransf->getInitialLength();
    double oneOverL = 1.0/L;
    
    P   = q(0);
    MZ1 = q(1);
    MZ2 = q(2);
    VY  = (MZ1+MZ2)*oneOverL;
    MY1 = q(3);
    MY2 = q(4);
    VZ  = (MY1+MY2)*oneOverL;
    T   = q(5);

    s << "FORCE\t" << eleTag << "\t" << counter << "\t0";
    s << "\t" << -P+p0[0] << "\t"  <<  VY+p0[1] << "\t"  << -VZ+p0[3]  << endln;
    s << "FORCE\t" << eleTag << "\t" << counter << "\t1";
    s << "\t"  << P  << ' '  << -VY+p0[2] << ' ' << VZ+p0[4] << endln;
    s << "MOMENT\t" << eleTag << "\t" << counter << "\t0";
    s << "\t" << -T << "\t"  << MY1 << "\t" << MZ1 << endln;
    s << "MOMENT\t" << eleTag << "\t" << counter << "\t1";
    s << "\t" << T << ' ' << MY2 << ' '  <<  MZ2 << endln;
    
   }

   else if (flag == 2){
     this->getResistingForce(); // in case linear algo

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

     this->getResistingForce(); // in case linear algo

    s << "\nElasticBeamWarping3d: " << this->getTag() << endln;
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
ElasticBeamWarping3d::displaySelf(Renderer &theViewer, int displayMode, float fact)
{
    // first determine the end points of the quad based on
    // the display factor (a measure of the distorted image)
    const Vector &end1Crd = theNodes[0]->getCrds();
    const Vector &end2Crd = theNodes[1]->getCrds();	

    static Vector v1(3);
    static Vector v2(3);

    if (displayMode >= 0) {
      const Vector &end1Disp = theNodes[0]->getDisp();
      const Vector &end2Disp = theNodes[1]->getDisp();
      
      for (int i = 0; i < 3; i++) {
	v1(i) = end1Crd(i) + end1Disp(i)*fact;
	v2(i) = end2Crd(i) + end2Disp(i)*fact;    
      }
    } else {
      int mode = displayMode * -1;
      const Matrix &eigen1 = theNodes[0]->getEigenvectors();
      const Matrix &eigen2 = theNodes[1]->getEigenvectors();
      if (eigen1.noCols() >= mode) {
	for (int i = 0; i < 3; i++) {
	  v1(i) = end1Crd(i) + eigen1(i,mode-1)*fact;
	  v2(i) = end2Crd(i) + eigen2(i,mode-1)*fact;    
	}    
      } else {
	for (int i = 0; i < 3; i++) {
	  v1(i) = end1Crd(i);
	  v2(i) = end2Crd(i);
	}    
      }
    }
    
    return theViewer.drawLine (v1, v2, 1.0, 1.0);
}

Response*
ElasticBeamWarping3d::setResponse(const char **argv, int argc, OPS_Stream &output)
{

  Response *theResponse = 0;

  output.tag("ElementOutput");
  output.attr("eleType","ElasticBeamWarping3d");
  output.attr("eleTag",this->getTag());
  output.attr("node1",connectedExternalNodes[0]);
  output.attr("node2",connectedExternalNodes[1]);
  
  // global forces
  if (strcmp(argv[0],"force") == 0 || strcmp(argv[0],"forces") == 0 ||
      strcmp(argv[0],"globalForce") == 0 || strcmp(argv[0],"globalForces") == 0) {


    output.tag("ResponseType","Px_1");
    output.tag("ResponseType","Py_1");
    output.tag("ResponseType","Pz_1");
    output.tag("ResponseType","Mx_1");
    output.tag("ResponseType","My_1");
    output.tag("ResponseType","Mz_1");
    output.tag("ResponseType","Px_2");
    output.tag("ResponseType","Py_2");
    output.tag("ResponseType","Pz_2");
    output.tag("ResponseType","Mx_2");
    output.tag("ResponseType","My_2");
    output.tag("ResponseType","Mz_2");

    theResponse =  new ElementResponse(this, 2, P);

	// local forces
  } else if (strcmp(argv[0],"localForce") == 0 || strcmp(argv[0],"localForces") == 0) {

    output.tag("ResponseType","N_ 1");
    output.tag("ResponseType","Vy_1");
    output.tag("ResponseType","Vz_1");
    output.tag("ResponseType","T_1");
    output.tag("ResponseType","My_1");
    output.tag("ResponseType","Tz_1");
    output.tag("ResponseType","N_2");
    output.tag("ResponseType","Py_2");
    output.tag("ResponseType","Pz_2");
    output.tag("ResponseType","T_2");
    output.tag("ResponseType","My_2");
    output.tag("ResponseType","Mz_2");

    theResponse =  new ElementResponse(this, 3, P);

  // basic forces
  } else if (strcmp(argv[0],"basicForce") == 0 || strcmp(argv[0],"basicForces") == 0) {

    output.tag("ResponseType","N");
    output.tag("ResponseType","Mz_1");
    output.tag("ResponseType","Mz_2");
    output.tag("ResponseType","My_1");
    output.tag("ResponseType","My_2");
    output.tag("ResponseType","T");
    
    theResponse = new ElementResponse(this, 4, Vector(6));
  }  

  output.endTag(); // ElementOutput

  return theResponse;
}

int
ElasticBeamWarping3d::getResponse (int responseID, Information &eleInfo)
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
    
  case 4: // basic forces
    return eleInfo.setVector(q);

  default:
    return -1;
  }
}
