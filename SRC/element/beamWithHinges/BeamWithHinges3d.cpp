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
                                                                        
// $Revision: 1.19 $
// $Date: 2004-06-07 23:20:55 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/beamWithHinges/BeamWithHinges3d.cpp,v $

#include <BeamWithHinges3d.h>
#include <Element.h>
#include <Domain.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Matrix.h>
#include <Vector.h>						
#include <Node.h>
#include <MatrixUtil.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include <SectionForceDeformation.h>
#include <CrdTransf3d.h>

#include <Information.h>
#include <ElementResponse.h>
#include <ElementalLoad.h>
#include <Renderer.h>

Matrix BeamWithHinges3d::theMatrix(12,12);
Vector BeamWithHinges3d::theVector(12);
double BeamWithHinges3d::workArea[200];

BeamWithHinges3d::BeamWithHinges3d(void)
  :Element(0, ELE_TAG_BeamWithHinges3d),
   E(0.0), A(0.0), Iz(0.0), Iy(0.0), G(0.0), J(0.0),
   beta1(0.0), beta2(0.0), rho(0.0),
   theCoordTransf(0),
   connectedExternalNodes(2),
   kb(6,6), q(6), load(12),
   kbCommit(6,6), qCommit(6),
   initialFlag(0), maxIter(0), tolerance(0.0), sp(0)
{
  section[0] = 0;
  section[1] = 0;

  p0[0] = 0.0;
  p0[1] = 0.0;
  p0[2] = 0.0;
  p0[3] = 0.0;
  p0[4] = 0.0;

  v0[0] = 0.0;
  v0[1] = 0.0;
  v0[2] = 0.0;
  v0[3] = 0.0;
  v0[4] = 0.0;

  theNodes[0] = 0;
  theNodes[1] = 0;
}

BeamWithHinges3d::BeamWithHinges3d(int tag, int nodeI, int nodeJ,
				   double e, double a, double iz,
				   double iy, double g, double j,
				   SectionForceDeformation &sectionRefI, double lpi,
				   SectionForceDeformation &sectionRefJ, double lpj,
				   CrdTransf3d &coordTransf,
				   double r, int max, double tol)
  :Element(tag, ELE_TAG_BeamWithHinges3d),
   E(e), A(a), Iz(iz), Iy(iy), G(g), J(j),
   beta1(lpi), beta2(lpj), rho(r),
   theCoordTransf(0),
   connectedExternalNodes(2),
   kb(6,6), q(6), load(12),
   kbCommit(6,6), qCommit(6),
   initialFlag(0), maxIter(max), tolerance(tol), sp(0)
{
  if (E <= 0.0)  {
    opserr << "BeamWithHinges3d::BeamWithHinges3d -- input parameter E is <= 0.0\n";
    exit(-1);
  }

  if (A <= 0.0)  {
    opserr << "BeamWithHinges3d::BeamWithHinges3d -- input parameter A is <= 0.0\n";
    exit(-1);
  }
  
  if (Iz <= 0.0)  {
    opserr << "BeamWithHinges3d::BeamWithHinges3d -- input parameter Iz is <= 0.0\n";
    exit(-1);
  }

  if (Iy <= 0.0)  {
    opserr << "BeamWithHinges3d::BeamWithHinges3d -- input parameter Iy is <= 0.0\n";
    exit(-1);
  }
  
  if (G <= 0.0)  {
    opserr << "BeamWithHinges3d::BeamWithHinges3d -- input parameter G is <= 0.0\n";
    exit(-1);
  }

  if (J <= 0.0)  {
    opserr << "BeamWithHinges3d::BeamWithHinges3d -- input parameter J is <= 0.0\n";
    exit(-1);
  }
  
  // Get copies of sections
  section[0] = sectionRefI.getCopy();
  
  if (section[0] == 0) {
    opserr << "BeamWithHinges3d::BeamWithHinges3d -- failed to get copy of section I\n";
    exit(-1);
  }
  
  section[1] = sectionRefJ.getCopy();
  
  if (section[1] == 0) {
    opserr << "BeamWithHinges3d::BeamWithHinges3d -- failed to get copy of section J\n";
    exit(-1);
  }
  
  theCoordTransf = coordTransf.getCopy();
  
  if (theCoordTransf == 0) {
    opserr << "BeamWithHinges3d::BeamWithHinges3d -- failed to get copy of coordinate transformation\n";
    exit(-1);
  }
  
  connectedExternalNodes(0) = nodeI;
  connectedExternalNodes(1) = nodeJ;

  theNodes[0] = 0;
  theNodes[1] = 0;

  p0[0] = 0.0;
  p0[1] = 0.0;
  p0[2] = 0.0;
  p0[3] = 0.0;
  p0[4] = 0.0;

  v0[0] = 0.0;
  v0[1] = 0.0;
  v0[2] = 0.0;
  v0[3] = 0.0;
  v0[4] = 0.0;
}

BeamWithHinges3d::~BeamWithHinges3d(void)
{
  for (int i = 0; i < 2; i++)
    if (section[i] != 0)
      delete section[i];
  
  if (theCoordTransf)
    delete theCoordTransf;

  if (sp != 0)
    delete sp;
}

int 
BeamWithHinges3d::getNumExternalNodes(void) const
{
  return 2;
}

const ID &
BeamWithHinges3d::getExternalNodes(void)
{
  return connectedExternalNodes;
}

Node **
BeamWithHinges3d::getNodePtrs()
{
    return theNodes;
}


int 
BeamWithHinges3d::getNumDOF(void)
{
  return 12;
}

void
BeamWithHinges3d::setDomain(Domain *theDomain)
{
  //This function may be called after a beam is constructed, so
  //geometry may change.  Therefore calculate all beam geometry here.
  
  if(theDomain == 0) {
    theNodes[0] = 0;
    theNodes[1] = 0;
    return;
  }
  
  // set the node pointers and verify them
  this->setNodePtrs(theDomain);
  
  // call the DomainComponent version of the function
  this->DomainComponent::setDomain(theDomain);
  
  if (theCoordTransf->initialize(theNodes[0], theNodes[1]) != 0) {
    opserr << "BeamWithHinges3d::setDomain() -- failed to initialize coordinate transformation\n";
    exit(-1);
  }
  
  // get element length
  double L = theCoordTransf->getInitialLength();
  if (L == 0.0) {
    opserr << "BeamWithHinges3d::setDomain() -- element has zero length\n";
    exit(-1);
  }

  // Set up section interpolation and hinge lengths
  this->setHinges();

  if (initialFlag == 2)
    theCoordTransf->update();
  else
    this->update();
}

int
BeamWithHinges3d::commitState(void)
{
  int err = 0;

  // call element commitState to do any base class stuff
  if ((err = this->Element::commitState()) != 0) {
    opserr << "BeamWithHinges3d::commitState () - failed in base class\n";
  }    
  
  for (int i = 0; i < 2; i++) {
    if (section[i] != 0)
      err += section[i]->commitState();
  }
  
  err += theCoordTransf->commitState();
  
  kbCommit = kb;
  qCommit = q;

  eCommit[0] = e[0];
  eCommit[1] = e[1];
  
  //initialFlag = 0;
  
  return err;
}

int
BeamWithHinges3d::revertToLastCommit(void)
{
  int err = 0;
  
  // Revert the sections and then get their last commited
  // deformations, stress resultants, and flexibilities
  for (int i = 0; i < 2; i++) {
    if (section[i] != 0) {
      err += section[i]->revertToLastCommit();
      section[i]->setTrialSectionDeformation(eCommit[i]);

      e[i] = eCommit[i];
      sr[i] = section[i]->getStressResultant();
      fs[i] = section[i]->getSectionFlexibility();
    }
  }
  
  // Commit the coordinate transformation
  err += theCoordTransf->revertToLastCommit();
  
  kb = kbCommit;
  q = qCommit;
  
  initialFlag = 0;

  return err;
}

int
BeamWithHinges3d::revertToStart(void)
{
  int err = 0;
  
  for (int i = 0; i < 2; i++) {
    if (section[i] != 0) {
      err += section[i]->revertToStart();
      fs[i].Zero();
      e[i].Zero();
      sr[i].Zero();
      eCommit[i].Zero();
    }
  }
  
  err += theCoordTransf->revertToStart();
  
  kb.Zero();
  q.Zero();
  
  initialFlag = 0;
  this->update();

  return err;
}

const Matrix &
BeamWithHinges3d::getTangentStiff(void)
{
  return theCoordTransf->getGlobalStiffMatrix(kb, q);
}


const Matrix &
BeamWithHinges3d::getInitialStiff(void)
{
  double L = theCoordTransf->getInitialLength();
  double oneOverL = 1.0/L;

  // Section locations along element length ...
  double xi[2];

  // and their integration weights
  double lp[2];
  
  lp[0] = beta1*L;
  lp[1] = beta2*L;

  xi[0] = 0.5*lp[0];
  xi[1] = L-0.5*lp[1];
  
  // element properties
  static Matrix f(6,6);	// element flexibility
  static Matrix Iden(6,6);   // an identity matrix for matrix inverse
  Iden.Zero();
  int i;
  for (i = 0; i < 6; i++)
    Iden(i,i) = 1.0;

  // Length of elastic interior
  double Le = L-lp[0]-lp[1];
  double LoverEA   = Le/(E*A);
  double Lover3EIz = Le/(3*E*Iz);
  double Lover6EIz = 0.5*Lover3EIz;
  double Lover3EIy = Le/(3*E*Iy);
  double Lover6EIy = 0.5*Lover3EIy;
  double LoverGJ   = Le/(G*J);

  // Elastic flexibility of element interior
  static Matrix fe(4,4);
  fe(0,0) = fe(1,1) =  Lover3EIz;
  fe(0,1) = fe(1,0) = -Lover6EIz;
  fe(2,2) = fe(3,3) =  Lover3EIy;
  fe(2,3) = fe(3,2) = -Lover6EIy;
  
  // Equilibrium transformation matrix
  static Matrix B(4,4);
  B(0,0) = B(2,2) = 1.0 - beta1;
  B(1,1) = B(3,3) = 1.0 - beta2;
  B(0,1) = B(2,3) = -beta1;
  B(1,0) = B(3,2) = -beta2;
  
  // Transform the elastic flexibility of the element
  // interior to the basic system
  static Matrix fElastic(4,4);
  fElastic.addMatrixTripleProduct(0.0, B, fe, 1.0);
    
  // Set element flexibility to flexibility of elastic region
  f.Zero();
  f(0,0) = LoverEA;
  f(1,1) = fElastic(0,0);
  f(2,2) = fElastic(1,1);
  f(1,2) = fElastic(0,1);
  f(2,1) = fElastic(1,0);
  f(3,3) = fElastic(2,2);
  f(4,4) = fElastic(3,3);
  f(3,4) = fElastic(2,3);
  f(4,3) = fElastic(3,2);
  f(5,5) = LoverGJ;    
  
  for (i = 0; i < 2; i++) {
    
    if (section[i] == 0 || lp[i] <= 0.0)
      continue;
    
    // Get section information
    int order = section[i]->getOrder();
    const ID &code = section[i]->getType();
    
    Vector s(workArea, order);
    Vector ds(&workArea[order], order);
    Vector de(&workArea[2*order], order);
    
    Matrix fb(&workArea[3*order], order, 6);
    
    double xL = xi[i]*oneOverL;
    double xL1 = xL-1.0;
    
    // get section flexibility matrix
    const Matrix &fSec = section[i]->getInitialFlexibility();
      
    // integrate section flexibility matrix
    // f += (b^ fs * b) * lp[i];
    //f.addMatrixTripleProduct(1.0, b, fSec, lp[i]);
    int ii, jj;
    fb.Zero();
    double tmp;
    for (ii = 0; ii < order; ii++) {
      switch(code(ii)) {
      case SECTION_RESPONSE_P:
	for (jj = 0; jj < order; jj++)
	  fb(jj,0) += fSec(jj,ii)*lp[i];
	break;
      case SECTION_RESPONSE_MZ:
	for (jj = 0; jj < order; jj++) {
	  tmp = fSec(jj,ii)*lp[i];
	  fb(jj,1) += xL1*tmp;
	  fb(jj,2) += xL*tmp;
	}
	break;
      case SECTION_RESPONSE_VY:
	for (jj = 0; jj < order; jj++) {
	  //tmp = oneOverL*fSec(jj,ii)*lp[i]*L/lp[i];
	  tmp = fSec(jj,ii);
	  fb(jj,1) += tmp;
	  fb(jj,2) += tmp;
	}
	break;
      case SECTION_RESPONSE_MY:
	for (jj = 0; jj < order; jj++) {
	  tmp = fSec(jj,ii)*lp[i];
	  fb(jj,3) += xL1*tmp;
	  fb(jj,4) += xL*tmp;
	  }
	break;
      case SECTION_RESPONSE_VZ:
	for (jj = 0; jj < order; jj++) {
	  //tmp = oneOverL*fSec(jj,ii)*lp[i]*L/lp[i];
	  tmp = fSec(jj,ii);
	  fb(jj,3) += tmp;
	  fb(jj,4) += tmp;
	}
	break;
      case SECTION_RESPONSE_T:
	for (jj = 0; jj < order; jj++)
	  fb(jj,5) += fSec(jj,ii)*lp[i];
	break;
      default:
	break;
      }
    }
    for (ii = 0; ii < order; ii++) {
      switch (code(ii)) {
      case SECTION_RESPONSE_P:
	for (jj = 0; jj < 6; jj++)
	  f(0,jj) += fb(ii,jj);
	break;
      case SECTION_RESPONSE_MZ:
	for (jj = 0; jj < 6; jj++) {
	  tmp = fb(ii,jj);
	  f(1,jj) += xL1*tmp;
	  f(2,jj) += xL*tmp;
	}
	break;
      case SECTION_RESPONSE_VY:
	for (jj = 0; jj < 6; jj++) {
	  tmp = oneOverL*fb(ii,jj);
	  f(1,jj) += tmp;
	  f(2,jj) += tmp;
	}
	break;
      case SECTION_RESPONSE_MY:
	for (jj = 0; jj < 6; jj++) {
	  tmp = fb(ii,jj);
	  f(3,jj) += xL1*tmp;
	  f(4,jj) += xL*tmp;
	}
	break;
      case SECTION_RESPONSE_VZ:
	for (jj = 0; jj < 6; jj++) {
	  tmp = oneOverL*fb(ii,jj);
	  f(3,jj) += tmp;
	  f(4,jj) += tmp;
	}
	break;
      case SECTION_RESPONSE_T:
	for (jj = 0; jj < 6; jj++)
	  f(5,jj) += fb(ii,jj);
	break;
      default:
	break;
      }
    }
    
  }
  // calculate element stiffness matrix
  static Matrix kbInit(6,6);
  if (f.Solve(Iden,kbInit) < 0)
    opserr << "BeamWithHinges3d::update() -- could not invert flexibility\n";
  
  return theCoordTransf->getInitialGlobalStiffMatrix(kbInit);
}

const Matrix &
BeamWithHinges3d::getMass(void)
{
  theMatrix.Zero();

  if (rho != 0.0) {
    double L = theCoordTransf->getInitialLength();  
    theMatrix(0,0) = theMatrix(1,1) = theMatrix(2,2) = 
      theMatrix(6,6) = theMatrix(7,7) = theMatrix(8,8) = 0.5*L*rho;
  }
  
  return theMatrix;
}

void 
BeamWithHinges3d::zeroLoad(void)
{
  if (sp != 0)
    sp->Zero();

  p0[0] = 0.0;
  p0[1] = 0.0;
  p0[2] = 0.0;
  p0[3] = 0.0;
  p0[4] = 0.0;

  v0[0] = 0.0;
  v0[1] = 0.0;
  v0[2] = 0.0;
  v0[3] = 0.0;
  v0[4] = 0.0;

  load.Zero();
}

int
BeamWithHinges3d::addLoad(ElementalLoad *theLoad, double loadFactor)
{
  int type;
  const Vector &data = theLoad->getData(type, loadFactor);
  
  if (sp == 0) {
    sp = new Matrix(5,2);
    if (sp == 0)
      opserr << "BeamWithHinges3d::addLoad() -- out of memory\n";
  }

  double L = theCoordTransf->getInitialLength();
  double oneOverL = 1.0/L;

  double lp1 = beta1*L;
  double lp2 = beta2*L;
  double Le = L-lp1-lp2;

  // Section locations along element length ...
  double xi[2];
  xi[0] = 0.5*lp1;
  xi[1] = L-0.5*lp2;

  if (type == LOAD_TAG_Beam3dUniformLoad) {
    double wy = data(0)*loadFactor;  // Transverse
    double wz = data(1)*loadFactor;  // Transverse
    double wx = data(2)*loadFactor;  // Axial

    Matrix &s_p = *sp;

    // Accumulate applied section forces due to uniform load
    for (int i = 0; i < 2; i++) {
      double x = xi[i];
      // Axial
      s_p(0,i) += wx*(L-x);
      // Moment
      s_p(1,i) += wy*0.5*x*(x-L);
      // Shear
      s_p(2,i) += wy*(x-0.5*L);
      // Moment
      s_p(3,i) += wz*0.5*x*(x-L);
      // Shear
      s_p(4,i) += wz*(x-0.5*L);
    }

    // Accumulate reactions in basic system
    p0[0] -= wx*L;
    double V;
    V = 0.5*wy*L;
    p0[1] -= V;
    p0[2] -= V;
    V = 0.5*wz*L;
    p0[3] -= V;
    p0[4] -= V;

    // Quick return
    if (Le == 0.0)
      return 0;

    // Accumulate basic deformations due to uniform load on interior
    // Midpoint rule for axial
    v0[0] += wx*0.5*(L-lp1+lp2)/(E*A)*Le;

    // Two point Gauss for bending ... will not be exact when
    // hinge lengths are not equal, but this is not a big deal!!!
    double x1 = lp1 + 0.5*Le*(1.0-1.0/sqrt(3.0));
    double x2 = lp1 + 0.5*Le*(1.0+1.0/sqrt(3.0));

    double Mz1 = 0.5*wy*x1*(x1-L);
    double Mz2 = 0.5*wy*x2*(x2-L);

    double My1 = 0.5*wz*x1*(x1-L);
    double My2 = 0.5*wz*x2*(x2-L);

    double Le2EIz = Le/(2*E*Iz);
    double Le2EIy = Le/(2*E*Iy);

    double b1, b2;
    b1 = x1*oneOverL;
    b2 = x2*oneOverL;
    v0[2] += Le2EIz*(b1*Mz1+b2*Mz2);
    v0[4] += Le2EIy*(b1*My1+b2*My2);

    b1 -= 1.0;
    b2 -= 1.0;
    v0[1] += Le2EIz*(b1*Mz1+b2*Mz2);
    v0[3] += Le2EIy*(b1*My1+b2*My2);
  }

  else if (type == LOAD_TAG_Beam3dPointLoad) {
    double Py = data(0)*loadFactor;
    double Pz = data(1)*loadFactor;
    double N  = data(2)*loadFactor;
    double aOverL = data(3);
    double a = aOverL*L;

    double Vy2 = Py*aOverL;
    double Vy1 = Py-Vy2;

    double Vz2 = Pz*aOverL;
    double Vz1 = Pz-Vz2;

    Matrix &s_p = *sp;

    // Accumulate applied section forces due to point load
    for (int i = 0; i < 2; i++) {
      double x = xi[i];
      if (x <= a) {
	s_p(0,i) += N;
	s_p(1,i) -= x*Vy1;
	s_p(2,i) -= Vy1;
	s_p(3,i) -= x*Vz1;
	s_p(4,i) -= Vz1;
      }
      else {
	s_p(1,i) -= (L-x)*Vy2;
	s_p(2,i) += Vy2;
	s_p(3,i) -= (L-x)*Vz2;
	s_p(4,i) += Vz2;
      }
    }

    // Accumulate reactions in basic system
    p0[0] -= N;
    p0[1] -= Vy1;
    p0[2] -= Vy2;
    p0[3] -= Vz1;
    p0[4] -= Vz2;

    // Quick return
    if (Le == 0.0)
      return 0;

    // Accumulate basic deformations of interior due to point load
    double M1, M2, M3;
    double b1, b2, b3;

    // Point load is on left hinge
    if (a < lp1) {
      M1 = (lp1-L)*Vy2;
      M2 = -lp2*Vy2;

      double Le_6EI;
      Le_6EI = Le/(6*E*Iz);

      b1 = lp1*oneOverL;
      b2 = 1.0-lp2*oneOverL;
      v0[2] += Le_6EI*(M1*(2*b1+b2)+M2*(b1+2*b2));

      b1 -= 1.0;
      b2 -= 1.0;
      v0[1] += Le_6EI*(M1*(2*b1+b2)+M2*(b1+2*b2));

      M1 = (lp1-L)*Vz2;
      M2 = -lp2*Vz2;

      Le_6EI = Le/(6*E*Iy);

      b1 = lp1*oneOverL;
      b2 = 1.0-lp2*oneOverL;
      v0[4] += Le_6EI*(M1*(2*b1+b2)+M2*(b1+2*b2));

      b1 -= 1.0;
      b2 -= 1.0;
      v0[3] += Le_6EI*(M1*(2*b1+b2)+M2*(b1+2*b2));

      // Nothing to do for axial
      //v0[0] += 0.0;
    }
    // Point load is on right hinge
    else if (a > L-lp2) {
      M1 = -lp1*Vy1;
      M2 = (lp2-L)*Vy1;

      double Le_6EI;
      Le_6EI = Le/(6*E*Iz);

      b1 = lp1*oneOverL;
      b2 = 1.0-lp2*oneOverL;
      v0[2] += Le_6EI*(M1*(2*b1+b2)+M2*(b1+2*b2));

      b1 -= 1.0;
      b2 -= 1.0;
      v0[1] += Le_6EI*(M1*(2*b1+b2)+M2*(b1+2*b2));

      M1 = -lp1*Vz1;
      M2 = (lp2-L)*Vz1;

      Le_6EI = Le/(6*E*Iy);

      b1 = lp1*oneOverL;
      b2 = 1.0-lp2*oneOverL;
      v0[4] += Le_6EI*(M1*(2*b1+b2)+M2*(b1+2*b2));

      b1 -= 1.0;
      b2 -= 1.0;
      v0[3] += Le_6EI*(M1*(2*b1+b2)+M2*(b1+2*b2));

      v0[0] += N*Le/(E*A);      
    }
    // Point load is on elastic interior
    else {
      M1 = -lp1*Vy1;
      M2 = -lp2*Vy2;
      M3 = -a*Vy1;

      double L1_6EI;
      double L2_6EI;

      L1_6EI = (a-lp1)/(6*E*Iz);
      L2_6EI = (Le-a+lp1)/(6*E*Iz);

      b1 = lp1*oneOverL;
      b2 = 1.0-lp2*oneOverL;
      b3 = a*oneOverL;
      v0[2] += L1_6EI*(M1*(2*b1+b3)+M3*(b1+2*b3));
      v0[2] += L2_6EI*(M2*(2*b2+b3)+M3*(b2+2*b3));

      b1 -= 1.0;
      b2 -= 1.0;
      b3 -= 1.0;
      v0[1] += L1_6EI*(M1*(2*b1+b3)+M3*(b1+2*b3));
      v0[1] += L2_6EI*(M2*(2*b2+b3)+M3*(b2+2*b3));

      M1 = -lp1*Vz1;
      M2 = -lp2*Vz2;
      M3 = -a*Vz1;

      L1_6EI = (a-lp1)/(6*E*Iy);
      L2_6EI = (Le-a+lp1)/(6*E*Iy);

      b1 = lp1*oneOverL;
      b2 = 1.0-lp2*oneOverL;
      b3 = a*oneOverL;
      v0[4] += L1_6EI*(M1*(2*b1+b3)+M3*(b1+2*b3));
      v0[4] += L2_6EI*(M2*(2*b2+b3)+M3*(b2+2*b3));

      b1 -= 1.0;
      b2 -= 1.0;
      b3 -= 1.0;
      v0[3] += L1_6EI*(M1*(2*b1+b3)+M3*(b1+2*b3));
      v0[3] += L2_6EI*(M2*(2*b2+b3)+M3*(b2+2*b3));

      v0[0] += N*(a-lp1)/(E*A);
    }
  }

  else {
    opserr << "BeamWithHinges3d::addLoad() -- load type unknown for element with tag: " << this->getTag() << endln;
    return -1;
  }

  return 0;  
}

int
BeamWithHinges3d::addInertiaLoadToUnbalance(const Vector &accel)
{
  if (rho == 0.0)
    return 0;
  
  const Vector &Raccel1 = theNodes[0]->getRV(accel);
  const Vector &Raccel2 = theNodes[1]->getRV(accel);
  
  double L = theCoordTransf->getInitialLength();
  double mass = 0.5*L*rho;
  
  int i,j;
  for (i = 0, j = 6; i < 3; i++, j++) {
    load(i) += mass*Raccel1(i);
    load(j) += mass*Raccel2(i);	// Yes, this should be 'i'
  }
  
  return 0;
}

const Vector &
BeamWithHinges3d::getResistingForce(void)
{
  Vector p0Vec(p0, 5);

  return theCoordTransf->getGlobalResistingForce(q, p0Vec);
}

const Vector &
BeamWithHinges3d::getResistingForceIncInertia(void)
{
  theVector = this->getResistingForce();

  if (rho != 0.0) {
  
    double ag[12];
    
    const Vector &accel1 = theNodes[0]->getTrialAccel();
    const Vector &accel2 = theNodes[1]->getTrialAccel();
    
    ag[0] = accel1(0);
    ag[1] = accel1(1);
    ag[2] = accel1(2);
    ag[6] = accel2(0);
    ag[7] = accel2(1);
    ag[8] = accel2(2);
    
    theVector = this->getResistingForce();
    
    double L = theCoordTransf->getInitialLength();
    double mass = 0.5*L*rho;
    
    int i,j;
    for (i = 0, j = 6; i < 3; i++, j++) {
      theVector(i) += mass*ag[i];
      theVector(j) += mass*ag[j];
    }
  
    // add the damping forces if rayleigh damping
    if (alphaM != 0.0 || betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0)
      theVector += this->getRayleighDampingForces();
  
  } else {

    // add the damping forces if rayleigh damping
    if (betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0)
      theVector += this->getRayleighDampingForces();

  }
  return theVector;
}

int
BeamWithHinges3d::sendSelf(int commitTag, Channel &theChannel)
{
  // place the integer data into an ID
  int dbTag = this->getDbTag();
  int i, j , k;
  int loc = 0;
  
  static ID idData(11);  
  idData(0) = this->getTag();
  idData(1) = connectedExternalNodes(0);
  idData(2) = connectedExternalNodes(1);
  idData(3) = theCoordTransf->getClassTag();
  int theCoordTransfDbTag  = theCoordTransf->getDbTag();
  if (theCoordTransfDbTag  == 0) {
    theCoordTransfDbTag = theChannel.getDbTag();
    if (theCoordTransfDbTag  != 0) 
      theCoordTransf->setDbTag(theCoordTransfDbTag);
  }
  idData(4) = theCoordTransfDbTag;

  idData(5) = initialFlag;
  idData(6) = maxIter;

  loc = 7;
  for (i = 0; i<2; i++) {
    int sectClassTag = section[i]->getClassTag();
    int sectDbTag = section[i]->getDbTag();
    if (sectDbTag == 0) {
      sectDbTag = theChannel.getDbTag();
      section[i]->setDbTag(sectDbTag);
    }

    idData(loc) = sectClassTag;
    idData(loc+1) = sectDbTag;
    loc += 2;
  }  

  if (theChannel.sendID(dbTag, commitTag, idData) < 0) {
    opserr << "NLBeamColumn3d::sendSelf() - failed to send ID data\n";
    return -1;
  }    

  // send the coordinate transformation
  
  if (theCoordTransf->sendSelf(commitTag, theChannel) < 0) {
    opserr << "NLBeamColumn3d::sendSelf() - failed to send crdTranf\n";
    return -1;
  }      

  //
  // send the sections
  //
  
  for (j = 0; j<2; j++) {
    if (section[j]->sendSelf(commitTag, theChannel) < 0) {
      opserr << "NLBeamColumn3d::sendSelf() - section " << j << " failed to send itself\n";
      return -1;
    }
  }
  
  // into a vector place distrLoadCommit, rho, UeCommit, Secommit and kvcommit
  int secDefSize = 0;
  for (i = 0; i < 2; i++) {
     int size = section[i]->getOrder();
     secDefSize   += size;
  }

  Vector dData(10+6+36+secDefSize); 
  loc = 0;

  // place double variables into Vector
  dData(loc++) = E;
  dData(loc++) = A;
  dData(loc++) = Iz;
  dData(loc++) = Iy;
  dData(loc++) = G;
  dData(loc++) = J;
  dData(loc++) = beta1;
  dData(loc++) = beta2;
  dData(loc++) = rho;
  dData(loc++) = tolerance;
  
  // place kvcommit into vector
  for (i=0; i<6; i++) 
    dData(loc++) = qCommit(i);

  // place kvcommit into vector
  for (i=0; i<6; i++) 
     for (j=0; j<6; j++)
        dData(loc++) = kbCommit(i,j);

  // place vscommit into vector
  for (k=0; k<2; k++)
    for (i=0; i<section[k]->getOrder(); i++) 
	dData(loc++) = (eCommit[k])(i);

  if (theChannel.sendVector(dbTag, commitTag, dData) < 0) {
     opserr << "NLBeamColumn3d::sendSelf() - failed to send Vector data\n";
     return -1;
  }    

  return 0;
}

int 
BeamWithHinges3d::recvSelf(int commitTag, Channel &theChannel,
			       FEM_ObjectBroker &theBroker)
{
  // place the integer data into an ID
  int dbTag = this->getDbTag();
  int i, j , k;
  int loc = 0;
  
  static ID idData(11);  

  if (theChannel.recvID(dbTag, commitTag, idData) < 0)  {
    opserr << "NLBeamColumn3d::recvSelf() - failed to recv ID data\n";
    return -1;
  }    

  this->setTag(idData(0));
  connectedExternalNodes(0) = idData(1);
  connectedExternalNodes(1) = idData(2);

  maxIter = idData(5);
  initialFlag = idData(6);
  
  int crdTransfClassTag = idData(3);
  int crdTransfDbTag = idData(4);

  // create a new crdTransf object if one needed
  if (theCoordTransf == 0 || theCoordTransf->getClassTag() != crdTransfClassTag) {
      if (theCoordTransf != 0)
	  delete theCoordTransf;

      theCoordTransf = theBroker.getNewCrdTransf3d(crdTransfClassTag);

      if (theCoordTransf == 0) {
	opserr << "NLBeamColumn3d::recvSelf() - failed to obtain a CrdTrans object with classTag " <<
				crdTransfClassTag << endln;
	  return -2;	  
      }
  }

  theCoordTransf->setDbTag(crdTransfDbTag);

  // invoke recvSelf on the crdTransf obkject
  if (theCoordTransf->recvSelf(commitTag, theChannel, theBroker) < 0)  
  {
     opserr << "NLBeamColumn3d::sendSelf() - failed to recv crdTranf\n";
     return -3;
  }      

  //
  // receive the sections
  //

  loc = 7;
  if (section[0] == 0) {

    // if no sections yet created, we create new ones and then do a recvSelf
    for (i=0; i<2; i++) {
      int sectClassTag = idData(loc);
      int sectDbTag = idData(loc+1);
      loc += 2;
      section[i] = theBroker.getNewSection(sectClassTag);
      if (section[i] == 0) {
	opserr << "NLBeamColumn3d::recvSelf() - Broker could not create Section of class type" << sectClassTag << endln;
	exit(-1);
      }
      section[i]->setDbTag(sectDbTag);
      if (section[i]->recvSelf(commitTag, theChannel, theBroker) < 0) {
	opserr << "NLBeamColumn3d::recvSelf() - section " << i << " failed to recv itself\n";
	return -1;
      }     
    }

    this->setHinges();

  } else {

    // if sections exist, we ensure of correct type and then do a recvSelf
    for (i=0; i<2; i++) {

      int sectClassTag = idData(loc);
      int sectDbTag = idData(loc+1);
      loc += 2;

      // check of correct type
      if (section[i]->getClassTag() !=  sectClassTag) {
	// delete the old section[i] and create a new one
	delete section[i];
	section[i] = theBroker.getNewSection(sectClassTag);
	if (section[i] == 0) {
	  opserr << "NLBeamColumn3d::recvSelf() - Broker could not create Section of class type" << sectClassTag << endln;
	  exit(-1);
	}
      }

      // recvvSelf on it
      section[i]->setDbTag(sectDbTag);
      if (section[i]->recvSelf(commitTag, theChannel, theBroker) < 0) {
	opserr << "NLBeamColumn3d::recvSelf() - section " <<
				i << "failed to recv itself\n";
	return -1;
      }     
    }
  }


  // into a vector place distrLoadCommit, rho, UeCommit, Secommit and kvcommit
  int secDefSize = 0;
  for (i = 0; i < 2; i++) {
     int size = section[i]->getOrder();
     secDefSize   += size;
  }

  Vector dData(10+6+36+secDefSize); 
  loc = 0;

  if (theChannel.recvVector(dbTag, commitTag, dData) < 0) {
     opserr << "NLBeamColumn3d::sendSelf() - failed to send Vector data\n";
     return -1;
  }    

  // place double variables into Vector
   E = dData(loc++);
   A = dData(loc++);
   Iz = dData(loc++);
   Iy = dData(loc++);
   G = dData(loc++);
   J = dData(loc++);
   beta1 = dData(loc++);
   beta2 = dData(loc++);
   rho = dData(loc++);
   tolerance = dData(loc++);
  
  // place kvcommit into vector
  for (i=0; i<6; i++) 
    qCommit(i) = dData(loc++);

  // place kvcommit into vector
  for (i=0; i<6; i++) 
     for (j=0; j<6; j++)
        kbCommit(i,j) = dData(loc++);;

  // place vscommit into vector
  for (k=0; k<2; k++) 
    for (i=0; i<section[k]->getOrder(); i++) 
       (eCommit[k])(i) = dData(loc++);

  initialFlag = 2;

  return 0;
}

void 
BeamWithHinges3d::Print(OPS_Stream &s, int flag)
{
  s << "\nBeamWithHinges3d, tag: " << this->getTag() << endln;
  s << "\tConnected Nodes: " << connectedExternalNodes;
  s << "\tE:  " << E << endln;
  s << "\tA:  " << A << endln;
  s << "\tIz: " << Iz << endln;
  s << "\tIy: " << Iy << endln;
  s << "\tG:  " << G << endln;
  s << "\tJ:  " << J << endln;
  
  double P, Mz1, Mz2, Vy, My1, My2, Vz, T;
  double L = theCoordTransf->getInitialLength();
  double oneOverL = 1.0/L;

  P   = qCommit(0);
  Mz1 = qCommit(1);
  Mz2 = qCommit(2);
  Vy  = (Mz1+Mz2)*oneOverL;
  My1 = qCommit(3);
  My2 = qCommit(4);
  Vz  = (My1+My2)*oneOverL;
  T   = qCommit(5);

  s << "\tEnd 1 Forces (P Mz Vy My Vz T): "
    << -P+p0[0] << ' ' << Mz1 << ' ' <<  Vy+p0[1] << ' ' << My1 << ' ' <<  Vz+p0[3] << ' ' << -T << endln;
  s << "\tEnd 2 Forces (P Mz Vy My Vz T): "
    <<  P << ' ' << Mz2 << ' ' << -Vy+p0[2] << ' ' << My2 << ' ' << -Vz+p0[4] << ' ' <<  T << endln;
  
  if (section[0] != 0) {
    s << "Hinge 1, section tag: " << section[0]->getTag() << 
      ", length: " << beta1*L << endln;
    section[0]->Print(s,flag);
  }
  
  if (section[1] != 0) {
    s << "Hinge 2, section tag: " << section[2]->getTag() << 
      ", length: " << beta2*L << endln;
    section[1]->Print(s,flag);
  }
}

//////////////////////////////
//Private Function Definitions

void 
BeamWithHinges3d::setNodePtrs(Domain *theDomain)
{
  theNodes[0] = theDomain->getNode(connectedExternalNodes(0));
  theNodes[1] = theDomain->getNode(connectedExternalNodes(1));
  
  if(theNodes[0] == 0) {
    opserr << "BeamWithHinges3d::setNodePtrs() -- node 1 does not exist\n";
    exit(-1);
  }
  
  if(theNodes[1] == 0) {
    opserr << "BeamWithHinges3d::setNodePtrs() -- node 2 does not exist\n";
    exit(-1);
  }
  
  // check for correct # of DOF's
  int dofNd1 = theNodes[0]->getNumberDOF();
  int dofNd2 = theNodes[1]->getNumberDOF();
  if ((dofNd1 != 6) || (dofNd2 != 6))  {
    opserr << "BeamWithHinges3d::setNodePtrs() -- nodal dof is not three\n";
    exit(-1);	       
  }
}

int
BeamWithHinges3d::update(void)
{
  // if have completed a recvSelf() - do a revertToLastCommit
  // to get e, kb, etc. set correctly
  if (initialFlag == 2)
    this->revertToLastCommit();

  // Update the coordinate transformation
  theCoordTransf->update();
  
  // Convert to basic system from local coord's (eliminate rb-modes)
  static Vector v(6);				// basic system deformations
  v = theCoordTransf->getBasicTrialDisp();
  
  static Vector dv(6);
  dv = theCoordTransf->getBasicIncrDeltaDisp();
  
  double L = theCoordTransf->getInitialLength();
  double oneOverL = 1.0/L;

  // Section locations along element length ...
  double xi[2];
  // and their integration weights
  double lp[2];
  
  lp[0] = beta1*L;
  lp[1] = beta2*L;

  xi[0] = 0.5*lp[0];
  xi[1] = L-0.5*lp[1];
  
  // element properties
  static Matrix f(6,6);	// element flexibility
  static Vector vr(6);	// Residual element deformations
  
  static Matrix Iden(6,6);   // an identity matrix for matrix inverse
  Iden.Zero();
  for (int i = 0; i < 6; i++)
    Iden(i,i) = 1.0;

  // Length of elastic interior
  double Le = L-lp[0]-lp[1];
  double LoverEA   = Le/(E*A);
  double Lover3EIz = Le/(3*E*Iz);
  double Lover6EIz = 0.5*Lover3EIz;
  double Lover3EIy = Le/(3*E*Iy);
  double Lover6EIy = 0.5*Lover3EIy;
  double LoverGJ   = Le/(G*J);

  // Elastic flexibility of element interior
  static Matrix fe(4,4);
  fe(0,0) = fe(1,1) =  Lover3EIz;
  fe(0,1) = fe(1,0) = -Lover6EIz;
  fe(2,2) = fe(3,3) =  Lover3EIy;
  fe(2,3) = fe(3,2) = -Lover6EIy;
  
  // Equilibrium transformation matrix
  static Matrix B(4,4);
  B(0,0) = B(2,2) = 1.0 - beta1;
  B(1,1) = B(3,3) = 1.0 - beta2;
  B(0,1) = B(2,3) = -beta1;
  B(1,0) = B(3,2) = -beta2;
  
  // Transform the elastic flexibility of the element
  // interior to the basic system
  static Matrix fElastic(4,4);
  fElastic.addMatrixTripleProduct(0.0, B, fe, 1.0);

  // calculate nodal force increments and update nodal forces
  static Vector dq(6);
  //dq = kb * dv;   // using previous stiff matrix k,i
  dq.addMatrixVector(0.0, kb, dv, 1.0);
  
  for (int j = 0; j < maxIter; j++) {
    
    // q += dq;
    q.addVector(1.0, dq, 1.0);
    
    // Set element flexibility to flexibility of elastic region
    f.Zero();
    f(0,0) = LoverEA;
    f(1,1) = fElastic(0,0);
    f(2,2) = fElastic(1,1);
    f(1,2) = fElastic(0,1);
    f(2,1) = fElastic(1,0);
    f(3,3) = fElastic(2,2);
    f(4,4) = fElastic(3,3);
    f(3,4) = fElastic(2,3);
    f(4,3) = fElastic(3,2);
    f(5,5) = LoverGJ;    

    // vr = fElastic * q;
    vr(0) = LoverEA*q(0) + v0[0];
    vr(1) = fElastic(0,0)*q(1) + fElastic(0,1)*q(2) + v0[1];
    vr(2) = fElastic(1,0)*q(1) + fElastic(1,1)*q(2) + v0[2];
    vr(3) = fElastic(2,2)*q(3) + fElastic(2,3)*q(4) + v0[3];
    vr(4) = fElastic(3,2)*q(3) + fElastic(3,3)*q(4) + v0[4];
    vr(5) = LoverGJ*q(5);

    for (int i = 0; i < 2; i++) {
      
      if (section[i] == 0 || lp[i] <= 0.0)
	continue;
      
      // Get section information
      int order = section[i]->getOrder();
      const ID &code = section[i]->getType();
      
      Vector s(workArea, order);
      Vector ds(&workArea[order], order);
      Vector de(&workArea[2*order], order);
      
      Matrix fb(&workArea[3*order], order, 6);
      
      double xL = xi[i]*oneOverL;
      double xL1 = xL-1.0;
      
      int ii;
      // Section forces
      // s = b*q + bp*currDistrLoad;
      //this->getForceInterpMatrix(b, xi[i], code);
      //s.addMatrixVector(0.0, b, q, 1.0);
      for (ii = 0; ii < order; ii++) {
	switch(code(ii)) {
	case SECTION_RESPONSE_P:
	  s(ii) = q(0);
	  break;
	case SECTION_RESPONSE_MZ:
	  s(ii) = xL1*q(1) + xL*q(2);
	  break;
	case SECTION_RESPONSE_VY:
	  s(ii) = oneOverL*(q(1)+q(2));
	  break;
	case SECTION_RESPONSE_MY:
	  s(ii) = xL1*q(3) + xL*q(4);
	  break;
	case SECTION_RESPONSE_VZ:
	  s(ii) = oneOverL*(q(3)+q(4));
	  break;
	case SECTION_RESPONSE_T:
	  s(ii) = q(5);
	  break;
	default:
	  s(ii) = 0.0;
	  break;
	}
      }

      // Add the effects of element loads, if present
      if (sp != 0) {
	const Matrix &s_p = *sp;
	for (ii = 0; ii < order; ii++) {
	  switch(code(ii)) {
	  case SECTION_RESPONSE_P:
	    s(ii) += s_p(0,i);
	    break;
	  case SECTION_RESPONSE_MZ:
	    s(ii) += s_p(1,i);
	    break;
	  case SECTION_RESPONSE_VY:
	    s(ii) += s_p(2,i);
	    break;
	  case SECTION_RESPONSE_MY:
	    s(ii) += s_p(3,i);
	    break;
	  case SECTION_RESPONSE_VZ:
	    s(ii) += s_p(4,i);
	    break;
	  default:
	    break;
	  }
	}
      }
      
      // Increment in section forces
      // ds = s - sr
      ds = s;
      ds.addVector(1.0, sr[i], -1.0);
      
      // compute section deformation increments and update current deformations
      // e += fs * ds;
      de.addMatrixVector(0.0, fs[i], ds, 1.0);
      if (initialFlag != 0)
	e[i].addVector(1.0, de, 1.0);
      
      // set section deformations
      section[i]->setTrialSectionDeformation(e[i]);
      
      // get section resisting forces
      sr[i] = section[i]->getStressResultant();
      
      // get section flexibility matrix
      fs[i] = section[i]->getSectionFlexibility();
      
      // ds = s - sr;
      ds = s;
      ds.addVector(1.0, sr[i], -1.0);
      
      de.addMatrixVector(0.0, fs[i], ds, 1.0);
      
      // integrate section flexibility matrix
      // f += (b^ fs * b) * lp[i];
      //f.addMatrixTripleProduct(1.0, b, fSec, lp[i]);
      int jj;
      fb.Zero();
      double tmp;
      const Matrix &fSec = fs[i];
      for (ii = 0; ii < order; ii++) {
	switch(code(ii)) {
	case SECTION_RESPONSE_P:
	  for (jj = 0; jj < order; jj++)
	    fb(jj,0) += fSec(jj,ii)*lp[i];
	  break;
	case SECTION_RESPONSE_MZ:
	  for (jj = 0; jj < order; jj++) {
	    tmp = fSec(jj,ii)*lp[i];
	    fb(jj,1) += xL1*tmp;
	    fb(jj,2) += xL*tmp;
	  }
	  break;
	case SECTION_RESPONSE_VY:
	  for (jj = 0; jj < order; jj++) {
	    //tmp = oneOverL*fSec(jj,ii)*lp[i]*L/lp[i];
	    tmp = fSec(jj,ii);
	    fb(jj,1) += tmp;
	    fb(jj,2) += tmp;
	  }
	  break;
	case SECTION_RESPONSE_MY:
	  for (jj = 0; jj < order; jj++) {
	    tmp = fSec(jj,ii)*lp[i];
	    fb(jj,3) += xL1*tmp;
	    fb(jj,4) += xL*tmp;
	  }
	  break;
	case SECTION_RESPONSE_VZ:
	  for (jj = 0; jj < order; jj++) {
	    //tmp = oneOverL*fSec(jj,ii)*lp[i]*L/lp[i];
	    tmp = fSec(jj,ii);
	    fb(jj,3) += tmp;
	    fb(jj,4) += tmp;
	  }
	  break;
	case SECTION_RESPONSE_T:
	  for (jj = 0; jj < order; jj++)
	    fb(jj,5) += fSec(jj,ii)*lp[i];
	  break;
	default:
	  break;
	}
      }
      for (ii = 0; ii < order; ii++) {
	switch (code(ii)) {
	case SECTION_RESPONSE_P:
	  for (jj = 0; jj < 6; jj++)
	    f(0,jj) += fb(ii,jj);
	  break;
	case SECTION_RESPONSE_MZ:
	  for (jj = 0; jj < 6; jj++) {
	    tmp = fb(ii,jj);
	    f(1,jj) += xL1*tmp;
	    f(2,jj) += xL*tmp;
	  }
	  break;
	case SECTION_RESPONSE_VY:
	  for (jj = 0; jj < 6; jj++) {
	    tmp = oneOverL*fb(ii,jj);
	    f(1,jj) += tmp;
	    f(2,jj) += tmp;
	  }
	  break;
	case SECTION_RESPONSE_MY:
	  for (jj = 0; jj < 6; jj++) {
	    tmp = fb(ii,jj);
	    f(3,jj) += xL1*tmp;
	    f(4,jj) += xL*tmp;
	  }
	  break;
	case SECTION_RESPONSE_VZ:
	  for (jj = 0; jj < 6; jj++) {
	    tmp = oneOverL*fb(ii,jj);
	    f(3,jj) += tmp;
	    f(4,jj) += tmp;
	  }
	  break;
	case SECTION_RESPONSE_T:
	  for (jj = 0; jj < 6; jj++)
	    f(5,jj) += fb(ii,jj);
	  break;
	default:
	  break;
	}
      }
      
      // UNCOMMENT WHEN DISTRIBUTED LOADS ARE ADDED TO INTERFACE
      // vr.addMatrixVector(1.0, vElastic, currDistrLoad, 1.0);
      
      // vr += (b^ (e+de)) * lp[i];
      de.addVector(1.0, e[i], 1.0);
      //vr.addMatrixTransposeVector(1.0, b, de, lp[i]);
      for (ii = 0; ii < order; ii++) {
	switch(code(ii)) {
	case SECTION_RESPONSE_P:
	  vr(0) += de(ii)*lp[i];
	  break;
	case SECTION_RESPONSE_MZ:
	  tmp = de(ii)*lp[i];
	  vr(1) += xL1*tmp; vr(2) += xL*tmp;
	  break;
	case SECTION_RESPONSE_VY:
	  //tmp = oneOverL*de(ii)*lp[i]*L/lp[i];
	  tmp = de(ii);
	  vr(1) += tmp; vr(2) += tmp;
	  break;
	case SECTION_RESPONSE_MY:
	  tmp = de(ii)*lp[i];
	  vr(3) += xL1*tmp; vr(4) += xL*tmp;
	  break;
	case SECTION_RESPONSE_VZ:
	  //tmp = oneOverL*de(ii)*lp[i]*L/lp[i];
	  tmp = de(ii);
	  vr(3) += tmp; vr(4) += tmp;
	  break;
	case SECTION_RESPONSE_T:
	  vr(5) += de(ii)*lp[i];
	  break;
	default:
	  break;
	}
      }
    }
    
    // calculate element stiffness matrix
    if (f.Solve(Iden,kb) < 0)
      opserr << "BeamWithHinges3d::update() -- could not invert flexibility\n";
			      

    // dv = v - vr;
    dv = v;
    dv.addVector(1.0, vr, -1.0);
    
    // determine resisting forces
    // dq = kb * dv;
    dq.addMatrixVector(0.0, kb, dv, 1.0);
    
    double dW = dv^ dq;
    
    if (fabs(dW) < tolerance)
      break;
  }
  
  // q += dq;
  q.addVector(1.0, dq, 1.0);
  
  initialFlag = 1;
  
  return 0;
}

void
BeamWithHinges3d::setHinges(void)
{
  for (int i = 0; i < 2; i++) {
    if (section[i] == 0)
      continue;
    
    // Get the number of section response quantities
    int order = section[i]->getOrder();
    
    fs[i] = Matrix(order,order);
    e[i]  = Vector(order);
    sr[i] = Vector(order);
    eCommit[i] = Vector(order);    
  }
}

void
BeamWithHinges3d::getForceInterpMatrix(Matrix &b, double x, const ID &code)
{			
  b.Zero();
  
  double L = theCoordTransf->getInitialLength();
  double xi = x/L;
  
  for (int i = 0; i < code.Size(); i++) {
    switch (code(i)) {
    case SECTION_RESPONSE_MZ:		// Moment, Mz, interpolation
      b(i,1) = xi - 1.0;
      b(i,2) = xi;
      break;
    case SECTION_RESPONSE_P:		// Axial, P, interpolation
      b(i,0) = 1.0;
      break;
    case SECTION_RESPONSE_VY:		// Shear, Vy, interpolation
      b(i,1) = b(i,2) = 1.0/L;
      break;
    case SECTION_RESPONSE_MY:		// Moment, My, interpolation
      b(i,3) = xi - 1.0;
      b(i,4) = xi;
      break;
    case SECTION_RESPONSE_VZ:		// Shear, Vz, interpolation
      b(i,3) = b(i,4) = 1.0/L;
      break;
    case SECTION_RESPONSE_T:		// Torsion, T, interpolation
      b(i,5) = 1.0;
      break;
    default:
      break;
    }
  }
}

void
BeamWithHinges3d::getDistrLoadInterpMatrix(Matrix &bp, double x, const ID & code)
{
  bp.Zero();

  double L = theCoordTransf->getInitialLength();
  double xi = x/L;
  
  for (int i = 0; i < code.Size(); i++) {
    switch (code(i)) {
    case SECTION_RESPONSE_MZ:		// Moment, Mz, interpolation
      bp(i,1) = 0.5*xi*(xi-1);
      break;
    case SECTION_RESPONSE_P:		// Axial, P, interpolation
      bp(i,0) = 1.0-xi;
      break;
    case SECTION_RESPONSE_VY:		// Shear, Vy, interpolation
      bp(i,1) = xi-0.5;
      break;
    case SECTION_RESPONSE_MY:		// Moment, My, interpolation
      bp(i,2) = 0.5*xi*(xi-1);
      break;
    case SECTION_RESPONSE_VZ:		// Shear, Vz, interpolation
      bp(i,2) = xi-0.5;
      break;
    case SECTION_RESPONSE_T:		// Torsion, T, interpolation
      bp(i,3) = 1.0-xi;
      break;
    default:
      break;
    }
  }
}

Response*
BeamWithHinges3d::setResponse(const char **argv, int argc, Information &info)
{
  // hinge rotations
  if (strcmp(argv[0],"plasticDeformation") == 0 ||
      strcmp(argv[0],"plasticRotation") == 0)
    return new ElementResponse(this, 1, Vector(3));
  
  // global forces
  else if (strcmp(argv[0],"force") == 0 || strcmp(argv[0],"forces") == 0 ||
	   strcmp(argv[0],"globalForce") == 0 || strcmp(argv[0],"globalForces") == 0)
    return new ElementResponse(this, 2, theVector);
  
  // stiffness
  else if (strcmp(argv[0],"stiffness") == 0)
    return new ElementResponse(this, 3, theMatrix);
  
  // local forces
  else if (strcmp(argv[0],"localForce") == 0 || strcmp(argv[0],"localForces") == 0)
    return new ElementResponse(this, 4, theVector);
  
  // section response
  else if (strcmp(argv[0],"section") == 0) {
    int sectionNum = atoi(argv[1]) - 1;
    
    if (sectionNum >= 0 && sectionNum < 2)
      if (section[sectionNum] != 0)
	return section[sectionNum]->setResponse(&argv[2], argc-2, info);
    return 0;
  }
  else
    return 0;
}

int
BeamWithHinges3d::getResponse(int responseID, Information &eleInfo)
{
  double V, N, T, M1, M2;
  double L = theCoordTransf->getInitialLength();
  static Vector force(12);
  static Vector def(6);
  
  switch (responseID) {
  case 1: {
    const Vector &v = theCoordTransf->getBasicTrialDisp();
    double LoverEA   = L/(E*A);
    double Lover3EIz = L/(3*E*Iz);
    double Lover6EIz = 0.5*Lover3EIz;
    double LoverGJ   = L/(G*J);
    double Lover3EIy = L/(3*E*Iy);
    double Lover6EIy = 0.5*Lover3EIy;

    double q1 = qCommit(1);
    double q2 = qCommit(2);
    double q3 = qCommit(3);
    double q4 = qCommit(4);

    def(0) = v(0) - LoverEA*qCommit(0);
    def(1) = v(1) - Lover3EIz*q1 + Lover6EIz*q2;
    def(2) = v(2) + Lover6EIz*q1 - Lover3EIz*q2;
    def(3) = v(3) - Lover3EIy*q3 + Lover6EIy*q4;
    def(4) = v(4) + Lover6EIy*q3 - Lover3EIy*q4;
    def(5) = v(5) - LoverGJ*qCommit(5);

    return eleInfo.setVector(def);
  }
  
  case 2: // global forces
    return eleInfo.setVector(this->getResistingForce());
    
  case 3: // stiffness
    return eleInfo.setMatrix(this->getTangentStiff());
    
  case 4: // local forces
    // Axial
    N = q(0);
    force(6) =  N;
    force(0) = -N+p0[0];
    
    // Torsion
    T = q(5);
    force(9) =  T;
    force(3) = -T;
    
    // Moments about z and shears along y
    M1 = q(1);
    M2 = q(2);
    force(5)  = M1;
    force(11) = M2;
    V = (M1+M2)/L;
    force(1) =  V+p0[1];
    force(7) = -V+p0[2];
    
    // Moments about y and shears along z
    M1 = q(3);
    M2 = q(4);
    force(4)  = M1;
    force(10) = M2;
    V = (M1+M2)/L;
    force(2) = -V+p0[3];
    force(8) =  V+p0[4];

    return eleInfo.setVector(force);
    
  default:
    return -1;
  }
}

int
BeamWithHinges3d::displaySelf(Renderer &theViewer, int displayMode, float fact)
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

int
BeamWithHinges3d::setParameter(const char **argv, int argc, Information &info)
{
  // E of the beam interior
  if (strcmp(argv[0],"E") == 0) {
    info.theType = DoubleType;
    return 1;
  }
  
  // A of the beam interior
  else if (strcmp(argv[0],"A") == 0) {
    info.theType = DoubleType;
    return 3;
  }
  
  // I of the beam interior
  else if (strcmp(argv[0],"Iz") == 0) {
    info.theType = DoubleType;
    return 4;
  }
  
  // Section parameter
  else if (strcmp(argv[0],"section") ==0) {
    if (argc <= 2)
      return -1;
    
    int sectionNum = atoi(argv[1]);
    
    int ok = -1;
    
    if (sectionNum == 1)
      ok = section[0]->setParameter (&argv[2], argc-2, info);
    if (sectionNum == 2)
      ok = section[1]->setParameter (&argv[2], argc-2, info);
    
    if (ok < 0)
      return -1;
    else if (ok < 100)
      return sectionNum*100 + ok;
    else 
      return -1;
  }
  
  // Unknown parameter
  else
    return -1;
}	

int
BeamWithHinges3d::updateParameter(int parameterID, Information &info)
{
  switch (parameterID) {
  case 1:
    this->E = info.theDouble;
    return 0;
  case 3:
    this->A = info.theDouble;
    return 0;
  case 4:
    this->Iz = info.theDouble;
    return 0;
  default:
    if (parameterID >= 100) { // section quantity
      int sectionNum = parameterID/100; 
      if (sectionNum == 1)
	return section[0]->updateParameter (parameterID-100, info);
      if (sectionNum == 2)
	return section[1]->updateParameter (parameterID-2*100, info);
      else
	return -1;
    }
    else // unknown
      return -1;
  }	
}	
