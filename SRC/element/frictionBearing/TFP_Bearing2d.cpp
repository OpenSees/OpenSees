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
** file 'COPYRIGHT'  in main dInrectory for information on usage and   **
** redistribution,  and for a DISCLAIMER OF ALL WARRANTIES.           **
**                                                                    **
** Developed by:                                                      **
**   Frank McKenna (fmckenna@ce.berkeley.edu)                         **
**   Gregory L. Fenves (fenves@ce.berkeley.edu)                       **
**   Filip C. Filippou (filippou@ce.berkeley.edu)                     **
**                                                                    **
** ****************************************************************** */
                                                                        
// Written: Tracy Becker
//
// What: "@(#) TFP_Bearing2d.C, revA"

// we specify what header files we need
#include "TFP_Bearing2d.h"
#include <elementAPI.h>
#include <G3Globals.h>
#include <math.h>
#include <float.h>

#include <Information.h>
#include <Domain.h>
#include <Node.h>
#include <Channel.h>
#include <Message.h>
#include <FEM_ObjectBroker.h>
#include <UniaxialMaterial.h>
#include <Renderer.h>
#include <ElementResponse.h>

#include <math.h>
#include <stdlib.h>
#include <string.h>


// typical constructor
TFP_Bearing2d::TFP_Bearing2d(int tag, 
			     int Nd1, int Nd2, 
			     double *R, 
			     double *DOUT,
			     double *DIN,
			     double *MU,
			     double *H,
			     double h0,
			     double a,
			     double k)
  :Element(tag, ELE_TAG_TFP_Bearing2d),     
   externalNodes(2),
   H0(h0), Ac(a), Ap(a),
   numDOF(0), theMatrix(0), theInitialMatrix(0), theVector(0),
   theMatrixC(0), theVectorC(0)
{	
  
  K = k;
  
  // fill in the ID containing external node info with node id's    
  if (externalNodes.Size() != 2) {
    opserr << "FATAL TFP_Bearing2d::TFP_Bearing2d() - out of memory, could not create an ID of size 2\n";
    exit(-1);
  }
  
  externalNodes(0) = Nd1;
  externalNodes(1) = Nd2;        
  
  theNodes[0] = 0; 
  theNodes[1] = 0;
  
  for (int i=0; i<4; i++) {
    r[i]   = R[i];
    dOut[i] = DOUT[i];
    dIn[i]  = DIN[i];
    mu[i]  = MU[i];
    h[i]   = H[i];
  }
  
  double dh = 0.0;
  
  for (int i=0; i<8; i++) {
    vpCommit[i] = 0.0;
    vpTrial[i] = 0.0;
    vCommit[i]  = 0.0;
    vTrial[i]  = 0.0;
    vs[i] = 0.0;
    FrCommit[i] = 0.0;
    FrTrial[i] = 0.0;
    d[i] = 0.0; // r[i] - sqrt(r[i]*r[i]-sqrt(vs[i]*vs[i]+vs[i+4]*vs[i+4]));
    dh += d[i];
  }
  
  for (int i=0; i<4; i++) {
    PCommit[i] = 0.0;
    PTrial[i] = 0.0;
    UCommit[i] = 0.0;
    UTrial[i] = 0.0;
    N[i] = a;
  }
  
  HTrial  = H0 + dh;
}

// constructor which should be invoked by an FE_ObjectBroker only
TFP_Bearing2d::TFP_Bearing2d()
 :Element(0, ELE_TAG_TFP_Bearing2d), 
  externalNodes(2),
  numDOF(0), theMatrix(0), theInitialMatrix(0), theVector(0),
  theMatrixC(0), theVectorC(0)
{
  theNodes[0] = 0; 
  theNodes[1] = 0;
}

//  destructor - provided to clean up any memory
TFP_Bearing2d::~TFP_Bearing2d()
{
  if (theMatrix != 0)
    delete theMatrix;
  if (theInitialMatrix != 0)
    delete theInitialMatrix;
  if (theVector != 0)
    delete theVector;
  if (theMatrixC != 0)
    delete theMatrixC;
  if (theVectorC != 0)
    delete theVectorC;
}

int
TFP_Bearing2d::getNumExternalNodes(void) const
{
    return 2;
}

const ID &
TFP_Bearing2d::getExternalNodes(void) 
{
  return externalNodes;
}

Node **
TFP_Bearing2d::getNodePtrs(void) 
{
  return theNodes;
}

int
TFP_Bearing2d::getNumDOF(void) {
  return numDOF;
}

// method: setDomain()
//    to set a link to the enclosing Domain, ensure nodes exist in Domain
//    and set pointers to these nodes, also determines the length and 
//    transformation Matrix.

static int setDomainFlag = 0;

void
TFP_Bearing2d::setDomain(Domain *theDomain)
{
  // check Domain is not null - invoked when object removed from a domain
  if (theDomain == 0) {
    exit(-1);    
    return;
  }
  
  // first ensure nodes exist in Domain and set the node pointers
  Node *end1Ptr, *end2Ptr;
  int Nd1 = externalNodes(0);
  int Nd2 = externalNodes(1);
  end1Ptr = theDomain->getNode(Nd1);
  end2Ptr = theDomain->getNode(Nd2);	
  if (end1Ptr == 0) {
    opserr << "WARNING TFP_Bearing2d::setDomain() - at truss " << this->getTag() << " node " <<
      Nd1 << "  does not exist in domain\n";
    exit(-1);    
    return;  // don't go any further - otherwise segemntation fault
  }
  if (end2Ptr == 0) {        
    opserr << "WARNING TFP_Bearing2d::setDomain() - at truss " << this->getTag() << " node " <<
      Nd2 << "  does not exist in domain\n";
    exit(-1);    
    return;  // don't go any further - otherwise segemntation fault
  }	
  theNodes[0] = end1Ptr;
  theNodes[1] = end2Ptr;
  // call the DomainComponent class method THIS IS VERY IMPORTANT
  this->DomainComponent::setDomain(theDomain);
  
  // ensure connected nodes have correct number of dof's
  int dofNd1 = end1Ptr->getNumberDOF();
  int dofNd2 = end2Ptr->getNumberDOF();	
  if ((dofNd1 != dofNd2) || ((dofNd2 != 2) && (dofNd2 != 3)) ) {
    opserr << "TFP_Bearing2d::setDomain(): 2 or 3 dof required at nodes\n";
    exit(-1);
    return;
  }	

  if (dofNd2 == 2) {
    theMatrix = new Matrix(4,4);
    theVector = new Vector(4);
    numDOF = 4;
  } else {
    theMatrix = new Matrix(6,6);
    theVector = new Vector(6);
    numDOF = 6;
  }

  setDomainFlag = 1;
  this->update();
  theInitialMatrix = new Matrix(*theMatrix);
  theMatrixC = new Matrix(*theMatrix);
  theVectorC = new Vector(*theVector);
  setDomainFlag = 0;
}   	 


int
TFP_Bearing2d::commitState()
{
  for (int i=0; i<8; i++) {
    vpCommit[i] = vpTrial[i];
    vCommit[i] = vTrial[i];
    FrCommit[i] = FrTrial[i];
  }
  
  for (int i=0; i<4; i++) {
    PCommit[i] = PTrial[i];
    UCommit[i] = UTrial[i];
  }
  HCommit = HTrial;

  Ac = Ap;

  Domain *theDomain = this->getDomain();

  //  theDomain->calculateNodalReactions(1);
  //  const Vector &nd2Reactions = theNodes[1]->getReaction();
  //  Ac = nd2Reactions(2);

  //  opserr << "Ac: " << Ac << endln;
  *theVectorC = *theVector;
  *theMatrixC = *theMatrix;

  return 0;
}

int
TFP_Bearing2d::revertToLastCommit()
{
  for (int i=0; i<8; i++) {
    vpTrial[i] = vpCommit[i];
    vTrial[i] = vCommit[i];
    FrTrial[i] = FrCommit[i];
  }
  for (int i=0; i<4; i++) {
    PTrial[i] = PCommit[i];
    UTrial[i] = UCommit[i];
  }
  HTrial = HCommit;

  Ac=Ap;

  *theVector = *theVectorC;
  *theMatrix = *theMatrixC;

  return 0;
}

int
TFP_Bearing2d::revertToStart()
{
  for (int i=0; i<8; i++) {
    vpTrial[i] = 0.0;
    vTrial[i] = 0.0; 
    FrTrial[i] = 0.0;
    vpCommit[i] = 0.0;
    vCommit[i] = 0.0; 
    FrCommit[i] = 0.0;
  }
  for (int i=0; i<4; i++) {
    PTrial[i] = 0.0;
    UTrial[i] = 0.0;
    PCommit[i] = 0.0;
    UCommit[i] = 0.0;
  }
  HTrial = H0;
  return 0;
}


static Matrix kthat(4,4);
static Matrix ksrest(8,8);
static Matrix kt(8,8);
static Matrix ks(8,8);
static Matrix Af(8,8);

static Vector d(4);

static Matrix Kee(4,4);
static Matrix kii(4,4);
static Matrix kei(4,4);
static Matrix keeI(4,4);
static Matrix kieT(4,4);


int
TFP_Bearing2d::kt3Drma(double *v, double *vp, double *Fr, double A, double *P, double *vpi) {

  Vector vF (v, 8); 
  Vector vpF (vp, 8); 
  Vector FrF (Fr, 8); 
  Vector PF (P,4);  

  /*
  opserr << "v: " << vF;
  opserr << "vp: " << vpF;
  opserr << "Fr: " << FrF;
  opserr << "A: " << A << endln;
  opserr << "P: " << PF;
  */

  static double Ri[8];
  //  static double R[8];
  //  static double N[4];

  static Matrix kcont(8,8);
  static Matrix krot(8,8);

  int cont = 0;

  kthat.Zero(); 
  kt.Zero(); 
  ks.Zero(); 
  Af.Zero(); 
  kcont.Zero(); 
  krot.Zero();
			
  
  for (int i=0; i<4; i++) {
    int z=4+i;
    Ri[i]=sqrt((r[i]-h[i])*(r[i]-h[i]) - v[z]*v[z]);
    Ri[z]=sqrt((r[i]-h[i])*(r[i]-h[i]) - v[i]*v[i]);
  }
  
  double dh =0;
  for (int i=0; i<4; i++) {
    dh += d[i];
  }  
  
  for (int i=0; i<4; i++) {
    int z=4+i;
    double vyield=0.01;
    double qYield=mu[i]*A;
    double k0=qYield/vyield;

    //    opserr << "qYield: " << qYield << " K0: " << k0 << endln;

    //get trial shear forces of hysteretic component
    double qTrialx = k0*(v[i] -vs[i]- vp[i]);
    double qTrialy = k0*(v[z] -vs[z]- vp[z]);

    // compute yield criterion of hysteretic component
    double qTrialNorm = sqrt(qTrialx*qTrialx+qTrialy*qTrialy);
    double Y = qTrialNorm - qYield;
 
    // elastic step -> no updates for pastic displacements required
    if (Y <= 0 ) {
      // set tangent stiffnesses
      ks(i,i) = k0 + A/Ri[i];
      ks(z,z) = k0 + A/Ri[z];

      vpi[i] = vp[i];
      vpi[z] = vp[z];

    // plastic step -> return mapping
    } else {    
      // compute consistency parameters
      double dGamma = Y/k0;
      // update plastic displacements
      vpi[i] = vp[i] + dGamma*qTrialx/qTrialNorm;
      vpi[z] = vp[z] + dGamma*qTrialy/qTrialNorm;

      //  set tangent stiffnesses
      double qTrialNorm3 = qTrialNorm*qTrialNorm*qTrialNorm;
      ks(i,i) =  qYield*k0*qTrialy*qTrialy/qTrialNorm3 + A/Ri[i];
      ks(i,z) = -qYield*k0*qTrialx*qTrialy/qTrialNorm3;
      ks(z,i) = -qYield*k0*qTrialx*qTrialy/qTrialNorm3;
      ks(z,z) =  qYield*k0*qTrialx*qTrialx/qTrialNorm3 + A/Ri[z];
    }

    // restrainer contact stiffness
    double vt=sqrt(v[i]*v[i]+v[z]*v[z]); //local displacment of surface
    double rt=(dOut[i]-dIn[i])/2.0*(r[i]-h[i])/r[i];  //restrainer distance
    double del=0.1;

    if (vt>rt) {
      cont=1;
      double krim=k0*2;
      // set restrainer stiffnesses
      double vi2 = v[i]*v[i];
      double vz2 = v[z]*v[z];
      kcont(i,i) =  krim*v[i]*v[i]/(vi2+vz2);
      kcont(i,z) =  krim*v[z]*v[i]/(vi2+vz2);
      kcont(z,i) =  krim*v[z]*v[i]/(vi2+vz2);
      kcont(z,z) =  krim*v[z]*v[z]/(vi2+vz2);
       
      //force rotation matrix
      double F=sqrt(Fr[i]*Fr[i]+Fr[z]*Fr[z]);

      krot(i,i) =  F* ((v[i]+del)/sqrt((v[i]+del)*(v[i]+del)+vz2) - (v[i]-del)/sqrt((v[i]-del)*(v[i]-del)+vz2));
      krot(i,z) =  F* (v[i]/sqrt(vi2+(v[z]+del)*(v[z]+del)) - v[i]/sqrt(vi2+(v[z]-del)*(v[z]-del)));
      krot(z,i) =  F* (v[z]/sqrt((v[i]+del)*(v[i]+del)+vz2) - v[z]/sqrt((v[i]-del)*(v[i]-del)+vz2));
      krot(z,z) =  F* ((v[z]+del)/sqrt(vi2+(v[z]+del)*v[z]+del) - (v[z]-del)/sqrt(vi2+(v[z]-del)*v[z]-del));
    }
  }

  double del = 0.1;
  for (int i=0; i<8; i++)
    for (int j=0; j<8; j++)
      ksrest(i,j)=kcont(i,j)+krot(i,j)/(del * 2.0);

  //  opserr << "ksrest: " << ksrest;
  Af.Zero();
  Af(0,4) = Ri[0];
  Af(1,5) = Ri[1];

  Af(2,0) = -Ri[2]/(Ri[2]+Ri[3]);
  Af(2,1) =  Ri[2]/(Ri[2]+Ri[3]);
  Af(2,4) = -Ri[2]*(Ri[0]+Ri[3])/(Ri[2]+Ri[3]);
  Af(2,5) = Ri[2]*(-Ri[1]+Ri[3])/(Ri[2]+Ri[3]);

  Af(3,0) = -Ri[3]/(Ri[2]+Ri[3]);
  Af(3,1) = Ri[3]/(Ri[2]+Ri[3]);
  Af(3,4) = Ri[3]*(-Ri[0]+Ri[2])/(Ri[2]+Ri[3]);
  Af(3,5) = Ri[3]*(-Ri[2]-Ri[1])/(Ri[2]+Ri[3]);

  Af(4,6) = Ri[4];
  Af(5,7) = Ri[5];

  Af(6,2) = -Ri[6]/(Ri[6]+Ri[7]);
  Af(6,3) = Ri[6]/(Ri[6]+Ri[7]);
  Af(6,6) = -Ri[6]*(Ri[4]+Ri[7])/(Ri[6]+Ri[7]);
  Af(6,7) = Ri[6]*(-Ri[5]+Ri[7])/(Ri[6]+Ri[7]);

  Af(7,2) = -Ri[7]/(Ri[6]+Ri[7]);
  Af(7,3) = Ri[7]/(Ri[6]+Ri[7]);
  Af(7,6) = Ri[7]*(-Ri[4]+Ri[7])/(Ri[6]+Ri[7]);
  Af(7,7) = Ri[7]*(-Ri[6]-Ri[5])/(Ri[6]+Ri[7]);

  /*
  opserr << "Af: " << Af;
  opserr << "ks: " << ks;
  opserr << "ksrest: " << ksrest;
  */

  static Matrix KsPlusKsrest(8,8);

  KsPlusKsrest = ks;
  KsPlusKsrest += ksrest;

  // add kg
  KsPlusKsrest(0,2) += A/Ri[2];
  KsPlusKsrest(1,3) += A/Ri[3];
  KsPlusKsrest(4,6) += A/Ri[6];
  KsPlusKsrest(5,7) += A/Ri[7];

  //  opserr << "ksPlusksrest: " << KsPlusKsrest;
 
  kt.addMatrixTripleProduct(0.0, Af, KsPlusKsrest,1.0);

  //  opserr << "kt:" << kt;

  for (int i=0; i<4; i++) {
    for (int j=0; j<4; j++) {
      kthat(i,j) = kt(i,j);
      Kee(i,j) = kt(i+4, j+4);
      kei(i,j) = kt(i+4, j);
      kieT(j,i) = kt(i, j+4); // transpose for matrix triple product routine
    }
  }

  /*
  opserr << "ktHat:\n" << kthat;
  opserr << "kee:\n" << Kee;
  opserr << "kei:\n" << kei;
  */  
  Kee.Invert(keeI);

  //  opserr << "kee^-1:\n" << keeI;

  kthat.addMatrixTripleProduct(1.0, kieT, keeI, kei, -1.0);

  return cont;
}


int
TFP_Bearing2d::update()
{
  static Vector delU(4);
  static Vector delP(4);

  const Vector &v1 = theNodes[0]->getIncrDisp();
  const Vector &v2 = theNodes[1]->getIncrDisp();	

  double vpi[8];

  delU(0)=v1(0);
  delU(1)=v2(0);
  delU(2)=0.0;
  delU(3)=0.0;

  /*
  if (delU.Norm() < 1.0e-18 && setDomainFlag == 0) {
    opserr << "BEARING::update - returning small delU\n";
    return 0;
  }
  */

  opserr << theNodes[1]->getTrialDisp();
  opserr << "TFP::update:: delU: " << delU;

  int contC = kt3Drma(vCommit, vpCommit, FrCommit, Ac, PCommit, vpi);

  //  Vector vpiF (vpi,8);  
  //  opserr << "vpiF: " << vpiF;

  static Matrix stiffCommit(8,8);
  stiffCommit = ks;
  stiffCommit += ksrest;

  Vector PC(PCommit,4);
  Vector PT(PTrial, 4);

  delP = kthat*delU;

  PT = PC;
  PT += delP;

  for (int i=0; i<4; i++)
    UTrial[i] = UCommit[i] + delU(i);

  static Vector delU58(4);
  static Vector tmp1(4);
  
  tmp1.addMatrixVector(0.0, kei, delU, 1.0);
  delU58.addMatrixVector(0.0, keeI, tmp1, -1.0);

  static double dvData[8];
  static Vector dv(dvData,8);
  static Vector dFr(8);
  static Vector tmp2(8);

  for (int i=0; i<4; i++) {
    tmp2(i)=delU(i);
    tmp2(i+4)=delU58(i);
  }

  dv = Af * tmp2;

  //  opserr << "dv: " << dv;
  //  Vector vC(vCommit, 8); opserr << "vCommit: " << vC;

  for (int i=0; i<8; i++) {
    vTrial[i] = vCommit[i] + dvData[i];
    FrTrial[i] = FrCommit[i] + dFr(i);
  }

  //  Vector vT(vTrial, 8); opserr << "vTrial: " << vT;

  HTrial = H0 + dh;
  double vpit[8];

  int contT = kt3Drma(vTrial, vpCommit, FrTrial, Ac, PTrial, vpit);

  //  opserr << "vTrial: " << vT;
  // Vector FT(FrTrial, 8); opserr << "FrTrial: " << FT;

  // opserr << "Ptrial 2:" << PT;
  // opserr << "kthat: " << kthat;

  //   Vector vpiO(vpi, 8); opserr << "VPI 0: " << vpiO;

  static Matrix stiffTrial(8,8);

  stiffTrial = ks;
  stiffTrial += ksrest;

  opserr << "stiffTrial:\n " << stiffTrial;
  
  int subDiv = 0;
  for (int j=0; j<8; j++) { // j=1:8
    double change = fabs(stiffCommit(j,j)/stiffTrial(j,j));
    //<stiffTrial(j,j) || stiffCommit(j,j)>1.25*stiffTrial(j,j));
    int f=(change < 0.75 || change > 1.25);
    if (f==1 || contT==1 || contC==1)
      subDiv=1;
  }

  opserr << "update::subDiv: " << subDiv << endln;

  //  opserr << "subDIV: " << subDiv << " contT: " << contT << " contC: " << contC << endln;

  if (subDiv==1) {

    double dumax = 0.0001; 
    //double dumax = 0.001; 
    double maxDelU = 0.0;

    for (int i=0; i<4; i++) {
      double delUi = fabs(delU(i));
      //      opserr << "delUi: " << delUi << " maxDelU: " << maxDelU << endln;

      if (delUi > maxDelU)
	maxDelU = delUi;
    }

    int n=ceil(maxDelU/dumax);
    
    opserr << "n: " << n << "maxDelU: " << maxDelU << " dumax: " << dumax << endln;

    static Vector delu(4);
    delu =  delU;
    delu /= 1.0*n;

    opserr << "n: " << n << "maxDelU: " << maxDelU << " delu: " << delu << endln;

    //    opserr << "delu: " << delu;

    static double padd[4];
    static double uadd[4];
    static double Ptemp[4];

    static double vadd[8];
    static double vpTemp[8];
    static double FrTemp[8];

    for (int i=0; i<4; i++) {
      padd[i] = 0.0; 
      uadd[i]=0.0; 
    }
    for (int i=0; i<8; i++) {
      vadd[i]=0.0;
      FrTemp[i]=FrCommit[i];
      vpTemp[i] = vpCommit[i];
    }
    
    for (int j=0; j<n; j++) {

      for (int i=0; i<4; i++) {
	vTrial[i] = vCommit[i] + vadd[i];
	Ptemp[i] = PCommit[i] + padd[i];
      }

      contT = kt3Drma(vTrial, vpTemp, FrTemp, Ac, Ptemp, vpi);    

      //      Vector vpiJ(vpi, 8); opserr << "vpiJ: " << vpiJ;
      
      static Vector delp(4);

      delp.addMatrixVector(0.0, kthat, delu, 1.0);

      //      opserr << "delp: " << delp;
      
      for (int i=0; i<4; i++) {
	padd[i] += delp(i);
	uadd[i] += delu[i];
      }

      // delu58=-kt(5:8,5:8)^-1*kt(5:8,1:4)*delu;

      tmp1.addMatrixVector(0.0, kei, delu, 1.0);
      delU58.addMatrixVector(0.0, keeI, tmp1, -1.0);

      //      opserr << "delU58: " << delU58;
      
      // dv=Af*[delu;delu58];
      static Vector tmp2(8);
      for (int i=0; i<4; i++) {
	tmp2(i)=delu[i];
	tmp2(i+4)=delU58(i);
      }
      
      dv.addMatrixVector(0.0, Af, tmp2, 1.0);

      ///      opserr << "tmp2: " << tmp2;
      // opserr << "dv: " << dv;

      dFr.addMatrixVector(0.0, ksrest, dv, 1.0);

      for (int i=0; i<8; i++) {
	vadd[i] += dv(i); 
	vpTemp[i]=vpi[i];
	FrTemp[i] = FrTemp[i] + dFr(i);
      }
    }
    
    for (int i=0; i<8; i++) {
      FrTrial[i] = FrTemp[i];
      vTrial[i] = vCommit[i] + vadd[i];
    }
  
    for (int i=0; i<4; i++) {
      PTrial[i] = PCommit[i]+padd[i];
      UTrial[i] = UCommit[i]+uadd[i];
      vTrial[i] = vCommit[i]+vadd[i];
    }
  }

  for (int i=0; i<8; i++) {
    vpTrial[i] = vpi[i];
  }
  HTrial=H0+dh;

  theMatrix->Zero();
  theVector->Zero();

  int numD = numDOF/2;

  (*theVector)(0) = PTrial[0];
  (*theVector)(numD) = PTrial[1];
  
  (*theMatrix)(0,0) = kthat(0,0);
  (*theMatrix)(0,numD) = kthat(0,1);
  (*theMatrix)(numD,numD) = kthat(1,1);
  (*theMatrix)(numD,0) = kthat(1,0);
  
  const Vector &d1 = theNodes[0]->getTrialDisp();
  const Vector &d2 = theNodes[1]->getTrialDisp();	
  
  double axialDefo = d1(1)-d2(1);

  if (axialDefo >= 0) {
    (*theMatrix)(1,1) = K;
    (*theMatrix)(1,1+numD) = -K;
    (*theMatrix)(1+numD,1) = -K;
    (*theMatrix)(1+numD,1+numD) = K;
    
    double force = axialDefo*K;
    (*theVector)(1) = force;
    (*theVector)(1+numD) = -force;
    Ap = force;
  } else {
    double Kmin = K*DBL_EPSILON;
    (*theMatrix)(1,1) = Kmin; // uisng Kmin to keep system stable
    (*theMatrix)(1,1+numD) = -Kmin;
    (*theMatrix)(1+numD,1) = -Kmin;
    (*theMatrix)(1+numD,1+numD) = Kmin;
    
    double force = 0.0;
    (*theVector)(1) = force;
    (*theVector)(1+numD) = -force;
    Ap = force;
  }

  return 0;
}


const Matrix &
TFP_Bearing2d::getTangentStiff(void)
{
  opserr << "TFP_Bearing2d::tangent " << *theMatrix;
  return *theMatrix;
}

const Matrix &
TFP_Bearing2d::getInitialStiff(void)
{
  opserr << "TFP_Bearing2d::initialtangent " << *theInitialMatrix;
  return *theInitialMatrix;
}

const Vector &
TFP_Bearing2d::getResistingForce()
{	
  opserr << "TFP_Bearing2d::resisting " << *theVector;
  return *theVector;
}

int
TFP_Bearing2d::sendSelf(int commitTag, Channel &theChannel)
{
  return 0;
}

int
TFP_Bearing2d::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  return 0;
}

void
TFP_Bearing2d::Print(OPS_Stream &s, int flag)
{
  s << "Element: " << this->getTag(); 
  s << " type: TFP_Bearing2d  iNode: " << externalNodes(0);
  s << " jNode: " << externalNodes(1) << endln;
}



Response *
TFP_Bearing2d::setResponse(const char **argv, int argc, OPS_Stream &output)
{
  Response *theResponse = 0;

  output.tag("ElementOutput");
  output.attr("eleType",this->getClassType());
  output.attr("eleTag",this->getTag());
  int numNodes = this->getNumExternalNodes();
  const ID &nodes = this->getExternalNodes();
  static char nodeData[32];

  for (int i=0; i<numNodes; i++) {
    sprintf(nodeData,"node%d",i+1);
    output.attr(nodeData,nodes(i));
  }

  if (strcmp(argv[0],"force") == 0 || strcmp(argv[0],"forces") == 0 ||
      strcmp(argv[0],"globalForce") == 0 || strcmp(argv[0],"globalForces") == 0) {
    const Vector &force = this->getResistingForce();
    int size = force.Size();
    for (int i=0; i<size; i++) {
      sprintf(nodeData,"P%d",i+1);
      output.tag("ResponseType",nodeData);
    }
    theResponse = new ElementResponse(this, 1, this->getResistingForce());
  }


  output.endTag();
  return theResponse;
}



int 
TFP_Bearing2d::getResponse(int responseID, Information &eleInfo)
{
 double strain;
 // Vector res(this->getResistingForce());
 // res(2) = Ac;

  switch (responseID) {
  case -1:
    return -1;
  case 1: // global forces
    return eleInfo.setVector(this->getResistingForce());
    //  return eleInfo.setVector(res);
    
  case 2:
    return eleInfo.setVector(this->getRayleighDampingForces());
  default:
    return 0;
  }
}

