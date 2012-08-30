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
// What: "@(#) TFP_Bearing.C, revA"

// we specify what header files we need
#include "TFP_Bearing.h"
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

// initialise the class wide variables

static int init = 1;

static int elastic = 1;

static int numMyBearing = 0;

void *
OPS_TFP_Bearing()
{
  // print out a message about who wrote this element & any copyright info wanted
  if (numMyBearing == 0) {
    opserr << "TFP_Bearing::4582 element - Written by Tracy Becker, UC Berkeley Copyright 2011\n";
    numMyBearing++;
  }

  Element *theEle = 0;

  int numRemainingArgs = OPS_GetNumRemainingInputArgs();
  if (numRemainingArgs == 0) { // parallel processing
    theEle = new TFP_Bearing();
    return theEle;
  }

  if (numRemainingArgs != 25 && numRemainingArgs != 24 && numRemainingArgs != 26) {
    opserr << "ERROR - TFP_Bearing incorrect # args provided, want: element TFP_Bearing tag? iNode? jNode? ";
    opserr << "$R1 $R2 $R3 $R4 $do1 $do2 $do3 $do4 $din1 $din2 $din3 $din4 $mu1 $mu2 $mu3 $mu4";
    opserr << " $h1 $h2 $h3 $h4 $H0 <$a> <$K>\n";
    return theEle;
  }

  // get the id and end nodes 
  int iData[3];
  double dData[23];
  int numData;

  numData = 3;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid element data\n";
    return 0;
  }

  int eleTag = iData[0];

  if (numRemainingArgs == 24) {
    numData = 21;
    dData[21] = 0.0; // initial Axial Load = 0.0
    dData[22] = 1.0e12;
  } else if (numRemainingArgs == 25) {
    numData = 22;
    dData[22] = 1.0e12;
  } else {
    numData = 23;
  }

  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "WARNING error reading element area for element" << eleTag << endln;
    return 0;
  }

  // now create the truss and add it to the Domain
  int ndm = OPS_GetNDM();
  if (ndm == 3) {
    theEle = new TFP_Bearing(eleTag, 
			     iData[1], 
			     iData[2], 
			     &dData[0],
			     &dData[4],
			     &dData[8],
			     &dData[12],
			     &dData[16],
			     dData[20],
			     dData[21],
			     dData[22]);
  } else {
    theEle = new TFP_Bearing2d(eleTag, 
			       iData[1], 
			       iData[2], 
			       &dData[0],
			       &dData[4],
			       &dData[8],
			       &dData[12],
			       &dData[16],
			       dData[20],
			       dData[21],
			       dData[22]);
  }
    
  if (theEle == 0) {
    opserr << "WARNING ran out of memory creating element with tag " << eleTag << endln;
    return 0;
  }

  return theEle;
}


// typical constructor
TFP_Bearing::TFP_Bearing(int tag, 
			 int Nd1, int Nd2, 
			 double *R, 
			 double *DOUT,
			 double *DIN,
			 double *MU,
			 double *H,
			 double h0,
			 double a,
			 double k)
  :Element(tag, ELE_TAG_TFP_Bearing),
  externalNodes(2),
  H0(h0), Ac(a), Ap(a),
   numDOF(0), theMatrix(0), theInitialMatrix(0),theVector(0),
   AfTrial(8,8),AfCommit(8,8),ksrestTrial(8,8),ksrestCommit(8,8),
   keeTrial(4,4), keeCommit(4,4), keiTrial(4,4), keiCommit(4,4),
   stiffTrial(8,8), stiffCommit(8,8), kthatTrial(4,4), kthatCommit(4,4)
{	

  K = k;

  // fill in the ID containing external node info with node id's    
  if (externalNodes.Size() != 2) {
    opserr << "FATAL TFP_Bearing::TFP_Bearing() - out of memory, could not create an ID of size 2\n";
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
TFP_Bearing::TFP_Bearing()
  :Element(0, ELE_TAG_TFP_Bearing),
   externalNodes(2),
   numDOF(0), theMatrix(0), theInitialMatrix(0),theVector(0),
   AfTrial(8,8),AfCommit(8,8),ksrestTrial(8,8),ksrestCommit(8,8),
   keeTrial(4,4), keeCommit(4,4), keiTrial(4,4), keiCommit(4,4),
   stiffTrial(8,8), stiffCommit(8,8), kthatTrial(4,4), kthatCommit(4,4)
{
  theNodes[0] = 0; 
  theNodes[1] = 0;
}

//  destructor - provided to clean up any memory
TFP_Bearing::~TFP_Bearing()
{
  if (theMatrix != 0)
    delete theMatrix;
  if (theInitialMatrix != 0)
    delete theInitialMatrix;
  if (theVector != 0)
    delete theVector;
}

int
TFP_Bearing::getNumExternalNodes(void) const
{
    return 2;
}

const ID &
TFP_Bearing::getExternalNodes(void) 
{
  return externalNodes;
}

Node **
TFP_Bearing::getNodePtrs(void) 
{
  return theNodes;
}

int
TFP_Bearing::getNumDOF(void) {
  return numDOF;
}

// method: setDomain()
//    to set a link to the enclosing Domain, ensure nodes exist in Domain
//    and set pointers to these nodes, also determines the length and 
//    transformation Matrix.
void
TFP_Bearing::setDomain(Domain *theDomain)
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
    opserr << "WARNING TFP_Bearing::setDomain() - at truss " << this->getTag() << " node " <<
      Nd1 << "  does not exist in domain\n";
    exit(-1);    
    return;  // don't go any further - otherwise segemntation fault
  }
  if (end2Ptr == 0) {        
    opserr << "WARNING TFP_Bearing::setDomain() - at truss " << this->getTag() << " node " <<
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
  if ((dofNd1 != dofNd2) || ((dofNd2 != 3) && (dofNd2 != 6)) ) {
    opserr << "TFP_Bearing::setDomain(): 3 or 6 dof required at nodes\n";
    exit(-1);
    return;
  }	

  if (dofNd2 == 3) {
    theMatrix = new Matrix(6,6);
    theVector = new Vector(6);
    numDOF = 6;
  } else {
    theMatrix = new Matrix(12,12);
    theVector = new Vector(12);
    numDOF = 12;
  }

  static Vector delU(4);
  static Vector delP(4);

  double vpi[8];

  delU(0)=0.0;
  delU(1)=0.0;
  delU(2)=0.0;
  delU(3)=0.0;

  int contC = kt3Drma(vCommit, vpCommit, FrCommit, Ac, PCommit, vpi);

  theMatrix->Zero();
  theVector->Zero();

  int numD = numDOF/2;

  for (int i=0; i<2; i++) {
    (*theVector)(i) = PTrial[i];
    (*theVector)(i+numD) = PTrial[i+2];
    
    for (int j=0; j<2; j++) {
      (*theMatrix)(i,j) = kthatTrial(i,j);
      (*theMatrix)(i+numD,j+numD) = kthatTrial(i+2,j+2);
      (*theMatrix)(i+numD,j) = kthatTrial(i+2,j);
      (*theMatrix)(i,j+numD) = kthatTrial(i,j+2);
    }
  }
  init = 1;
  this->update();
  Ap = Ac;
  theInitialMatrix = new Matrix(*theMatrix);
  this->commitState();
  init = 0;
}   	 


int
TFP_Bearing::commitState()
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
  //  opserr << "commitState: " << Ac << endln;

  Ac = Ap;

  Domain *theDomain = this->getDomain();

  //  theDomain->calculateNodalReactions(1);
  //  const Vector &nd2Reactions = theNodes[1]->getReaction();
  //  Ac = nd2Reactions(2);

  //  opserr << "Ac: " << Ac << endln;
  AfCommit=AfTrial;
  ksrestCommit=ksrestTrial;
  keeCommit = keeTrial;
  keiCommit = keiTrial;
  stiffCommit=stiffTrial;
  kthatCommit=kthatTrial;

  return 0;
}

int
TFP_Bearing::revertToLastCommit()
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

  Ap=Ac;

  AfTrial=AfCommit;
  ksrestTrial=ksrestCommit;
  keeTrial = keeCommit;
  keiTrial = keiCommit;
  stiffTrial=stiffCommit;
  kthatTrial=kthatCommit;

  return 0;
}

int
TFP_Bearing::revertToStart()
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


//static Matrix kthat(4,4);
static Matrix ksrest(8,8);

static Matrix kt(8,8);
static Matrix ks(8,8);
static Matrix Af(8,8);

static Vector d(4);

static Matrix kei(4,4);
static Matrix kie(4,4);
static Matrix kee(4,4);

static Vector f(8);
static Vector F(8);

int
TFP_Bearing::kt3Drma(double *v, double *vp, double *Fr, double A, double *P, double *vpi) {

  //  opserr << "TFP::kt3Drma A: " << A << endln;
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
  static double R[8];
  static double N[4];

  static Matrix kcont(8,8);
  static Matrix krot(8,8);

  int cont = 0;

  kthatTrial.Zero(); 
  kt.Zero(); 
  ks.Zero(); 
  Af.Zero(); 
  kcont.Zero(); 
  krot.Zero();
  f.Zero();
  F.Zero();
			
  for (int i=0; i<4; i++)
    N[i] = A;
  
  for (int i=0; i<4; i++) {
    int z=4+i;
    Ri[i]=sqrt((r[i]-h[i])*(r[i]-h[i]) - v[z]*v[z]);
    Ri[z]=sqrt((r[i]-h[i])*(r[i]-h[i]) - v[i]*v[i]);

    /*
    d[i] = r[i] - sqrt(r[i]*r[i]) - sqrt(v[i]*v[i]+v[z]*v[z]);
    N[i] = A + sqrt((P[0]-P[2])*(P[0]-P[2]) + (P[1]-P[3])*(P[1]-P[3])) * 
      sqrt(v[i]*v[i]+v[z]*v[z])/r[i];
    */

    R[i] = Ri[i];
    R[z] = Ri[z];
  }
  
  double dh =0;
  for (int i=0; i<4; i++) {
    dh += d[i];
  }

  /*  
  R[0] = (Ri[0]*Ri[2])/(Ri[2]+fabs(v[2])*Ri[0]);
  R[1]=(Ri[1]*Ri[3])/(Ri[3]+fabs(v[3])*Ri[1]);
  R[4]=(Ri[4]*Ri[6])/(Ri[6]+fabs(v[4])*Ri[6]);
  R[5]=(Ri[5]*Ri[7])/(Ri[7]+fabs(v[5])*Ri[7]);

  double PNorm = 0.0;
  for (int i=0; i<4; i++) {
    PNorm += P[i]*P[i];
  }
  PNorm = sqrt(PNorm);
  N[0]=A+PNorm*(sqrt(v[0]*v[0]+v[4]*v[4])/r[0]+sqrt(v[2]*v[2]+v[6]*v[6])/r[2]);
  N[1]=A+PNorm*(sqrt(v[1]*v[1]+v[5]*v[5])/r[1]+sqrt(v[3]*v[3]+v[7]*v[7])/r[3]);
  */

  N[0] = A;
  N[1] = A;
  N[2] = A;
  N[3] = A;

  for (int i=0; i<4; i++) {
    int z=4+i;
    double vyield=0.01;
    double qYield=mu[i]*N[i];
    double k0=qYield/vyield;
    
    //get trial shear forces of hysteretic component
    double qTrialx = k0*(v[i] -vs[i]- vp[i]);
    double qTrialy = k0*(v[z] -vs[z]- vp[z]);

    // compute yield criterion of hysteretic component
    double qTrialNorm = sqrt(qTrialx*qTrialx+qTrialy*qTrialy);
    double Y = qTrialNorm - qYield;

    // elastic step -> no updates for pastic displacements required
    //    opserr << "Y: " << Y << endln;

    if (Y <= 0 ) {
      //      opserr <<  "ELASTIC STEP:\n";
      // set tangent stiffnesses
      //      opserr << "k0: " << k0 << " A: " << A << " R[i]: " << R[i] << endln;

      ks(i,i) = k0 + A/Ri[i]; // fmk last A/R[i];
      ks(z,z) = k0 + A/Ri[z]; // A/R[z];
      vpi[i] = vp[i];
      vpi[z] = vp[z];

      f(i)=qTrialx + A/Ri[i]*v[i];
      f(z)=qTrialy + A/Ri[z]*v[z];

    // plastic step -> return mapping
    } else {    
      /*
      opserr <<  "PLASTIC STEP:\n";
      opserr << "Y: " << Y << " qTriaxx: " << qTrialx << " qTrialy: " << qTrialy << " " << qYield << endln;

      opserr << k0 << " " << v[i] << " " << vs[i] << " " << vp[i] << endln;
      opserr << k0 << " " << v[z] << " " << vs[z] << " " << vp[z] << endln;
      */
      // compute consistency parameters
      double dGamma = Y/k0;
      // update plastic displacements
      vpi[i] = vp[i] + dGamma*qTrialx/qTrialNorm;
      vpi[z] = vp[z] + dGamma*qTrialy/qTrialNorm;
      //  set tangent stiffnesses
      double qTrialNorm3 = qTrialNorm*qTrialNorm*qTrialNorm;

      ks(i,i) =  qYield*k0*qTrialy*qTrialy/qTrialNorm3 + N[i]/Ri[i];
      ks(i,z) = -qYield*k0*qTrialx*qTrialy/qTrialNorm3;
      ks(z,i) = -qYield*k0*qTrialx*qTrialy/qTrialNorm3;
      ks(z,z) =  qYield*k0*qTrialx*qTrialx/qTrialNorm3 + N[i]/Ri[z];

      f(i)=qYield*qTrialx/qTrialNorm + A/Ri[i]*v[i];
      f(z)=qYield*qTrialy/qTrialNorm  + A/Ri[z]*v[z];
    }
    
    if (i==0) {
      f(i) = f(i)+A*v[2]/Ri[2];
      f(z) = f(z)+A*v[6]/Ri[6];
    } else if (i == 1) {
      f(i) = f(i)+A*v[3]/Ri[3];
      f(z) = f(z)+A*v[7]/Ri[7];
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

      f(i)=f(i)+krim*(vt-rt)*v[i]/vt;
      f(z)=f(i)+krim*(vt-rt)*v[z]/vt;
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
  Af(2,0) = Ri[2]/(Ri[2]+Ri[3]); 
  Af(2,2) = -Ri[2]/(Ri[2]+Ri[3]);
  Af(2,4) = -Ri[2]*(Ri[0]+Ri[3])/(Ri[2]+Ri[3]);
  Af(2,5) = Ri[2]*(-Ri[1]+Ri[3])/(Ri[2]+Ri[3]);
  Af(3,0) = Ri[3]/(Ri[2]+Ri[3]);
  Af(3,2) = -Ri[3]/(Ri[2]+Ri[3]);
  Af(3,4) = Ri[3]*(-Ri[0]+Ri[2])/(Ri[2]+Ri[3]);
  Af(3,5) = Ri[3]*(-Ri[2]-Ri[1])/(Ri[2]+Ri[3]);
  Af(4,6) = Ri[4];
  Af(5,7) = Ri[5];
  Af(6,1) = Ri[6]/(Ri[6]+Ri[7]);
  Af(6,3) = -Ri[6]/(Ri[6]+Ri[7]);
  Af(6,6) = -Ri[6]*(Ri[4]+Ri[7])/(Ri[6]+Ri[7]);
  Af(6,7) = Ri[6]*(-Ri[5]+Ri[7])/(Ri[6]+Ri[7]);
  Af(7,1) = Ri[7]/(Ri[6]+Ri[7]);
  Af(7,3) = -Ri[7]/(Ri[6]+Ri[7]);
  Af(7,6) = Ri[7]*(-Ri[4]+Ri[6])/(Ri[6]+Ri[7]);
  Af(7,7) = Ri[7]*(-Ri[6]-Ri[5])/(Ri[6]+Ri[7]);

  /*
    opserr << "Af: " << Af;
    opserr << "ks: " << ks;
    opserr << "ksrest: " << ksrest;    
  */
  static Matrix KsPlusKsrest(8,8);


  KsPlusKsrest = ks;
  KsPlusKsrest += ksrest;
  KsPlusKsrest(0,2) = KsPlusKsrest(0,2) + A/Ri[2];
  KsPlusKsrest(1,3) = KsPlusKsrest(1,3) + A/Ri[3];
  KsPlusKsrest(4,6) = KsPlusKsrest(4,6) + A/Ri[6];
  KsPlusKsrest(5,7) = KsPlusKsrest(5,7) + A/Ri[7];

  //  opserr << "ksPlusKsrest: " << KsPlusKsrest;

  kt.addMatrixTripleProduct(0.0, Af, KsPlusKsrest,1.0);

  //  opserr << "kt:" << kt;

  static Matrix Kee(4,4);

  for (int i=0; i<4; i++) {
    for (int j=0; j<4; j++) {
      kthatTrial(i,j) = kt(i,j);
      Kee(i,j) = kt(i+4, j+4);
      kei(i,j) = kt(i+4, j);
      kie(i,j) = kt(i, j+4);
    }
  }
  //  opserr << "K11: " << kthat;
  //  opserr << "K12: " << kie;
  //  opserr << "K21: " << kei;
  //  opserr << "K22: " << Kee;

  Kee.Invert(kee);
  Matrix tmp1 = kee * kei;
  Matrix tmp2 = kie* tmp1;

  kthatTrial.addMatrix(1.0, tmp2, -1.0);

  //  opserr << "kthat: " << kthat;

  /*
  Vector vpiE(vpi, 8); opserr << "VPI:kt3: " << vpiE;
  Vector p1(P, 4); opserr << "P:kt3: " << p1;
  opserr << "f:kt3: " << f;
  */

  return cont;
}


int
TFP_Bearing::update()
{
  //  opserr << "TFP_Bearing::update(): Ap: " << Ap << " Ac: " << Ac << endln;

  static Vector delU(4);
  static Vector delP(4);

  const Vector &v1 = theNodes[0]->getIncrDisp(); 
  const Vector &v2 = theNodes[1]->getIncrDisp();	

  const Vector &d1 = theNodes[0]->getTrialDisp();
  const Vector &d2 = theNodes[1]->getTrialDisp();	

  double axialDefo = d1(2)-d2(2);

  if (axialDefo > 0)
    Ap = axialDefo*K;

  double vpi[8];

  delU(0)=v1(0);
  delU(1)=v1(1);
  delU(2)=v2(0);
  delU(3)=v2(1);

  Vector PC(PCommit,4);
  Vector PT(PTrial, 4);
  delP = kthatCommit*delU;

  for (int i=0; i<8; i++) {
    vpTrial[i] = vpCommit[i];
    vTrial[i] = vCommit[i];
    FrTrial[i] = FrCommit[i];
  }

  for (int i=0; i<4; i++) {
    UTrial[i] = UCommit[i] + delU(i);
    PTrial[i] = PCommit[i] + delP(i);
  }

  static Vector delU58(4);
  static Vector tmp1(4);
  
  tmp1.addMatrixVector(0.0, kei, delU, 1.0);
  delU58.addMatrixVector(0.0, kee, tmp1, -1.0);
    
  static double dvData[8];
  static Vector dv(dvData,8);
  static Vector dFr(8);
  static Vector tmp2(8);

  for (int i=0; i<4; i++) {
    tmp2(i)=delU(i);
    tmp2(i+4)=delU58(i);
  }
  
  //  opserr << "delU delU58: " << tmp2;

  Af = AfCommit; 
  dv = Af * tmp2;

  dFr = ksrestCommit*dv; 

  for (int i=0; i<8; i++) {
    vTrial[i] = vCommit[i] + dvData[i];
    FrTrial[i] = FrCommit[i] + dFr(i);
  }

  HTrial = H0 + dh;
  //  double vpit[8];

  int contT = kt3Drma(vTrial, vpCommit, FrTrial, Ap, PTrial, vpi);

  static Vector tracyP(8);
  tracyP = Af^f;
  
  stiffTrial = ks;
  stiffTrial += ksrest;
  
  int subDiv = 0;
  for (int j=0; j<8; j++) { // j=1:8
    int f=(stiffTrial(j,j)*2<stiffCommit(j,j) || stiffTrial(j,j)>2*stiffCommit(j,j));
    if (f==1 || contT==1) // fmk || contC==1)
      subDiv=1;
  }

  if (init == 1)
    subDiv = 0;

  //  opserr << "subDIV: " << subDiv << " contT: " << contT << endln;
  if (subDiv==1) {

    double dumax = 0.0001; 
    double maxDelU = 0.0;

    for (int i=0; i<4; i++) {
      double delUi = fabs(delU(i));
      if (delUi > maxDelU)
	maxDelU = delUi;
    }

    int n=ceil(maxDelU/dumax);
    //    opserr << "n: " << n << "maxDelU: " << maxDelU << " dumax: " << dumax << endln;

    if (n < 50)
      n = 50;
    
    static Vector delu(4);
    static Vector delp(4);

    delu = delU;

    if (n != 0.0)
      delu /= 1.0*n;
    
    //    opserr << "n: " << n; 
    //    opserr << " delu: " << delu;
    
    static double padd[4];
    static double uadd[4];
    static double Ptemp[4];
    
    static double vadd[8];
    static double vpTemp[8];
    static double FrTemp[8];
    
    for (int i=0; i<4; i++) {
      padd[i] = 0.0; 
      uadd[i]=0.0; 
      PTrial[i] = PCommit[i];
      UTrial[i] = UCommit[i];
    }
    
    for (int i=0; i<8; i++) {
      vadd[i]=0.0;
      vpTrial[i] = vpCommit[i];
      FrTrial[i] = FrCommit[i];
      vTrial[i] = vCommit[i];
    }

    for (int j=0; j<n; j++) {
      //      opserr << "ITERATION: " << j << endln; 
      contT = kt3Drma(vTrial, vpTrial, FrTrial, Ap, PTrial, vpi);    

      // delu58=-kt(5:8,5:8)^-1*kt(5:8,1:4)*delu;
      //      opserr << "kei: " << kei;
      //      opserr << "kee: " << kee;

      tmp1.addMatrixVector(0.0, kei, delu, 1.0);
      delU58.addMatrixVector(0.0, kee, tmp1, -1.0);
      
      // dv=Af*[delu;delu58];
      static Vector tmp2(8);
      for (int i=0; i<4; i++) {
	tmp2(i)=delu[i];
	tmp2(i+4)=delU58(i);
      }
          
      dv.addMatrixVector(0.0, Af, tmp2, 1.0);

      //      opserr << "delu delu58: " << tmp2;
      //      opserr << "dv: " << dv;

      dFr.addMatrixVector(0.0, ksrest, dv, 1.0);

      static Vector delp(4);
      delp.addMatrixVector(0.0, kthatTrial, delu, 1.0);
      
      for (int i=0; i<4; i++) {
	padd[i] += delp(i);
	uadd[i] += delu[i];
	PTrial[i] += delp(i);
      }      

      for (int i=0; i<8; i++) {
	vTrial[i] += dv(i); 
	vpTrial[i]=vpi[i];
	FrTrial[i] += dFr(i);
      }
    }
  }

  /*      OLD FMK*/ 
  tracyP = Af^f;
  static Vector F14(4);
  static Vector F58(4);
  for (int i=0; i<4; i++) {
    F14(i)=tracyP(i);
    F58(i)=tracyP(i+4);
  }

  static Vector FinalP(4);
  FinalP = F14 - kie*kee*F58;

  //  opserr << "FINALP: " << FinalP << "PTRIAL: ";
  for (int i=0; i<4; i++) {
    //    opserr << PTrial[i] << " ";
    PTrial[i] = FinalP(i);
  }
  //  opserr << endln;

  for (int i=0; i<8; i++) {
    vpTrial[i] = vpi[i];
  }


  HTrial=H0+dh;

  theMatrix->Zero();
  theVector->Zero();

  int numD = numDOF/2;

  for (int i=0; i<2; i++) {
    (*theVector)(i) = PTrial[i];
    (*theVector)(i+numD) = PTrial[i+2];

    for (int j=0; j<2; j++) {
      (*theMatrix)(i,j) = kthatTrial(i,j);
      (*theMatrix)(i+numD,j+numD) = kthatTrial(i+2,j+2);
      (*theMatrix)(i+numD,j) = kthatTrial(i+2,j);
      (*theMatrix)(i,j+numD) = kthatTrial(i,j+2);
    }
  }

  if (axialDefo >= 0) {
    (*theMatrix)(2,2) = K;
    (*theMatrix)(2,2+numD) = -K;
    (*theMatrix)(2+numD,2) = -K;
    (*theMatrix)(2+numD,2+numD) = K;
    

    double force = axialDefo*K;
    (*theVector)(2) = force;
    (*theVector)(2+numD) = -force;
    Ap = force;
  } else {
    double Kmin = K*DBL_EPSILON;
    (*theMatrix)(2,2) = Kmin; // uisng Kmin to keep system stable
    (*theMatrix)(2,2+numD) = -Kmin;
    (*theMatrix)(2+numD,2) = -Kmin;
    (*theMatrix)(2+numD,2+numD) = Kmin;
    
    double force = 0.0;
    (*theVector)(2) = force;
    (*theVector)(2+numD) = -force;
    Ap = force;
  }

  AfTrial = Af;
  ksrestTrial = ksrest;
  keeTrial = kee;
  keiTrial=kei;

  return 0;
}


const Matrix &
TFP_Bearing::getTangentStiff(void)
{
  //  opserr << "TFP_Bearing::getTangentStiff(void)\n";
  //  opserr << *theMatrix;
  return *theMatrix;
}

const Matrix &
TFP_Bearing::getInitialStiff(void)
{
  return *theInitialMatrix;  
}

const Vector &
TFP_Bearing::getResistingForce()
{	
  //  opserr << "TFP_Bearing::getResistingForce(void)\n";
  //  opserr << *theVector;
  return *theVector;
}

int
TFP_Bearing::sendSelf(int commitTag, Channel &theChannel)
{
  return 0;
}

int
TFP_Bearing::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  return 0;
}

void
TFP_Bearing::Print(OPS_Stream &s, int flag)
{
  s << "Element: " << this->getTag(); 
  s << " type: TFP_Bearing  iNode: " << externalNodes(0);
  s << " jNode: " << externalNodes(1) << endln;
}



Response *
TFP_Bearing::setResponse(const char **argv, int argc, OPS_Stream &output)
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
TFP_Bearing::getResponse(int responseID, Information &eleInfo)
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

