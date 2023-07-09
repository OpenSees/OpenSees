



// Written: Quan Gu, Yichao Gao and Zhijian Qiu  
// Created: 2015/01/25 
// Sensitivity analysis of absorbing-transmitting element for boundaries of water
//------------------------------------------


#ifndef AV3D4QuadWithSensitivity_CPP
#define AV3D4QuadWithSensitivity_CPP

#include "AV3D4QuadWithSensitivity.h"
#include <math.h>

#define FixedOrder 2

// static variables
const int AV3D4QuadWithSensitivity::numDOF = 4;
const int AV3D4QuadWithSensitivity::nodes_in_elem = 4;
const int AV3D4QuadWithSensitivity::nodes_in_quad = 4;
const int AV3D4QuadWithSensitivity::r_integration_order = FixedOrder;
const int AV3D4QuadWithSensitivity::s_integration_order = FixedOrder;
const int AV3D4QuadWithSensitivity::dim = 3;
const int AV3D4QuadWithSensitivity::numGP = 4;

ID AV3D4QuadWithSensitivity::integFlags(0,numDOF);
ID AV3D4QuadWithSensitivity::actDOFs(0,numDOF);

Matrix AV3D4QuadWithSensitivity::K(numDOF,numDOF);
Matrix AV3D4QuadWithSensitivity::C(numDOF,numDOF);
Matrix AV3D4QuadWithSensitivity::M(numDOF,numDOF);
Vector AV3D4QuadWithSensitivity::P(numDOF);

Matrix ** AV3D4QuadWithSensitivity::H =0;
Matrix ** AV3D4QuadWithSensitivity::DH =0;
Matrix ** AV3D4QuadWithSensitivity::HH =0;

////Added by Qiu

Matrix AV3D4QuadWithSensitivity::CSensitivity(numDOF,numDOF);
#include <elementAPI.h>

void *
OPS_AV3D4QuadWithSensitivity(void){

  int eleID, numNodes, matTag;
  int nodes[8];
  int i;

  static int idData[6];

  //if the number of arguments is less than the minimum, throw an error
  int argc = OPS_GetNumRemainingInputArgs();
  if (argc != 6) {
    opserr << "element AV3D4Quad incorrect num args .. 6 expected\n";
    return 0;
  }

  if (OPS_GetIntInput(&argc, idData) != 0) {
    opserr << "element AV3D4Quad error reading integers\n";
    return 0;
  }  


  matTag = idData[5];
  NDMaterial *theMaterial = OPS_getNDMaterial(matTag);

  if (theMaterial == 0) {
    opserr << "command: element AC3D8Hex " << idData[0] << 
      " - no NDMaterial with tag " << matTag << " exists\n";
    return 0;      
  }

  Element *theEle = new AV3D4QuadWithSensitivity(idData[0],idData[1],idData[2],idData[3],idData[4], theMaterial);
  return theEle;
}



// 
// /////////////////////////////////////////////////////////////////////////////
//
// constructors and de-constructors
//

AV3D4QuadWithSensitivity::AV3D4QuadWithSensitivity(int element_number,
						   int node_numb_1, int node_numb_2, int node_numb_3, int node_numb_4)
:Element(element_number, ELE_TAG_AV3D4QuadWithSensitivity),
connectedExternalNodes(nodes_in_elem), Ki(0), hasConstrained(0)
{
  // Set connected external node IDs
  connectedExternalNodes(0) = node_numb_1;
  connectedExternalNodes(1) = node_numb_2;
  connectedExternalNodes(2) = node_numb_3;
  connectedExternalNodes(3) = node_numb_4;
  
  // zero node pointers
  for (int i = 0; i < nodes_in_elem; i++)
    theNodes[i] = 0;

  detJ = 0;

  // AddingSensitivity:BEGIN ////////////////////////////////
	parameterID = 0;
// AddingSensitivity:END /////////////////////////////////


}


AV3D4QuadWithSensitivity::AV3D4QuadWithSensitivity(int element_number,
						   int node_numb_1, int node_numb_2, int node_numb_3, int node_numb_4, 
						   NDMaterial * Globalmmodel)
:Element(element_number, ELE_TAG_AV3D4QuadWithSensitivity),
connectedExternalNodes(nodes_in_elem), Ki(0), hasConstrained(0)
{
  // Set connected external node IDs
  connectedExternalNodes(0) = node_numb_1;
  connectedExternalNodes(1) = node_numb_2;
  connectedExternalNodes(2) = node_numb_3;
  connectedExternalNodes(3) = node_numb_4;
  
  // zero node pointers
  for (int i = 0; i < nodes_in_elem; i++)
    theNodes[i] = 0;

  // check the material type
  const char *type = Globalmmodel->getType();
  if(strcmp(type, "AcousticMedium") != 0){
    opserr << "AV3D4QuadWithSensitivity::AV3D4QuadWithSensitivity - incompatible material model\n";
    exit(-1);
  }
  
  theMaterial = Globalmmodel;
  detJ = 0;
// AddingSensitivity:BEGIN ////////////////////////////////
	parameterID = 0;
// AddingSensitivity:END /////////////////////////////////

}

AV3D4QuadWithSensitivity::AV3D4QuadWithSensitivity ()
:Element(0, ELE_TAG_AV3D4QuadWithSensitivity ),
connectedExternalNodes(nodes_in_elem), Ki(0), hasConstrained(0)
{
  // zero node pointers
  for (int i = 0; i < nodes_in_elem; i++)
    theNodes[i] = 0;
  detJ = 0;

}

AV3D4QuadWithSensitivity::~AV3D4QuadWithSensitivity ()
{
  if (Ki != 0)
    delete Ki;

}

//
// /////////////////////////////////////////////////////////////////////////////
//

AV3D4QuadWithSensitivity & 
AV3D4QuadWithSensitivity::operator[](int subscript)
{
  return ( *(this+subscript) );
}

int 
AV3D4QuadWithSensitivity::getNumExternalNodes () const
{
  return nodes_in_elem;
}

const ID& 
AV3D4QuadWithSensitivity::getExternalNodes ()
{
  return connectedExternalNodes;
}

Node ** 
AV3D4QuadWithSensitivity::getNodePtrs(void)
{
  return theNodes;
}

int 
AV3D4QuadWithSensitivity::getNumDOF ()
{
  return numDOF;
}

//
// /////////////////////////////////////////////////////////////////////////////
//

void 
AV3D4QuadWithSensitivity::setDomain (Domain *theDomain)
{
  // Check Domain is not null - invoked when object removed from a domain
  if (theDomain == 0) {
    for (int i = 0; i < nodes_in_elem; i++) {
      theNodes[i] = 0;
    }
  }
  else {
    for (int i = 0; i < nodes_in_elem; i++) {
      theNodes[i] = theDomain->getNode(connectedExternalNodes(i));
      if (theNodes[i] == 0) {
        opserr << "FATAL ERROR AV3D4QuadWithSensitivity (tag: " << this->getTag();
        opserr << " ), node not found in domain\n";
        exit(-1);
      }
    }
    
    this->DomainComponent::setDomain(theDomain);
  }
}

int
AV3D4QuadWithSensitivity::update()
{
  // do nothing
  return 0;
}

int 
AV3D4QuadWithSensitivity::commitState ()
{
  int retVal = 0;

  // call element commitState to do any base class stuff
  if ((retVal = this->Element::commitState()) != 0) {
    opserr << "AV3D4QuadWithSensitivity::commitState () - failed in base class";
  }
  
  return retVal;
}

int 
AV3D4QuadWithSensitivity::revertToLastCommit ()
{ 
  // do nothing
  return 0;
}

int 
AV3D4QuadWithSensitivity::revertToStart () 
{
  // do nothing
  return 0;
} 

//
// /////////////////////////////////////////////////////////////////////////////
//

void 
AV3D4QuadWithSensitivity::zeroLoad(void)
{
  // Q.Zero();
}

int 
AV3D4QuadWithSensitivity::addLoad(ElementalLoad *theLoad, double loadFactor)
{
  // do nothing
  return 0;
}

int 
AV3D4QuadWithSensitivity::addInertiaLoadToUnbalance(const Vector &accel)
{
  // do nothing
  return 0;
}

//
// /////////////////////////////////////////////////////////////////////////////
//

const Vector &
AV3D4QuadWithSensitivity::getResistingForce(void) 
{
  P.Zero();
    Vector v(numDOF);
  int i;
  for(i = 0; i < nodes_in_elem; i++) {
    const Vector &vel = theNodes[i]->getTrialVel();
    v(i) = vel(0);
  }
  //opserr<<v<<endln;
  this->getDamp();
  
  P.addMatrixVector(1.0, C, v, 1.0);
    //  opserr<<"v is "<<v<<endln;
     // opserr<<"P IS "<<P<<endln;
  return P;

}

const Vector &
AV3D4QuadWithSensitivity::getResistingForceIncInertia(void) 
{
 /* P.Zero();
  
  Vector v(numDOF);
  int i;
  for(i = 0; i < nodes_in_elem; i++) {
    const Vector &vel = theNodes[i]->getTrialVel();
    v(i) = vel(0);
  }
  //opserr<<v<<endln;
  this->getDamp();
  
  P.addMatrixVector(1.0, C, v, 1.0);
 
  return P;
  */
	return this->getResistingForce();
}

//
// /////////////////////////////////////////////////////////////////////////////
//

int 
AV3D4QuadWithSensitivity::sendSelf (int commitTag, Channel &theChannel) 
{ 
  // Not implemtented yet
  return 0;
}

int 
AV3D4QuadWithSensitivity::recvSelf (int commitTag, Channel &theChannel, 
  FEM_ObjectBroker &theBroker) 
{   
  // Not implemtented yet
  return 0;
}

int AV3D4QuadWithSensitivity::displaySelf (Renderer &theViewer, int displayMode, float fact, const char **modes, int numMode)
{
  int error = 0;
  
  return error;

}     

void AV3D4QuadWithSensitivity::Print(OPS_Stream &s, int flag)
{
  if(flag == 1) {
    s << "AV3D4QuadWithSensitivity, element id:  " << this->getTag() << endln;
    s << "Connected external nodes:  " << connectedExternalNodes;
    s << this->getResistingForce();
  } 
  else {
    s << "AV3D4QuadWithSensitivity, element id:  " << this->getTag() << endln;
    s << "Connected external nodes:  " << connectedExternalNodes;
    for(int i = 0; i < numGP; i++) {
      theNodes[i]->Print(s);
    }
  }
}

//
// /////////////////////////////////////////////////////////////////////////////
//

Response * 
AV3D4QuadWithSensitivity::setResponse (const char **argv, int argc, OPS_Stream &output)
{
  Response *theResponse = 0;

  char outputData[32];

  output.tag("ElementOutput");
  output.attr("eleType","AV3D4QuadWithSensitivity");
  output.attr("eleTag",this->getTag());
  for(int i = 1; i <= nodes_in_elem; i++) {
    sprintf(outputData,"node%d",i);
    output.attr(outputData,theNodes[i-1]->getTag());
  }
  
  // // sample
  // if(_strcmpi(argv[0], "Stress") == 0) {
  //   theResponse =  new ElementResponse(this, 4, InfoS);
  // }
 
  // otherwise response quantity is unknown for the quad class
  
  output.endTag(); // ElementOutput
  return theResponse;
}

int 
AV3D4QuadWithSensitivity::getResponse (int responseID, Information &eleInfo)
{
  switch (responseID) {
    case 1:
      return eleInfo.setVector(this->getResistingForce());
      
    case 2:
      return eleInfo.setMatrix(this->getTangentStiff());
    
    default: 
      return -1;
  }
}

//
// /////////////////////////////////////////////////////////////////////////////
//

ID *
AV3D4QuadWithSensitivity::getActiveDofs(void)
{
  if( actDOFs.Size() == 0){
    for(int i = 0; i < nodes_in_elem; i++){
     actDOFs[i] = 8;
    }
  }

  return &actDOFs;
}

int 
AV3D4QuadWithSensitivity::getIntegrateFlag(void)
{
  return 1;
}

ID *
AV3D4QuadWithSensitivity::getIntegrateFlags(void)
{
  if (integFlags.Size() == 0){
    for (int i = 0; i < nodes_in_elem; i++){
      integFlags[i] = 1;
    }
  }

  return &integFlags;
}

//
// /////////////////////////////////////////////////////////////////////////////
//

Matrix 
AV3D4QuadWithSensitivity::getNodalCoords(void)
{
  Matrix N_Coord(nodes_in_quad,3);
  
  for(int i = 0; i < nodes_in_quad; i++) {
    const Vector &ndCrds = theNodes[i]->getCrds();
    N_Coord(i,0) = ndCrds(0);
    N_Coord(i,1) = ndCrds(1);
    N_Coord(i,2) = ndCrds(2);
  }

  return N_Coord;
}

//
// /////////////////////////////////////////////////////////////////////////////
//
// get stiffness matrix and mass matrix
//

const Matrix &
AV3D4QuadWithSensitivity::getTangentStiff(void)
{
  K.Zero();
  
  return K;
}

const Matrix &
AV3D4QuadWithSensitivity::getInitialStiff(void)
{
  if (Ki == 0) {
    Ki = new Matrix(this->getTangentStiff());
  }
  
  if (Ki == 0) {
    opserr << "FATAL fElement::getInitialStiff() -";
    opserr << "ran out of memory\n";
    exit(-1);
  }
    
  return *Ki;
}

const Matrix &
AV3D4QuadWithSensitivity::getDamp(void)
{
  // zero [C] first
  C.Zero();
  
  double Kf = 0.0;
  double rho = 0.0;
  double rw = 0.0;
  double sw = 0.0;
  short where = 0;
  
  // 
  const Matrix &D = theMaterial->getTangent();
  Kf = D(0,0);
  if(Kf == 0.0){
    opserr << "ERROR: The Kf is zero!\n";
    exit(-1);
  }
  
  rho = theMaterial->getRho();
  if(rho == 0.0){
    opserr << "ERROR: The rho is zero!\n";
    exit(-1);
  }
  
  // printf("<Kf,rho> is <%g,%g>;\n", Kf, rho);
  
  this->computeHH();
  this->computeDetJ();
  //opserr<<*HH[0]<<endln;
  
  double c1 = 1.0/sqrt(Kf*rho);
  // printf("c1 is %g;\n", c1);
  
  for(short GP_c_r = 1; GP_c_r <= r_integration_order; GP_c_r++) {
    rw = get_Gauss_p_w(r_integration_order, GP_c_r);
    for(short GP_c_s = 1; GP_c_s <= s_integration_order; GP_c_s++) {
      sw = get_Gauss_p_w(s_integration_order, GP_c_s);
      
      C.addMatrix(1.0, *HH[where], rw*sw*c1*detJ[where]);
      
      where++;
    }
  }
  //opserr<<"C IS "<<C<<endln;
  return C;
}

const Matrix &
AV3D4QuadWithSensitivity::getDampSensitivity(int gradNumber)
{

  // zero [C] first
  CSensitivity.Zero();
  double Kf = 0.0;
  double rho = 0.0; 
  double KfSensitivity = 0.0;
  double rhoSensitivity = 0.0;
  double rw = 0.0;
  double sw = 0.0;
  short where = 0;

  const Matrix &D = theMaterial->getTangent();
  Kf = D(0,0);

  const Matrix &DSensitivity = theMaterial->getDampTangentSensitivity(gradNumber);
  KfSensitivity = DSensitivity(0,0);
  
  rho = theMaterial->getRho();
  rhoSensitivity = theMaterial->getRhoSensitivity(gradNumber);
  
  // printf("<Kf,rho> is <%g,%g>;\n", Kf, rho);
  
  this->computeHH();
  this->computeDetJ();
  //opserr<<*HH[0]<<endln;
  
  double c1 = 1.0/sqrt(Kf*rho);
  double c1Sensitivity = -1.0/(2*sqrt(Kf*Kf*Kf*rho*rho*rho))*(KfSensitivity*rho+Kf*rhoSensitivity);

  // printf("c1 is %g;\n", c1);
  
  for(short GP_c_r = 1; GP_c_r <= r_integration_order; GP_c_r++) {
    rw = get_Gauss_p_w(r_integration_order, GP_c_r);
    for(short GP_c_s = 1; GP_c_s <= s_integration_order; GP_c_s++) {
      sw = get_Gauss_p_w(s_integration_order, GP_c_s);
      
      CSensitivity.addMatrix(1.0, *HH[where], rw*sw*c1Sensitivity*detJ[where]);
      
      where++;
    }
  }
  
  return CSensitivity;

}

const Matrix &
AV3D4QuadWithSensitivity::getMass(void)
{
  // zero [M] first
  M.Zero();
  
  return M;
}

const Matrix &
AV3D4QuadWithSensitivity::getConsMass(void)
{
  // zero [M] first
  M.Zero();
  
  return M;
}

//
// /////////////////////////////////////////////////////////////////////////////
//

// shape function
Matrix 
AV3D4QuadWithSensitivity::interp_fun(double r1, double r2)
{
  Matrix h(1,nodes_in_quad);

  // influence of the node number 4
  h(0,3)=(1.0-r1)*(1.0+r2)*0.25;
  // influence of the node number 3
  h(0,2)=(1.0+r1)*(1.0+r2)*0.25;
  // influence of the node number 2
  h(0,1)=(1.0+r1)*(1.0-r2)*0.25;
  // influence of the node number 1
  h(0,0)=(1.0-r1)*(1.0-r2)*0.25;
  
  return h;
}

Matrix 
AV3D4QuadWithSensitivity::diff_interp_fun(double r1, double r2)
{
  Matrix dh(2,nodes_in_quad);

  // influence of the node number 4
  dh(0,3)=-(1.0+r2)*0.25;
  dh(1,3)= (1.0-r1)*0.25;
  // influence of the node number 3
  dh(0,2)= (1.0+r2)*0.25;
  dh(1,2)= (1.0+r1)*0.25;
  // influence of the node number 2
  dh(0,1)= (1.0-r2)*0.25;
  dh(1,1)=-(1.0+r1)*0.25;
  // influence of the node number 1
  dh(0,0)=-(1.0-r2)*0.25;
  dh(1,0)=-(1.0-r1)*0.25;
  
  return dh;
}

int
AV3D4QuadWithSensitivity::computeH(void)
{
  if(H == 0 || DH == 0) {
    H = new Matrix*[numGP];
    DH = new Matrix*[numGP];
    if (H == 0 || DH == 0) {
      opserr << "AV3D4QuadWithSensitivity::computeH - out of memory!\n";
      return -3;
    }

    double r = 0.0;
    double s = 0.0;
    short where = 0;
    
    for(short GP_c_r = 1; GP_c_r <= r_integration_order; GP_c_r++) {
      r = get_Gauss_p_c(r_integration_order, GP_c_r);
      for(short GP_c_s = 1; GP_c_s <= s_integration_order; GP_c_s++) {
        s = get_Gauss_p_c(s_integration_order, GP_c_s);
        
        H[where] = new Matrix(1, nodes_in_quad);
        DH[where] = new Matrix(2, nodes_in_quad);
        if(H[where] == 0 || DH[where] == 0) {
          opserr << "AV3D4QuadWithSensitivity::computeH - out of memory!\n";
          return -3;
        }
        
        *H[where] = interp_fun(r, s);
        *DH[where] = diff_interp_fun(r, s);
        
        where++;
      } // for GP_c_s
    } // for GP_c_r
  }
  
  return 0;
}

int
AV3D4QuadWithSensitivity::computeHH(void)
{
  if (HH == 0) {
    HH = new Matrix*[numGP];
    if (HH == 0) {
      opserr << "AV3D4QuadWithSensitivity::computeHH - out of memory!\n";
      return -3;
    }
    
    // compute H first
    this->computeH();
    
    for(int i = 0; i < numGP; i++) {
      HH[i] = new Matrix(nodes_in_elem, nodes_in_elem);
      if (HH[i] == 0) {
        opserr << "AV3D4QuadWithSensitivity::computeHH - out of memory!\n";
        return -3;
      }
      
      HH[i]->addMatrixTransposeProduct(0.0, *H[i], *H[i], 1.0);
    }
  }
  
  return 0;
}

int 
AV3D4QuadWithSensitivity::computeDetJ(void)
{
  if(detJ == 0) {
    detJ = new double[numGP];
    if(detJ == 0) {
      opserr << "AV3D4QuadWithSensitivity::computeDetJ - out of memory!\n";
      return -3;
    }
    
    Matrix Jacobian(2,3);
    double x1,y1,z1,area;
    
    this->computeH();
    Matrix NC = getNodalCoords();
    
    for(int i = 0; i < numGP; i++) {
      Jacobian = (*DH[i])*NC;
      
      // pointer to the fluid
      x1 = Jacobian(0,1)*Jacobian(1,2) - Jacobian(0,2)*Jacobian(1,1);
      y1 = Jacobian(0,2)*Jacobian(1,0) - Jacobian(0,0)*Jacobian(1,2);
      z1 = Jacobian(0,0)*Jacobian(1,1) - Jacobian(0,1)*Jacobian(1,0);
      
      // length of the bar tangent
      area = sqrt(x1*x1 + y1*y1 + z1*z1);
      if(area == 0.0){
        opserr <<"The length of tangent should not be 0!\n";
        exit(-1);
      }
      
      detJ[i] = area;
      printf("detJ[%d] = %g;\n", i+1, detJ[i]);
    }
  }
  
  return 0;
}


/*
int
AV3D4QuadWithSensitivity::computeDiff(void)
{
  if (L == 0 || detJ == 0) {
    L = new Matrix*[numGP];
    detJ = new double[numGP];
    if (L == 0 || detJ == 0) {
      opserr << "AV3D4QuadWithSensitivity::computeDiff - out of memory!\n";
      return -3;
    }
    
    Matrix Jacobian(3,3);
    
    this->computeH();
    Matrix NC = getNodalCoords();
    
    for(int i = 0; i < numGP; i++) {
      L[i] = new Matrix(3, nodes_in_elem);
      if(L[i] == 0) {
        opserr << "AV3D4QuadWithSensitivity::computDiff() - out of memory!\n";
        return -3;
      }
      
      Matrix &dh = *DH[i];
      Jacobian = dh*NC;
      detJ[i] = Jacobian_det(Jacobian);
      
      // compute [Jac]^-1*[DH]
      Jacobian.Solve(dh, *L[i]);
    }
  }
  
  return 0;
}
*/

//
// /////////////////////////////////////////////////////////////////////////////
//




//
// /////////////////////////////////////////////////////////////////////////////
//

double 
AV3D4QuadWithSensitivity::get_Gauss_p_c(short order, short point_numb)
{
// Abscissae coefficient of the Gaussian quadrature formula
// starting from 1 not from 0
  static double Gauss_coordinates[7][7];
  
  Gauss_coordinates[1][1] = 0.0 ;
  Gauss_coordinates[2][1] = -0.577350269189626;
  Gauss_coordinates[2][2] = -Gauss_coordinates[2][1];
  Gauss_coordinates[3][1] = -0.774596669241483;
  Gauss_coordinates[3][2] = 0.0;
  Gauss_coordinates[3][3] = -Gauss_coordinates[3][1];
  Gauss_coordinates[4][1] = -0.861136311594053;
  Gauss_coordinates[4][2] = -0.339981043584856;
  Gauss_coordinates[4][3] = -Gauss_coordinates[4][2];
  Gauss_coordinates[4][4] = -Gauss_coordinates[4][1];
  Gauss_coordinates[5][1] = -0.906179845938664;
  Gauss_coordinates[5][2] = -0.538469310105683;
  Gauss_coordinates[5][3] = 0.0;
  Gauss_coordinates[5][4] = -Gauss_coordinates[5][2];
  Gauss_coordinates[5][5] = -Gauss_coordinates[5][1];
  Gauss_coordinates[6][1] = -0.932469514203152;
  Gauss_coordinates[6][2] = -0.661209386466265;
  Gauss_coordinates[6][3] = -0.238619186083197;
  Gauss_coordinates[6][4] = -Gauss_coordinates[6][3];
  Gauss_coordinates[6][5] = -Gauss_coordinates[6][2];
  Gauss_coordinates[6][6] = -Gauss_coordinates[6][1];
  
  return Gauss_coordinates[order][point_numb];
 }

double 
AV3D4QuadWithSensitivity::get_Gauss_p_w(short order, short point_numb)
{
// Weight coefficient of the Gaussian quadrature formula
// starting from 1 not from 0
  static double Gauss_weights[7][7]; // static data ??

  Gauss_weights[1][1] = 2.0;
  Gauss_weights[2][1] = 1.0;
  Gauss_weights[2][2] = 1.0;
  Gauss_weights[3][1] = 0.555555555555556;
  Gauss_weights[3][2] = 0.888888888888889;
  Gauss_weights[3][3] = Gauss_weights[3][1];
  Gauss_weights[4][1] = 0.347854845137454;
  Gauss_weights[4][2] = 0.652145154862546;
  Gauss_weights[4][3] = Gauss_weights[4][2];
  Gauss_weights[4][4] = Gauss_weights[4][1];
  Gauss_weights[5][1] = 0.236926885056189;
  Gauss_weights[5][2] = 0.478628670499366;
  Gauss_weights[5][3] = 0.568888888888889;
  Gauss_weights[5][4] = Gauss_weights[5][2];
  Gauss_weights[5][5] = Gauss_weights[5][1];
  Gauss_weights[6][1] = 0.171324492379170;
  Gauss_weights[6][2] = 0.360761573048139;
  Gauss_weights[6][3] = 0.467913934572691;
  Gauss_weights[6][4] = Gauss_weights[6][3];
  Gauss_weights[6][5] = Gauss_weights[6][2];
  Gauss_weights[6][6] = Gauss_weights[6][1];

  return Gauss_weights[order][point_numb];
}

//
// end
//
// /////////////////////////////////////////////////////////////////////////////


#endif

//////////////////////// Add sensitivity /////////////////////////////////////////////
int
AV3D4QuadWithSensitivity::setParameter(const char **argv, int argc, Parameter &param)
{


	int numberGauss=4;

    if (strstr(argv[0],"material") != 0) {
		int ok;
		//for (int i=0; i<numberGauss; i++) {
			ok = theMaterial->setParameter(&argv[1], argc-1, param);
			if (ok < 0 ) {
				//opserr<<"AV3D4QuadWithSensitivity::setParameter() can not setParameter for "<<i<<"th Gauss Point\n";
				return -1;
			//}
		}
		return ok;  
	}
	else{
		opserr<<"AV3D4QuadWithSensitivity can not setParameter!"<<endln;
    // otherwise parameter is unknown for the Truss class
	return -1;
	}


}
 int
AV3D4QuadWithSensitivity::updateParameter(int parameterID, Information &info)
{

  opserr<<"warnning: AV3D4QuadWithSensitivity can not updateParameter!"<<endln;
  return -1;
}

 int
AV3D4QuadWithSensitivity::activateParameter(int passedParameterID)
{

	int numberGauss = 4;
	parameterID = passedParameterID;

	if (parameterID == 1) {
		// Don't treat this case for now
	}

	else if (parameterID==0) {
		// Go down to the materials and zero out the identifier
		int ok;
		//for (int i=0; i<numberGauss; i++) {
			ok = theMaterial->activateParameter(parameterID);
			if (ok < 0 ) {
				return -1;  
			}
		//}
	}
	else if (parameterID > 100) {
		// In this case the parameter belongs to the material
		int ok;
		//for (int i=0; i<numberGauss; i++) {
			ok = theMaterial->activateParameter(parameterID-100);
			if (ok < 0 ) {
				return -1;
			}
		//}
	}
	else {
		opserr << "AV3D4QuadWithSensitivity::activateParameter() -- unknown parameter " << endln;
	}

	return 0;
}
 const Matrix &
AV3D4QuadWithSensitivity::getMassSensitivity(int gradNumber)
{
	static Matrix mass(4,4);
	mass.Zero();
	return mass; 
}

 const Vector &
AV3D4QuadWithSensitivity::getResistingForceSensitivity(int gradNumber)
 {
    

  static Vector PSensitivity(numDOF);
  PSensitivity.Zero();

  static Vector v(numDOF);

  static Vector  vSensitivity(numDOF);
  vSensitivity.Zero();

  int i;
  for(i = 0; i < nodes_in_elem; i++) {
    const Vector &vel = theNodes[i]->getTrialVel();
    v(i) = vel(0);
  }
//opserr<<"v is "<<v<<endln;
 /*
  vSensitivity(0) = theNodes[0]->getVelSensitivity(1,gradNumber);
  vSensitivity(1) = theNodes[1]->getVelSensitivity(1,gradNumber);
  vSensitivity(2) = theNodes[2]->getVelSensitivity(1,gradNumber);
  vSensitivity(3) = theNodes[3]->getVelSensitivity(1,gradNumber);
*/

  this->getDampSensitivity(gradNumber);

   // opserr<<"C  is "<<C<<endln;
    //opserr<<"CSensitivity is "<<CSensitivity<<endln;
    //opserr<<"v is "<<v<<endln;
   // opserr<<"vSensitivity is "<<vSensitivity<<endln;
  
  PSensitivity.addMatrixVector(1.0, CSensitivity, v, 1.0);
 // PSensitivity.addMatrixVector(1.0, C, vSensitivity, 1.0);
    
 // opserr<<"PSensitivity is "<<PSensitivity<<endln;

  return PSensitivity;
 
 }


 int
AV3D4QuadWithSensitivity::commitSensitivity(int gradNumber, int numGrads)
 {

   return 0;
 
 }
