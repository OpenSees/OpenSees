


// Written: Quan Gu, Yichao Gao and Zhijian Qiu  
// Created: 2015/01/25 
// Sensitivity analysis of absorbing-transmitting element for boundaries of solid
//------------------------------------------


#ifndef VS3D4QuadWithSensitivity_CPP
#define VS3D4QuadWithSensitivity_CPP

#include "VS3D4QuadWithSensitivity.h"
#include <math.h>
  
#define FixedOrder 2

// static variables
const int VS3D4QuadWithSensitivity::numDOF = 12;
const int VS3D4QuadWithSensitivity::nodes_in_elem = 4;
const int VS3D4QuadWithSensitivity::nodes_in_quad = 4;
const int VS3D4QuadWithSensitivity::r_integration_order = FixedOrder;
const int VS3D4QuadWithSensitivity::s_integration_order = FixedOrder;
const int VS3D4QuadWithSensitivity::dim = 3;
const int VS3D4QuadWithSensitivity::numGP = 4;

ID VS3D4QuadWithSensitivity::integFlags(0,numDOF);
ID VS3D4QuadWithSensitivity::actDOFs(0,numDOF);
  
Matrix VS3D4QuadWithSensitivity::K(numDOF,numDOF);
Matrix VS3D4QuadWithSensitivity::C(numDOF,numDOF);
Matrix VS3D4QuadWithSensitivity::M(numDOF,numDOF);
Vector VS3D4QuadWithSensitivity::P(numDOF);

Matrix ** VS3D4QuadWithSensitivity::H =0;
Matrix ** VS3D4QuadWithSensitivity::DH =0;
Matrix ** VS3D4QuadWithSensitivity::HH =0;


Matrix  VS3D4QuadWithSensitivity::mass(12,12) ;


#include <elementAPI.h>

void *
OPS_VS3D4WuadWithSensitivity(void) {

  double _rho = 1.0;
  double _R = 1.0;
  double _alphaN = 1.33;
  double _alphaT = 0.67;

  static int idData[5];
  static double dData[6];

  dData[2] = _rho;
  dData[3] = _R;
  dData[4] = _alphaN;
  dData[5] = _alphaT;

  //if the number of arguments is less than the minimum, throw an error
  int argc = OPS_GetNumRemainingInputArgs();
  if (argc < 9 || argc > 11) {
    opserr << "element Vs3D4 incorrect num args .. between 9 and 11 expected\n";
    return 0;
  }

  int numData = 5;
  if (OPS_GetIntInput(&numData, idData) != 0) {
    opserr << "element Vs3D4 error reading first 5 integers\n";
    return 0;
  }

  numData = argc-5;
  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "element Vs3D4 error reading last few doubles for element" << idData[0] << endln;
    return 0;
  }
  
  return new VS3D4QuadWithSensitivity(idData[0],idData[1],idData[2],idData[3],idData[4],
				      dData[0],dData[1],dData[2],dData[3],dData[4],dData[5]);
}

// 
// /////////////////////////////////////////////////////////////////////////////
//
// constructors and de-constructors
//

VS3D4QuadWithSensitivity::VS3D4QuadWithSensitivity(int element_number,
  int node_numb_1, int node_numb_2, int node_numb_3, int node_numb_4, 
  double _E, double _G, double _rho, double _R, double _alphaN, 
  double _alphaT)
:Element(element_number, ELE_TAG_VS3D4QuadWithSensitivity),
connectedExternalNodes(nodes_in_elem), Ki(0), E(_E), G(_G), rho(_rho), 
R(_R), alphaN(_alphaN), alphaT(_alphaT), area(0.0), NdotN(3,3), 
hasConstrained(0)
{
  // Set connected external node IDs
  connectedExternalNodes(0) = node_numb_1;
  connectedExternalNodes(1) = node_numb_2;
  connectedExternalNodes(2) = node_numb_3;
  connectedExternalNodes(3) = node_numb_4;
  
  // zero node pointers
  for (int i = 0; i < nodes_in_elem; i++)
    theNodes[i] = 0;

  // AddingSensitivity:BEGIN ////////////////////////////////
	parameterID = 0;
// AddingSensitivity:END /////////////////////////////////
}


VS3D4QuadWithSensitivity::VS3D4QuadWithSensitivity ()
:Element(0, ELE_TAG_VS3D4QuadWithSensitivity ),
connectedExternalNodes(nodes_in_elem), Ki(0), E(0.0), G(0.0), rho(1.0), 
R(1.0), alphaN(1.33), alphaT(0.67), area(0.0), NdotN(3,3), 
hasConstrained(0)
{
  // zero node pointers
  for (int i = 0; i < nodes_in_elem; i++)
    theNodes[i] = 0;
  // AddingSensitivity:BEGIN ////////////////////////////////
	parameterID = 0;
// AddingSensitivity:END /////////////////////////////////
}

VS3D4QuadWithSensitivity::~VS3D4QuadWithSensitivity ()
{
  if (Ki != 0)
    delete Ki;

}

//
// /////////////////////////////////////////////////////////////////////////////
//

VS3D4QuadWithSensitivity & 
VS3D4QuadWithSensitivity::operator[](int subscript)
{
  return ( *(this+subscript) );
}

int 
VS3D4QuadWithSensitivity::getNumExternalNodes () const
{
  return nodes_in_elem;
}

const ID& 
VS3D4QuadWithSensitivity::getExternalNodes ()
{
  return connectedExternalNodes;
}

Node ** 
VS3D4QuadWithSensitivity::getNodePtrs(void)
{
  return theNodes;
}

int 
VS3D4QuadWithSensitivity::getNumDOF ()
{
  return numDOF;
}

//
// /////////////////////////////////////////////////////////////////////////////
//

void 
VS3D4QuadWithSensitivity::setDomain (Domain *theDomain)
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
        opserr << "FATAL ERROR VS3D4QuadWithSensitivity (tag: " << this->getTag();
        opserr << " ), node not found in domain\n";
        exit(-1);
      }
    }
    
    this->DomainComponent::setDomain(theDomain);
  }
}

int
VS3D4QuadWithSensitivity::update()
{
  // do nothing
  return 0;
}

int 
VS3D4QuadWithSensitivity::commitState ()
{
  int retVal = 0;

  // call element commitState to do any base class stuff
  if ((retVal = this->Element::commitState()) != 0) {
    opserr << "VS3D4QuadWithSensitivity::commitState () - failed in base class";
  }
  
  return retVal;
}

int 
VS3D4QuadWithSensitivity::revertToLastCommit ()
{ 
  // do nothing
  return 0;
}

int 
VS3D4QuadWithSensitivity::revertToStart () 
{
  // do nothing
  return 0;
} 

//
// /////////////////////////////////////////////////////////////////////////////
//

void 
VS3D4QuadWithSensitivity::zeroLoad(void)
{
  // Q.Zero();
}

int 
VS3D4QuadWithSensitivity::addLoad(ElementalLoad *theLoad, double loadFactor)
{
  // do nothing
  return 0;
}

int 
VS3D4QuadWithSensitivity::addInertiaLoadToUnbalance(const Vector &accel)
{
  // do nothing
  return 0;
}

//
// /////////////////////////////////////////////////////////////////////////////
//

const Vector &
VS3D4QuadWithSensitivity::getResistingForce(void) 
{
 // P.Zero();
  
//  return P;

P.Zero();
  
  this->computeCoef();
  
  Vector subFD(3);
  Vector subFV(3);
  int i, pos;
  Matrix subMat(3,3);
  
  // Matrix h(1,nodes_in_quad);
  
  // normal and tangent spring coefficient
  double KN = alphaN*G/R*area*0.25;
  double KT = alphaT*G/R*area*0.25;
  
  // normal and tangent damping coefficient
  double CN = sqrt(E*rho)*area*0.25;
  double CT = sqrt(G*rho)*area*0.25;
  
  // int counter = 0;
  pos = 0;
  for(i = 0; i < nodes_in_elem; i++) {
    const Vector &disp = theNodes[i]->getTrialDisp();
    const Vector &vel = theNodes[i]->getTrialVel();
    
    subFD.addMatrixVector(0.0, NdotN, disp, KN-KT);
    subFD.addVector(1.0, disp, KT);
    
    subFV.addMatrixVector(0.0, NdotN, vel, CN-CT);
    subFV.addVector(1.0, vel, CT);
    
    // P(counter) = P(counter) + subFD(0) + subFV(0);
    // counter++;
    // P(counter) = P(counter) + subFD(1) + subFV(1);
    // counter++;
    // P(counter) = P(counter) + subFD(2) + subFV(2);
    // counter++;
    
    P.Assemble(subFD, pos, 1.0);
    P.Assemble(subFV, pos, 1.0);
    pos += 3;
    
  }
 
  return P;

}

const Vector &
VS3D4QuadWithSensitivity::getResistingForceIncInertia(void) 
{

  return this->getResistingForce();
  
}

//
// /////////////////////////////////////////////////////////////////////////////
//

int 
VS3D4QuadWithSensitivity::sendSelf (int commitTag, Channel &theChannel) 
{ 
  // Not implemtented yet
  return 0;
}

int 
VS3D4QuadWithSensitivity::recvSelf (int commitTag, Channel &theChannel, 
  FEM_ObjectBroker &theBroker) 
{   
  // Not implemtented yet
  return 0;
}

int VS3D4QuadWithSensitivity::displaySelf (Renderer &theViewer, int displayMode, float fact, const char **modes, int numMode)
{
  int error = 0;
  
  return error;

}     

void VS3D4QuadWithSensitivity::Print(OPS_Stream &s, int flag)
{
  if(flag == 1) {
    s << "VS3D4QuadWithSensitivity, element id:  " << this->getTag() << endln;
    s << "Connected external nodes:  " << connectedExternalNodes;
    s << this->getResistingForce();
  } 
  else {
    s << "VS3D4QuadWithSensitivity, element id:  " << this->getTag() << endln;
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
VS3D4QuadWithSensitivity::setResponse (const char **argv, int argc, OPS_Stream &output)
{
  Response *theResponse = 0;

  char outputData[32];

  output.tag("ElementOutput");
  output.attr("eleType","VS3D4QuadWithSensitivity");
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
VS3D4QuadWithSensitivity::getResponse (int responseID, Information &eleInfo)
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
VS3D4QuadWithSensitivity::getActiveDofs(void)
{
  if( actDOFs.Size() == 0){
    int counter = 0;
    for(int i = 0; i < nodes_in_elem; i++){
     actDOFs[counter++] = 1;
     actDOFs[counter++] = 2;
     actDOFs[counter++] = 3;
    }
  }

  return &actDOFs;
}

int 
VS3D4QuadWithSensitivity::getIntegrateFlag(void)
{
  return 0;
}

ID *
VS3D4QuadWithSensitivity::getIntegrateFlags(void)
{
  // if (integFlags.Size() == 0){
  //   for (int i = 0; i < nodes_in_elem; i++){
  //     integFlags[i] = 1;
  //   }
  // }

  return &integFlags;
}

//
// /////////////////////////////////////////////////////////////////////////////
//

Matrix 
VS3D4QuadWithSensitivity::getNodalCoords(void)
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
VS3D4QuadWithSensitivity::getTangentStiff(void)
{
  K.Zero();
  
  this->computeCoef();
  
  int i;
  Matrix subMat(3,3);
  
  // Matrix h(1,nodes_in_quad);
  
  // normal and tangent spring coefficient
  double KN = alphaN*G/R*area*0.25;
  double KT = alphaT*G/R*area*0.25;
  
  // subMat.Zero();
  // subMat.addMatrix(0.0, NdotN, KN);
  // subMat.addMatrix(1.0, NdotN, -KT);
  
  subMat.addMatrix(0.0, NdotN, KN-KT);
  subMat(0,0) = subMat(0,0) + KT;
  subMat(1,1) = subMat(1,1) + KT;
  subMat(2,2) = subMat(2,2) + KT;
  
  for(i = 0; i < 12; i += 3) {
    K.Assemble(subMat, i, i, 1.0);
  }
  
  return K;
}

const Matrix &
VS3D4QuadWithSensitivity::getInitialStiff(void)
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
VS3D4QuadWithSensitivity::getDamp(void)
{
  // zero [C] first
  C.Zero();
  
  this->computeCoef();
  
  int i;
  Matrix subMat(3,3);
  
  // normal and tangent damping coefficient
  double CN = sqrt(E*rho)*area*0.25;
  double CT = sqrt(G*rho)*area*0.25;
  
  // subMat.Zero();
  // subMat.addMatrix(0.0, NdotN, CN);
  // subMat.addMatrix(1.0, NdotN, -CT);
  
  subMat.addMatrix(0.0, NdotN, CN-CT);
  subMat(0,0) = subMat(0,0) + CT;
  subMat(1,1) = subMat(1,1) + CT;
  subMat(2,2) = subMat(2,2) + CT;
  
  for(i = 0; i < 12; i += 3) {
    C.Assemble(subMat, i, i, 1.0);
  }
  
  return C;
}

const Matrix &
VS3D4QuadWithSensitivity::getMass(void)
{
  // zero [M] first
  M.Zero();
  
  return M;
}

const Matrix &
VS3D4QuadWithSensitivity::getConsMass(void)
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
VS3D4QuadWithSensitivity::interp_fun(double r1, double r2)
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
VS3D4QuadWithSensitivity::diff_interp_fun(double r1, double r2)
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
VS3D4QuadWithSensitivity::computeH(void)
{
  if(H == 0 || DH == 0) {
    H = new Matrix*[numGP];
    DH = new Matrix*[numGP];
    if (H == 0 || DH == 0) {
      opserr << "VS3D4QuadWithSensitivity::computeH - out of memory!\n";
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
          opserr << "VS3D4QuadWithSensitivity::computeH - out of memory!\n";
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


// compute area and normal vector
// assume: four nodes are on the same plane

int 
VS3D4QuadWithSensitivity::computeCoef(void)
{
  if(area > 0.0) {
    return 0;
  }
  
  if(area < 0.0) area = 0.0;
  
  Matrix Jacobian(2,3);
  Matrix NC = getNodalCoords();
  
  this->computeH();
  
  short where = 0;
  double rw, sw;
  double x1, y1, z1;
  // double nx, ny, nz;  // normal vector
  double len;
  
  Matrix nVec(1,3);  // matrix hold normal vector
  
  
  // compute normal vector
  Jacobian = (*DH[0])*NC;
  
  x1 = Jacobian(0,1)*Jacobian(1,2) - Jacobian(0,2)*Jacobian(1,1);
  y1 = Jacobian(0,2)*Jacobian(1,0) - Jacobian(0,0)*Jacobian(1,2);
  z1 = Jacobian(0,0)*Jacobian(1,1) - Jacobian(0,1)*Jacobian(1,0);
  
  len = sqrt(x1*x1 + y1*y1 + z1*z1);
  if(len == 0.0){
    opserr <<"The length of tangent should not be 0!\n";
    exit(-1);
  }   
  
  // nx = x1/len;
  // ny = y1/len;
  // nz = z1/len;
  
  nVec(0,0) = x1/len;
  nVec(0,1) = y1/len;
  nVec(0,2) = z1/len;
  
  //printf("Normal vector is <%g,%g,%g>!\n", nVec(0,0), nVec(0,1), nVec(0,2));
  
  NdotN.addMatrixTransposeProduct(0.0, nVec, nVec, 1.0);
  
  for(short GP_c_r = 1; GP_c_r <= r_integration_order; GP_c_r++) {
    rw = get_Gauss_p_w(r_integration_order, GP_c_r);
    for(short GP_c_s = 1; GP_c_s <= s_integration_order; GP_c_s++) {
      sw = get_Gauss_p_w(s_integration_order, GP_c_s);
      
      Jacobian = (*DH[where])*NC;
      
      // pointer to the inward direction
      x1 = Jacobian(0,1)*Jacobian(1,2) - Jacobian(0,2)*Jacobian(1,1);
      y1 = Jacobian(0,2)*Jacobian(1,0) - Jacobian(0,0)*Jacobian(1,2);
      z1 = Jacobian(0,0)*Jacobian(1,1) - Jacobian(0,1)*Jacobian(1,0);
      
      // length of the normal vectio
      len = sqrt(x1*x1 + y1*y1 + z1*z1);
      if(len == 0.0){
        opserr <<"The length of tangent should not be 0!\n";
        exit(-1);
      }
      
      area += rw*sw*len;
      
      where++;
    }
  }
  
  return 0;
}






/*

int
VS3D4QuadWithSensitivity::computeHH(void)
{
  if (HH == 0) {
    HH = new Matrix*[numGP];
    if (HH == 0) {
      opserr << "VS3D4QuadWithSensitivity::computeHH - out of memory!\n";
      return -3;
    }
    
    // compute H first
    this->computeH();
    
    for(int i = 0; i < numGP; i++) {
      HH[i] = new Matrix(nodes_in_elem, nodes_in_elem);
      if (HH[i] == 0) {
        opserr << "VS3D4QuadWithSensitivity::computeHH - out of memory!\n";
        return -3;
      }
      
      HH[i]->addMatrixTransposeProduct(0.0, *H[i], *H[i], 1.0);
    }
  }
  
  return 0;
}

*/

/*
int 
VS3D4QuadWithSensitivity::computeDetJ(void)
{
  if(detJ == 0) {
    detJ = new double[numGP];
    if(detJ == 0) {
      opserr << "VS3D4QuadWithSensitivity::computeDetJ - out of memory!\n";
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

*/

/*
int
VS3D4QuadWithSensitivity::computeDiff(void)
{
  if (L == 0 || detJ == 0) {
    L = new Matrix*[numGP];
    detJ = new double[numGP];
    if (L == 0 || detJ == 0) {
      opserr << "VS3D4QuadWithSensitivity::computeDiff - out of memory!\n";
      return -3;
    }
    
    Matrix Jacobian(3,3);
    
    this->computeH();
    Matrix NC = getNodalCoords();
    
    for(int i = 0; i < numGP; i++) {
      L[i] = new Matrix(3, nodes_in_elem);
      if(L[i] == 0) {
        opserr << "VS3D4QuadWithSensitivity::computDiff() - out of memory!\n";
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
VS3D4QuadWithSensitivity::get_Gauss_p_c(short order, short point_numb)
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
VS3D4QuadWithSensitivity::get_Gauss_p_w(short order, short point_numb)
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
VS3D4QuadWithSensitivity::setParameter(const char **argv, int argc, Parameter &param)
{


	int numberGauss=4;

    if (strstr(argv[0],"material") != 0) {
		int ok;
		for (int i=0; i<numberGauss; i++) {
			ok = theMaterial[i]->setParameter(&argv[1], argc-1, param);
			if (ok < 0 ) {
				opserr<<"VS3D4QuadWithSensitivity::setParameter() can not setParameter for "<<i<<"th Gauss Point\n";
				return -1;
			}
		}
		return ok;
	}
	else{
		opserr<<"VS3D4QuadWithSensitivity can not setParameter!"<<endln;
    // otherwise parameter is unknown for the Truss class
	return -1;
	}


}
 int
VS3D4QuadWithSensitivity::updateParameter(int parameterID, Information &info)
{

  opserr<<"warnning: VS3D4QuadWithSensitivity can not updateParameter!"<<endln;
  return -1;
}

 int
VS3D4QuadWithSensitivity::activateParameter(int passedParameterID)
{

	int numberGauss = 4;
	parameterID = passedParameterID;

	if (parameterID == 1) {
		// Don't treat this case for now
	}

	else if (parameterID==0) {
		// Go down to the materials and zero out the identifier
		int ok;
		for (int i=0; i<numberGauss; i++) {
			ok = theMaterial[i]->activateParameter(parameterID);
			if (ok < 0 ) {
				return -1;
			}
		}
	}
	else if (parameterID > 100) {
		// In this case the parameter belongs to the material
		int ok;
		for (int i=0; i<numberGauss; i++) {
			ok = theMaterial[i]->activateParameter(parameterID-100);
			if (ok < 0 ) {
				return -1;
			}
		}
	}
	else {
		opserr << "VS3D4QuadWithSensitivity::activateParameter() -- unknown parameter " << endln;
	}

	return 0;
}
 const Matrix &
VS3D4QuadWithSensitivity::getMassSensitivity(int gradNumber)
{
	mass.Zero();
	return mass; 
}

 const Vector &
VS3D4QuadWithSensitivity::getResistingForceSensitivity(int gradNumber)
 {
 
  P.Zero();

  return P;
 
 
 }

 int
VS3D4QuadWithSensitivity::commitSensitivity(int gradNumber, int numGrads)
 {

   return 0;
 
 }
