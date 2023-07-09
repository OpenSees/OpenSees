
// Written: Quan Gu, Yichao Gao and Zhijian Qiu  
// Created: 2015/01/25 
// Sensitivity analysis of contact element for fluid-structure coupling
//------------------------------------------

#ifndef ASI3D8QuadWithSensitivity_CPP
#define ASI3D8QuadWithSensitivity_CPP

#include "ASI3D8QuadWithSensitivity.h"

#define FixedOrder 2

// static variables
const int ASI3D8QuadWithSensitivity::numDOF = 16;
const int ASI3D8QuadWithSensitivity::nodes_in_elem = 8;
const int ASI3D8QuadWithSensitivity::nodes_in_quad = 4;
const int ASI3D8QuadWithSensitivity::r_integration_order = FixedOrder;
const int ASI3D8QuadWithSensitivity::s_integration_order = FixedOrder;
const int ASI3D8QuadWithSensitivity::dim = 3;
const int ASI3D8QuadWithSensitivity::numGP = 4;
const int ASI3D8QuadWithSensitivity::numSDOF = 12;
const int ASI3D8QuadWithSensitivity::numFDOF = 4;

ID ASI3D8QuadWithSensitivity::integFlags(0,numDOF);
ID ASI3D8QuadWithSensitivity::actDOFs(0,numDOF);

Matrix ASI3D8QuadWithSensitivity::K(numDOF,numDOF);
Matrix ASI3D8QuadWithSensitivity::C(numDOF,numDOF);
Matrix ASI3D8QuadWithSensitivity::M(numDOF,numDOF);
Vector ASI3D8QuadWithSensitivity::P(numDOF);
Matrix ASI3D8QuadWithSensitivity::QMAT(numSDOF,numFDOF);

Matrix ** ASI3D8QuadWithSensitivity::H =0;
Matrix ** ASI3D8QuadWithSensitivity::DH =0;
// Matrix ** ASI3D8QuadWithSensitivity::HH =0;

Vector ASI3D8QuadWithSensitivity::VecS(numSDOF);  
Vector ASI3D8QuadWithSensitivity::VecF(numFDOF);

Matrix  ASI3D8QuadWithSensitivity::mass(16,16) ;


#include <elementAPI.h>

void *
OPS_ASID8QuadWithSensitivity(void){

  static int idData[9];

  //if the number of arguments is less than the minimum, throw an error
  int argc = OPS_GetNumRemainingInputArgs();
  if (argc != 9) {
    opserr << "element ASI3D8Quad incorrect num args .. 9 expected\n";
    return 0;
  }

  if (OPS_GetIntInput(&argc, idData) != 0) {
    opserr << "element ASI3D8Quad error reading first 9 integers\n";
    return 0;
  }  

  Element *theEle = new ASI3D8QuadWithSensitivity(idData[0],idData[1],idData[2],idData[3],idData[4],idData[5],idData[6],idData[7],idData[8]);
  return theEle;
}




// 
// /////////////////////////////////////////////////////////////////////////////
//
// constructors and de-constructors
//

ASI3D8QuadWithSensitivity::ASI3D8QuadWithSensitivity(int element_number,
  int node_numb_1, int node_numb_2, int node_numb_3, int node_numb_4, 
  int node_numb_5, int node_numb_6, int node_numb_7, int node_numb_8)
:Element(element_number, ELE_TAG_ASI3D8QuadWithSensitivity),
connectedExternalNodes(nodes_in_elem), Ki(0), hasConstrained(0)
{
  // Set connected external node IDs
  connectedExternalNodes(0) = node_numb_1;
  connectedExternalNodes(1) = node_numb_2;
  connectedExternalNodes(2) = node_numb_3;
  connectedExternalNodes(3) = node_numb_4;
  
  connectedExternalNodes(4) = node_numb_5;
  connectedExternalNodes(5) = node_numb_6;
  connectedExternalNodes(6) = node_numb_7;
  connectedExternalNodes(7) = node_numb_8;
  
  // zero node pointers
  for (int i = 0; i < nodes_in_elem; i++)
    theNodes[i] = 0;

  // AddingSensitivity:BEGIN ////////////////////////////////
	parameterID = 0;
// AddingSensitivity:END /////////////////////////////////

}

ASI3D8QuadWithSensitivity::ASI3D8QuadWithSensitivity ()
:Element(0, ELE_TAG_ASI3D8QuadWithSensitivity ),
connectedExternalNodes(nodes_in_elem), Ki(0), hasConstrained(0)
{
  // zero node pointers
  for (int i = 0; i < nodes_in_elem; i++)
    theNodes[i] = 0;
  // AddingSensitivity:BEGIN ////////////////////////////////
	parameterID = 0;
// AddingSensitivity:END /////////////////////////////////
}

ASI3D8QuadWithSensitivity::~ASI3D8QuadWithSensitivity ()
{
  if (Ki != 0)
    delete Ki;

}

//
// /////////////////////////////////////////////////////////////////////////////
//

ASI3D8QuadWithSensitivity & 
ASI3D8QuadWithSensitivity::operator[](int subscript)
{
  return ( *(this+subscript) );
}

int 
ASI3D8QuadWithSensitivity::getNumExternalNodes () const
{
  return nodes_in_elem;
}

const ID& 
ASI3D8QuadWithSensitivity::getExternalNodes ()
{
  return connectedExternalNodes;
}

Node ** 
ASI3D8QuadWithSensitivity::getNodePtrs(void)
{
  return theNodes;
}

int 
ASI3D8QuadWithSensitivity::getNumDOF ()
{
  return numDOF;
}

//
// /////////////////////////////////////////////////////////////////////////////
//

void 
ASI3D8QuadWithSensitivity::setDomain (Domain *theDomain)
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
        opserr << "FATAL ERROR ASI3D8QuadWithSensitivity (tag: " << this->getTag();
        opserr << " ), node not found in domain\n";
        exit(-1);
      }
    }
    
    this->DomainComponent::setDomain(theDomain);
  }
}

int
ASI3D8QuadWithSensitivity::update()
{
  // do nothing
  return 0;
}

int 
ASI3D8QuadWithSensitivity::commitState ()
{
  int retVal = 0;

  // call element commitState to do any base class stuff
  if ((retVal = this->Element::commitState()) != 0) {
    opserr << "ASI3D8QuadWithSensitivity::commitState () - failed in base class";
  }
  
  return retVal;
}

int 
ASI3D8QuadWithSensitivity::revertToLastCommit ()
{ 
  // do nothing
  return 0;
}

int 
ASI3D8QuadWithSensitivity::revertToStart () 
{
  // do nothing
  return 0;
} 

//
// /////////////////////////////////////////////////////////////////////////////
//

void 
ASI3D8QuadWithSensitivity::zeroLoad(void)
{
  // Q.Zero();
}

int 
ASI3D8QuadWithSensitivity::addLoad(ElementalLoad *theLoad, double loadFactor)
{
  // do nothing
  return 0;
}

int 
ASI3D8QuadWithSensitivity::addInertiaLoadToUnbalance(const Vector &accel)
{
  // do nothing
  return 0;
}

//
// /////////////////////////////////////////////////////////////////////////////
//

const Vector &
ASI3D8QuadWithSensitivity::getResistingForce(void) 
{
  P.Zero();
  
  Vector &pa = VecF;
  pa.Zero();
  
  int i, loc;
  loc = 0;
  for(i = nodes_in_quad; i < nodes_in_elem; i++) {
    const Vector &TotDisp = theNodes[i]->getTrialDisp();
    pa(loc++) = TotDisp(0);

  }
  
  Matrix matQ = getQMatrix();
  Vector Qp(numSDOF);
  
  Qp.addMatrixVector(0.0, matQ, pa, 1.0);

  for(i = 0; i < numSDOF; i++) {
    P(i) = Qp(i);
  }
  
  // P(0) = Qp(0);
  // P(1) = Qp(1);
  // P(3) = Qp(2);
  // P(4) = Qp(3);
  // P(6) = Qp(4);
  // P(7) = Qp(5);
  
  return P;
}

const Vector &
ASI3D8QuadWithSensitivity::getResistingForceIncInertia(void) 
{
  P.Zero();

  Vector &pa = VecF;
  Vector &a = VecS;
  
  pa.Zero();
  a.Zero();

  int i, counter;
  counter = 0;
  for(i = 0; i < nodes_in_quad; i++) {
    const Vector &accel = theNodes[i]->getTrialAccel();
    a(counter++) = accel(0);
    a(counter++) = accel(1);
    a(counter++) = accel(2);
  }

  counter = 0;
  for(i = nodes_in_quad; i < nodes_in_elem; i++) {
    const Vector &TotDisp = theNodes[i]->getTrialDisp();
    pa(counter++) = TotDisp(0);
	//opserr<<"TotDisp(0) is "<<TotDisp(0)<<endln;
  }
  
  Matrix matQ = getQMatrix();
  Vector Qp(numSDOF);
  Qp.addMatrixVector(0.0, matQ, pa, 1.0);
  
  Vector Qa(numFDOF);
  Qa.addMatrixTransposeVector(0.0, matQ, a, -1.0);
  
  for(i = 0; i < numSDOF; i++) {
    P(i) = Qp(i);
  }

  counter = 0;
  for(i = numSDOF; i < numDOF; i++) {
    P(i) = Qa(counter++);
  }
  
  // printf("Element tag is %d!\n", this->getTag());
  // for(i = 0; i < numDOF; i++) {
  //   printf("P(%d) = %g;\n", i+1, P(i));
  // }
 
  return P;
}

//
// /////////////////////////////////////////////////////////////////////////////
//

int 
ASI3D8QuadWithSensitivity::sendSelf (int commitTag, Channel &theChannel) 
{ 
  // Not implemtented yet
  return 0;
}

int 
ASI3D8QuadWithSensitivity::recvSelf (int commitTag, Channel &theChannel, 
  FEM_ObjectBroker &theBroker) 
{   
  // Not implemtented yet
  return 0;
}

int ASI3D8QuadWithSensitivity::displaySelf (Renderer &theViewer, int displayMode, float fact, const char **modes, int numMode)
{
  int error = 0;
  
  return error;

}     

void ASI3D8QuadWithSensitivity::Print(OPS_Stream &s, int flag)
{
  if(flag == 1) {
    s << "ASI3D8QuadWithSensitivity, element id:  " << this->getTag() << endln;
    s << "Connected external nodes:  " << connectedExternalNodes;
    s << this->getResistingForce();
  } 
  else {
    s << "ASI3D8QuadWithSensitivity, element id:  " << this->getTag() << endln;
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
ASI3D8QuadWithSensitivity::setResponse (const char **argv, int argc, OPS_Stream &output)
{
  Response *theResponse = 0;

  char outputData[32];

  output.tag("ElementOutput");
  output.attr("eleType","ASI3D8QuadWithSensitivity");
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
ASI3D8QuadWithSensitivity::getResponse (int responseID, Information &eleInfo)
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
ASI3D8QuadWithSensitivity::getActiveDofs(void)
{
  if( actDOFs.Size() == 0){
    int counter = 0;
    for (int i = 0; i < nodes_in_elem; i++){
     actDOFs[counter++] = 1;
     actDOFs[counter++] = 2;
     actDOFs[counter++] = 3;
     actDOFs[counter++] = 8;
    }
  }

  return &actDOFs;
}

int 
ASI3D8QuadWithSensitivity::getIntegrateFlag(void)
{
  return 2;
}

ID *
ASI3D8QuadWithSensitivity::getIntegrateFlags(void)
{
  if (integFlags.Size() == 0){
    int counter = 0;
    for (int i = 0; i < nodes_in_elem; i++){
      integFlags[counter++] = 0;
      integFlags[counter++] = 0;
      integFlags[counter++] = 0;
      integFlags[counter++] = 1;
    }
  }

  return &integFlags;
}

//
// /////////////////////////////////////////////////////////////////////////////
//

Matrix 
ASI3D8QuadWithSensitivity::getNodalCoords(void)
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
ASI3D8QuadWithSensitivity::getTangentStiff(void)
{
  Matrix &matQ = getQMatrix();
  K.Zero();
  
  // use matrix assemble
  ID rows(numSDOF);
  ID cols(numFDOF);
  
  int i, pos = 0;
  int counter = 0;
  for(i = 0; i < numSDOF; i++) {
    rows(i) = i;
  }
  
  counter = numSDOF;
  for(i = 0; i < numFDOF; i++) {
    cols(i) = counter;
    counter++;
  }
  
  // for(i = 0; i < nodes_in_elem; i++) {
  //   rows(pos++) = counter++;
  //   rows(pos++) = counter++;
  //   rows(pos++) = counter++;
  //   cols(i) = counter++;
  // }
  
  K.Assemble(matQ, rows, cols, 1.0);
  
  // print [K] matrix
  // static int kFlag = 0;
  // if(kFlag == 0) {
    // printf("Element tag is %d!\n", this->getTag());
    // for(int i = 0; i < numDOF; i++) {
    //   for(int j = 0; j < numDOF; j++) {
    //     printf("K(%d,%d) = %g;\n", i+1, j+1, K(i,j));
    //   }
    // }
    // kFlag++;
  // }
  
  return K;
}

const Matrix &
ASI3D8QuadWithSensitivity::getInitialStiff(void)
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
ASI3D8QuadWithSensitivity::getMass(void)
{
  Matrix &matQ = getQMatrix();
  // zero [M] first
  M.Zero();
  
  // use matrix assemble
  ID rows(numFDOF);
  ID cols(numSDOF);
  
  int i, j, pos = 0;
  int counter = 0;
  
  for(i = 0; i < numSDOF; i++) {
    cols(i) = i;
  }
  
  counter = numSDOF;
  for(i = 0; i < numFDOF; i++) {
    rows(i) = counter;
    counter++;
  }
  
  
  // for(i = 0; i < nodes_in_elem; i++) {
  //   cols(pos++) = counter++;
  //   cols(pos++) = counter++;
  //   cols(pos++) = counter++;
  //   rows(i) = counter++;
  // }
  
  int pos_Rows, pos_Cols;
  for(i = 0; i < numSDOF; i++) { // column number
    pos_Cols = cols(i);
    for(j = 0; j < numFDOF; j++) { // row number
      pos_Rows = rows(j);
      
      M(pos_Rows,pos_Cols) = -matQ(i,j);
    }
  }
  
  // printf [M] matrix
  // static int mFlag = 0;
  // if(mFlag == 0) {
    // printf("Element tag is %d!\n", this->getTag());
    // for(int i = 0; i < numDOF; i++) {
    //   for(int j = 0; j < numDOF; j++) {
    //     printf("M(%d,%d) = %g;\n", i+1, j+1, M(i,j));
    //   }
    // }
    // mFlag++;
  // }
  
  return M;
}

const Matrix &
ASI3D8QuadWithSensitivity::getConsMass(void)
{
  return M;
}

Matrix &
ASI3D8QuadWithSensitivity::getQMatrix(void)
{
  double r  = 0.0;
  double rw = 0.0;
  double s  = 0.0;
  double sw = 0.0;
  double x1, y1, z1;
  
  short where = 0;
  int i, counter;
  
  Matrix Jacobian(2,3);
  Matrix hdotN(numSDOF,1);
  
  // zero QMAT first
  QMAT.Zero();
  
  // get nodal coordinates
  Matrix N_C = getNodalCoords();
  
  this->computeH();
  
  for(short GP_c_r = 1; GP_c_r <= r_integration_order; GP_c_r++) {
    r = get_Gauss_p_c(r_integration_order, GP_c_r);
    rw = get_Gauss_p_w(r_integration_order, GP_c_r);
    for(short GP_c_s = 1; GP_c_s <= s_integration_order; GP_c_s++) {
      s = get_Gauss_p_c(s_integration_order, GP_c_s);
      sw = get_Gauss_p_w(s_integration_order, GP_c_s);
      
      Matrix &dh = *DH[where];
      Jacobian = dh*N_C;
      
      // pointer to the fluid
      x1 = Jacobian(0,1)*Jacobian(1,2) - Jacobian(0,2)*Jacobian(1,1);
      y1 = Jacobian(0,2)*Jacobian(1,0) - Jacobian(0,0)*Jacobian(1,2);
      z1 = Jacobian(0,0)*Jacobian(1,1) - Jacobian(0,1)*Jacobian(1,0);
      
      // printf("<x,y,z> = <%g,%g,%g>;\n", x1, y1, z1);
      
      Matrix &h = *H[where];
      
      counter = 0;
      for(i = 0; i < nodes_in_quad; i++) {
        hdotN(counter++,0) = h(0,i)*x1;
        hdotN(counter++,0) = h(0,i)*y1;
        hdotN(counter++,0) = h(0,i)*z1;
      }
      
      QMAT.addMatrixProduct(1.0, hdotN, h, rw*sw);
      
      where++;
    }
  }
  
  return QMAT;
}

//
// /////////////////////////////////////////////////////////////////////////////
//

// shape function
Matrix 
ASI3D8QuadWithSensitivity::interp_fun(double r1, double r2)
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
ASI3D8QuadWithSensitivity::diff_interp_fun(double r1, double r2)
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
ASI3D8QuadWithSensitivity::computeH(void)
{
  if(H == 0 || DH == 0) {
    H = new Matrix*[numGP];
    DH = new Matrix*[numGP];
    if (H == 0 || DH == 0) {
      opserr << "ASI3D8QuadWithSensitivity::computeH - out of memory!\n";
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
          opserr << "ASI3D8QuadWithSensitivity::computeH - out of memory!\n";
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

//
// /////////////////////////////////////////////////////////////////////////////
//




//
// /////////////////////////////////////////////////////////////////////////////
//

double 
ASI3D8QuadWithSensitivity::get_Gauss_p_c(short order, short point_numb)
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
ASI3D8QuadWithSensitivity::get_Gauss_p_w(short order, short point_numb)
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
ASI3D8QuadWithSensitivity::setParameter(const char **argv, int argc, Parameter &param)
{


	int numberGauss=8;

    if (strstr(argv[0],"material") != 0) {
		int ok;
		for (int i=0; i<numberGauss; i++) {
			ok = theMaterial[i]->setParameter(&argv[1], argc-1, param);
			if (ok < 0 ) {
				opserr<<"ASI3D8QuadWithSensitivityWithSensitivity::setParameter() can not setParameter for "<<i<<"th Gauss Point\n";
				return -1;
			}
		}
		return ok;
	}
	else{
		opserr<<"ASI3D8QuadWithSensitivityWithSensitivity can not setParameter!"<<endln;
    // otherwise parameter is unknown for the Truss class
	return -1;
	}


}
 int
ASI3D8QuadWithSensitivity::updateParameter(int parameterID, Information &info)
{

  opserr<<"warnning: ASI3D8QuadWithSensitivityWithSensitivity can not updateParameter!"<<endln;
  return -1;
}

 int
ASI3D8QuadWithSensitivity::activateParameter(int passedParameterID)
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
			ok = materialPointers[i]->activateParameter(parameterID-100);
			if (ok < 0 ) {
				return -1;
			}
		}
	}
	else {
		opserr << "ASI3D8QuadWithSensitivityWithSensitivity::activateParameter() -- unknown parameter " << endln;
	}

	return 0;
}

 const Matrix &
ASI3D8QuadWithSensitivity::getMassSensitivity(int gradNumber)
{
	mass.Zero();
	return mass; 
}

 const Vector &
ASI3D8QuadWithSensitivity::getResistingForceSensitivity(int gradNumber)
 {
 
  P.Zero();
 /* 
  Vector &pa = VecF;
  pa.Zero();
  
  int i, loc;
  loc = 0;


   static Vector TotDisp(4);
   TotDisp.Zero();	
		
   TotDisp(0) = theNodes[0]->getDispSensitivity(1,gradNumber);
   TotDisp(1) = theNodes[1]->getDispSensitivity(2,gradNumber);
   TotDisp(2) = theNodes[2]->getDispSensitivity(3,gradNumber);
   TotDisp(3) = theNodes[3]->getDispSensitivity(4,gradNumber);

   pa(0) = TotDisp(0);
   pa(1) = TotDisp(1);
   pa(2) = TotDisp(2);
   pa(3) = TotDisp(3);
 
  
  Matrix matQ = getQMatrix();
  Vector Qp(numSDOF);
  
  Qp.addMatrixVector(0.0, matQ, pa, 1.0);

  for(i = 0; i < numSDOF; i++) {
    P(i) = Qp(i);
  }
*/
  return P;
 
 
 }

 int
ASI3D8QuadWithSensitivity::commitSensitivity(int gradNumber, int numGrads)
 {

   return 0;
 
 }
