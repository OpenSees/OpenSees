
// Written: Quan Gu, Yichao Gao and Zhijian Qiu  
// Created: 2015/01/25 
// Fluid Element
//------------------------------------------

#ifndef AC3D8HexWithSensitivity_CPP
#define AC3D8HexWithSensitivity_CPP

// #include <DOF_Group.h>
#include "AC3D8HexWithSensitivity.h"
#include <math.h>
#include <Vector.h>
#include <Matrix.h>
  

#define FixedOrder 2

// static variables
const int AC3D8HexWithSensitivity::numDOF = 8;
const int AC3D8HexWithSensitivity::nodes_in_elem = 8;
const int AC3D8HexWithSensitivity::r_integration_order = FixedOrder;
const int AC3D8HexWithSensitivity::s_integration_order = FixedOrder;
const int AC3D8HexWithSensitivity::t_integration_order = FixedOrder;
const int AC3D8HexWithSensitivity::dim = 3;
const int AC3D8HexWithSensitivity::numGP = 8;


Matrix AC3D8HexWithSensitivity::K(numDOF,numDOF);
Matrix AC3D8HexWithSensitivity::C(numDOF,numDOF);
Matrix AC3D8HexWithSensitivity::M(numDOF,numDOF);
Vector AC3D8HexWithSensitivity::P(numDOF);
ID AC3D8HexWithSensitivity::actDOFs(0,numDOF);

Matrix ** AC3D8HexWithSensitivity::H =0;
Matrix ** AC3D8HexWithSensitivity::DH =0;
Matrix ** AC3D8HexWithSensitivity::HH =0;

Vector AC3D8HexWithSensitivity::VecA(numDOF);
Vector AC3D8HexWithSensitivity::VecV(numDOF);
Vector AC3D8HexWithSensitivity::VecD(numDOF);

Matrix  AC3D8HexWithSensitivity::mass(8,8) ;


#include <elementAPI.h>

void *
OPS_AC3D8HexWithSensitivity(void){

  int matTag;

  static int idData[10];

  //if the number of arguments is less than the minimum, throw an error
  int argc = OPS_GetNumRemainingInputArgs();
  if (argc != 10) {
    opserr << "element AC3D8Hex incorrect num args .. 11 expected\n";
    return 0;
  }

  if (OPS_GetIntInput(&argc, idData) != 0) {
    opserr << "element AC3D8Hex error reading integers\n";
    return 0;
  }  

  matTag = idData[9];
  NDMaterial *theMaterial = OPS_getNDMaterial(matTag);

  if (theMaterial == 0) {
    opserr << "command: element AC3D8Hex " << idData[0] << 
      " - no NDMaterial with tag " << matTag << " exists\n";
    return 0;      
  }
  
  Element *theEle = new AC3D8HexWithSensitivity(idData[0],idData[1],idData[2],idData[3],idData[4],idData[5],idData[6],idData[7],idData[8],theMaterial);
  return theEle;
}




// 
// /////////////////////////////////////////////////////////////////////////////
//
// constructors and de-constructors
//   

AC3D8HexWithSensitivity::AC3D8HexWithSensitivity(int element_number,
  int node_numb_1,  int node_numb_2,  int node_numb_3,  int node_numb_4,
  int node_numb_5,  int node_numb_6,  int node_numb_7,  int node_numb_8,
  NDMaterial * Globalmmodel)
  :Element(element_number, ELE_TAG_AC3D8HexWithSensitivity), hasConstrained(0), 
  connectedExternalNodes(nodes_in_elem), Ki(0), Q(numDOF), impVals(0), 
  L(0), detJ(0), theMaterial(0)
{
  // Set connected external node IDs 
  connectedExternalNodes( 0) = node_numb_1;
  connectedExternalNodes( 1) = node_numb_2;
  connectedExternalNodes( 2) = node_numb_3;
  connectedExternalNodes( 3) = node_numb_4;
  connectedExternalNodes( 4) = node_numb_5;
  connectedExternalNodes( 5) = node_numb_6;
  connectedExternalNodes( 6) = node_numb_7;
  connectedExternalNodes( 7) = node_numb_8;

  const char *type = Globalmmodel->getType();
  if(strcmp(type, "AcousticMedium") != 0){
    opserr << "AC3D8HexWithSensitivity::AC3D8HexWithSensitivity - incompatible material model\n";
    exit(-1);
  }
  
  // Allocate arrays of pointers to NDMaterials
  theMaterial = new NDMaterial *[numGP];
  if (theMaterial == 0) {
    opserr << "AC3D8HexWithSensitivity::AC3D8HexWithSensitivity - failed allocate material model pointer\n";
    exit(-1);
  }
  
  for (int i = 0; i < numGP; i++) {
    // Get copies of the material model for each integration point
    theMaterial[i] = Globalmmodel->getCopy();
      
    // Check allocation
    if (theMaterial[i] == 0) {
      opserr << "AC3D8HexWithSensitivity::AC3D8HexWithSensitivity -- failed to get a copy of material model\n";
      exit(-1);
    }
  }
  
  // zero node pointers
  for (int i = 0; i < nodes_in_elem; i++)
    theNodes[i] = 0;

  // AddingSensitivity:BEGIN ////////////////////////////////
 	parameterID = 0;  //1031
// AddingSensitivity:END /////////////////////////////////

}

AC3D8HexWithSensitivity::AC3D8HexWithSensitivity(int element_number,
  int node_numb_1,  int node_numb_2,  int node_numb_3,  int node_numb_4,
  int node_numb_5,  int node_numb_6,  int node_numb_7,  int node_numb_8)
  :Element(element_number, ELE_TAG_AC3D8HexWithSensitivity),
  connectedExternalNodes(nodes_in_elem), Ki(0), Q(numDOF), theMaterial(0), impVals(0), 
  L(0), detJ(0), hasConstrained(0)
{
  // Set connected external node IDs
  connectedExternalNodes( 0) = node_numb_1;
  connectedExternalNodes( 1) = node_numb_2;
  connectedExternalNodes( 2) = node_numb_3;
  connectedExternalNodes( 3) = node_numb_4;
  connectedExternalNodes( 4) = node_numb_5;
  connectedExternalNodes( 5) = node_numb_6;
  connectedExternalNodes( 6) = node_numb_7;
  connectedExternalNodes( 7) = node_numb_8;
  
  // zero node pointers
  for (int i = 0; i < nodes_in_elem; i++)
    theNodes[i] = 0;
}

AC3D8HexWithSensitivity::AC3D8HexWithSensitivity ()
  :Element(0, ELE_TAG_AC3D8HexWithSensitivity), connectedExternalNodes(nodes_in_elem), 
  Ki(0), Q(numDOF), impVals(0), L(0), detJ(0), hasConstrained(0), theMaterial(0)
{
  // zero node pointers
  for (int i = 0; i < nodes_in_elem; i++)
    theNodes[i] = 0;

  // AddingSensitivity:BEGIN ////////////////////////////////
  parameterID = 0; 
// AddingSensitivity:END /////////////////////////////////

}

AC3D8HexWithSensitivity::~AC3D8HexWithSensitivity ()
{
  if (Ki != 0)
    delete Ki;
  
  for (int i = 0; i < numGP; i++) {
    if (theMaterial[i])
      delete theMaterial[i];
    
    if (L[i] != 0)
      delete L[i];
    
  }
  
  if (impVals != 0)
    delete impVals;
  
  if (theMaterial)
    delete [] theMaterial;

  if (L != 0)
    delete L;
  
  if (detJ != 0)
    delete detJ;

}

//
// /////////////////////////////////////////////////////////////////////////////
//

AC3D8HexWithSensitivity & 
AC3D8HexWithSensitivity::operator[](int subscript)
{
  return ( *(this+subscript) );
}

int 
AC3D8HexWithSensitivity::getNumExternalNodes () const
{
  return nodes_in_elem;
}

const ID& 
AC3D8HexWithSensitivity::getExternalNodes ()
{
  return connectedExternalNodes;
}

Node ** 
AC3D8HexWithSensitivity::getNodePtrs(void)
{
  return theNodes;
}

int 
AC3D8HexWithSensitivity::getNumDOF ()
{
  return numDOF;
}

//
// /////////////////////////////////////////////////////////////////////////////
//

void 
AC3D8HexWithSensitivity::setDomain (Domain *theDomain)
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
        opserr << "FATAL ERROR AC3D8HexWithSensitivity (tag: " << this->getTag();
        opserr << " ), node not found in domain\n";
        exit(-1);
      }
    }
    
    this->DomainComponent::setDomain(theDomain);
  }
}

int
AC3D8HexWithSensitivity::update()
{
  int i;
  double r  = 0.0;
  double s  = 0.0;
  
  Vector epsilon(3);
  Matrix sstrain(3,1);
  
  Matrix trial_disp = this->getTotalDisp();
  
  //opserr<<"trial_disp"<<trial_disp<<endln;
  this->computeDiff();
  
  for(i = 0; i < numGP; i++) {
    const Matrix &dhGlobal = *L[i];
    
    sstrain.addMatrixProduct(0.0, dhGlobal, trial_disp, 1.0);
    epsilon(0) = sstrain(0,0);
    epsilon(1) = sstrain(1,0);
    epsilon(2) = sstrain(2,0);
    
    theMaterial[i]->setTrialStrain(epsilon);
    
    // printf("Integration point: %d\n", i+1);
    // printf("epsilon is <%g,%g,%g>\n", epsilon(0), epsilon(1), epsilon(2));
  }
  
  return 0;
}

int 
AC3D8HexWithSensitivity::commitState ()
{
  int retVal = 0;

  // call element commitState to do any base class stuff
  if ((retVal = this->Element::commitState()) != 0) {
    opserr << "AC3D8HexWithSensitivity::commitState () - failed in base class";
  }
  
  // Loop over the integration points and commit the material states
  for (int i = 0; i < numGP; i++)
    retVal += theMaterial[i]->commitState();
  
  return retVal;
}

int 
AC3D8HexWithSensitivity::revertToLastCommit ()
{ 
  int retVal = 0;

  // Loop over the integration points and revert to last committed material states
  for (int i = 0; i < numGP; i++)      
    retVal += theMaterial[i]->revertToLastCommit(); 
  
  return retVal;
}

int 
AC3D8HexWithSensitivity::revertToStart () 
{
  int retVal = 0;

  // Loop over the integration points and revert to initial material states
  for (int i = 0; i < numGP; i++)      
    retVal +=theMaterial[i]->revertToStart(); 
  
  return retVal;
} 

//
// /////////////////////////////////////////////////////////////////////////////
//

void 
AC3D8HexWithSensitivity::zeroLoad(void)
{
  Q.Zero();
}

int 
AC3D8HexWithSensitivity::addLoad(ElementalLoad *theLoad, double loadFactor)
{
  
  return 0;
}

int AC3D8HexWithSensitivity::addInertiaLoadToUnbalance(const Vector &accel)
{
  
  return 0;
}

//
// /////////////////////////////////////////////////////////////////////////////
//

const Vector &
AC3D8HexWithSensitivity::getResistingForce(void) 
{
  Matrix nodalforces = getNodalForces();
  
  int i, counter = 0;
  //converting nodalforce matrix to vector
  for (i = 0; i < nodes_in_elem; i++){
    P(i) = nodalforces(0,i);
  }

  // P.addVector(1.0, Q, -1.0);
  
  return P;
}

const Vector &AC3D8HexWithSensitivity::getResistingForceIncInertia(void) 
{
  // form resisting force
  // this->getResistingForce();
  
  int i;

  // printf("Resting-Force: Element tag is %d!\n", this->getTag());
  // for(i = 0; i < numDOF; i++) {
  //   printf("P(%d) = %g;\n", i+1, P(i));
  // }
  
  //
  // now add dynamic terms
  // P += M * a + C * v
  //
  
  Vector &a = VecA;
  Vector &v = VecV;
  Vector &d = VecD;
  
  a.Zero();
  v.Zero();
  
  this->getMass();
  this->getDamp();
  this->getTangentStiff();
  
  for(i = 0; i < nodes_in_elem; i++) {
    const Vector &accel = theNodes[i]->getTrialAccel();
    const Vector &vel = theNodes[i]->getTrialVel();
    const Vector &disp = theNodes[i]->getTrialDisp();
    
    a(i) = accel(0);
    v(i) = vel(0);
    d(i) = disp(0);
  }
  
  P.Zero();
  
  P.addMatrixVector(1.0, K, d, 1.0);
  
  // printf("Element tag is %d!\n", this->getTag());
  // for(i = 0; i < numDOF; i++) {
  //   printf("P(%d) = %g;\n", i+1, P(i));
  // }
  
  P.addMatrixVector(1.0, M, a, 1.0);
  // P.addMatrixVector(1.0, C, v, 1.0);
  
  // // add the damping forces if rayleigh damping
  // if (alphaM != 0.0 || betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0)
  //   P += this->getRayleighDampingForces();
  // 
  // } else {
  // 
  //   // add the damping forces if rayleigh damping
  //   if (betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0) 
  //     P += this->getRayleighDampingForces();
  
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
AC3D8HexWithSensitivity::sendSelf (int commitTag, Channel &theChannel) 
{ 
  // Not implemtented yet
  return 0;
}

int 
AC3D8HexWithSensitivity::recvSelf (int commitTag, Channel &theChannel, 
  FEM_ObjectBroker &theBroker) 
{   
  // Not implemtented yet
  return 0;
}

int 
AC3D8HexWithSensitivity::displaySelf (Renderer &theViewer, int displayMode, float fact, const char **modes, int numMode)
{
  int error  = 0;
  // Not implemtented yet
  return error;

}     

void AC3D8HexWithSensitivity::Print(OPS_Stream &s, int flag)
{
  if(flag == 1) {
    s << "AC3D8HexWithSensitivity, element id:  " << this->getTag() << endln;
    s << "Connected external nodes:  " << connectedExternalNodes;
    s << this->getResistingForce();
  } 
  else {
    s << "AC3D8HexWithSensitivity, element id:  " << this->getTag() << endln;
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
AC3D8HexWithSensitivity::setResponse (const char **argv, int argc, OPS_Stream &output)
{
  Response *theResponse = 0;

  char outputData[32];

  output.tag("ElementOutput");
  output.attr("eleType","AC3D8HexWithSensitivity");
  output.attr("eleTag",this->getTag());
  
  for(int i = 1; i <= nodes_in_elem; i++) {
    sprintf(outputData,"node%d",i);
    output.attr(outputData,theNodes[i-1]->getTag());
  }
  
  // // sample
  // if(_strcmpi(argv[0], "Stress") == 0) {
  //   theResponse =  new ElementResponse(this, 4, InfoS);
  // }
  
  output.endTag(); // ElementOutput
  return theResponse;
}

int 
AC3D8HexWithSensitivity::getResponse (int responseID, Information &eleInfo)
{
  switch(responseID) {
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
AC3D8HexWithSensitivity::getActiveDofs(void)
{
  if (actDOFs.Size() == 0){
    for (int i = 0; i < nodes_in_elem; i++){
     actDOFs[i] = 8;
    }
  }

  return &actDOFs;
}

int 
AC3D8HexWithSensitivity::setNDMaterial(NDMaterial *Globalmmodel)
{
  // printf("AC3D8HexWithSensitivity::setNDMaterial - tag is %d!\n", this->getTag());
  if(theMaterial != 0) {
    printf("AC3D8HexWithSensitivity -- NDMaterial has been set!\n");
    return -1;
  }
  
  const char *type = Globalmmodel->getType();
  if (strcmp(type, "AcousticMedium") != 0){
    opserr << "AC3D8HexWithSensitivity::setNDMaterial - incompatible material model\n";
    return -4;
  }
  
  theMaterial = new NDMaterial*[numGP];
  if (theMaterial == 0) {
    opserr << "AC3D8HexWithSensitivity::setNDMaterial - out of memory!\n";
    return -2;
  }
  
  for (int i = 0; i < numGP; i++) {
    theMaterial[i] = Globalmmodel->getCopy();
    if (theMaterial[i] == 0) {
      opserr << "AC3D8HexWithSensitivity::setNDMaterial -- failed to get a copy of material model\n";
      return -3;
    }
  }
  
  return 0;
}

int 
AC3D8HexWithSensitivity::getIntegrateFlag(void)
{
  return 1;
}

//
// /////////////////////////////////////////////////////////////////////////////
//

Matrix 
AC3D8HexWithSensitivity::getNodalCoords(void)
{
  Matrix N_Coord(nodes_in_elem,dim);
  
  for(int i = 0; i < nodes_in_elem; i++) {
    const Vector &ndCrds = theNodes[i]->getCrds();
    N_Coord(i,0) = ndCrds(0);
    N_Coord(i,1) = ndCrds(1);
    N_Coord(i,2) = ndCrds(2);
  }

  return N_Coord;
}

Matrix 
AC3D8HexWithSensitivity::getTotalDisp(void)
{
  Matrix total_disp(nodes_in_elem, 1);

  int i;
  for(i = 0; i < nodes_in_elem; i++) {
    const Vector &TotDis = theNodes[i]->getTrialDisp();
    total_disp(i,0) = TotDis(0);
  }
  
  // for(i = 0; i < nodes_in_elem; i++) {
  //   printf("total_disp(%d) = %g\n", i+1, total_disp(i,0));
  // }

  return total_disp;
}


Matrix  
AC3D8HexWithSensitivity::getNodalForces(void)
{
  double r  = 0.0;
  double rw = 0.0;
  double s  = 0.0;
  double sw = 0.0;
  double t  = 0.0;
  double tw = 0.0;
  double weight = 0.0;
  double det_of_Jacobian = 0.0;
  
  short where = 0;
  
  Matrix sigma(1,3);
  Matrix NF(1,nodes_in_elem);
  
  this->computeDiff();
  NF.Zero();
  
  for(short GP_c_r = 1; GP_c_r <= r_integration_order; GP_c_r++) {
    r = get_Gauss_p_c(r_integration_order, GP_c_r);
    rw = get_Gauss_p_w(r_integration_order, GP_c_r);
    for( short GP_c_s = 1; GP_c_s <= s_integration_order; GP_c_s++ ) {
      s = get_Gauss_p_c(s_integration_order, GP_c_s);
      sw = get_Gauss_p_w(s_integration_order, GP_c_s);
      for(short GP_c_t = 1; GP_c_t <= t_integration_order; GP_c_t++) {
        t = get_Gauss_p_c(t_integration_order, GP_c_t);
        tw = get_Gauss_p_w(t_integration_order, GP_c_t);
        
        det_of_Jacobian = detJ[where];
        Matrix &dhGlobal = *L[where];
        
        weight = sw * rw * tw * det_of_Jacobian;
        
        // assemble stress tensor
        const Vector &stressvec = theMaterial[where]->getStress();
        sigma(0,0) = stressvec(0);
        sigma(0,1) = stressvec(1);
        sigma(0,2) = stressvec(2);
        
        NF.addMatrixProduct(1.0, sigma, dhGlobal, weight);
        
        // nodal forces See Zienkievicz part 1 pp 108
        // nodal_forces = nodal_forces + dhGlobal("ib")*stress_at_GP("ab")*weight;
        
        where++;
      }
    }
  }
  
  return NF;
}


//
// /////////////////////////////////////////////////////////////////////////////
//
// get stiffness matrix and mass matrix
//

const Matrix &
AC3D8HexWithSensitivity::getTangentStiff(void)
{
  double r  = 0.0;
  double rw = 0.0;
  double s  = 0.0;
  double sw = 0.0;
  double t  = 0.0;
  double tw = 0.0;
  double weight = 0.0;
  double det_of_Jacobian = 0.0;
  
  short where = 0;
  
  this->computeDiff();
  
  // zero [K] first
  K.Zero();
  
  // get mass density of fluid
  double RHO = theMaterial[0]->getRho();
  if(RHO == 0.0){
    opserr << "ERROR: The mass density is zero!\n";
    exit(-1);
  }
  
  for(short GP_c_r = 1; GP_c_r <= r_integration_order; GP_c_r++) {
    r = get_Gauss_p_c(r_integration_order, GP_c_r);
    rw = get_Gauss_p_w(r_integration_order, GP_c_r);
    for(short GP_c_s = 1; GP_c_s <= s_integration_order; GP_c_s++) {
      s = get_Gauss_p_c(s_integration_order, GP_c_s);
      sw = get_Gauss_p_w(s_integration_order, GP_c_s);
      for(short GP_c_t = 1; GP_c_t <= t_integration_order; GP_c_t++) {
        t = get_Gauss_p_c(t_integration_order, GP_c_t);
        tw = get_Gauss_p_w(t_integration_order, GP_c_t);
        
        det_of_Jacobian = detJ[where];
        Matrix &dhGlobal = *L[where];
        weight = rw * sw * tw * det_of_Jacobian / RHO;
        
        K.addMatrixTransposeProduct(1.0, dhGlobal, dhGlobal, weight);
        
        where++;
      }
    }
  }
  
  // printf [K] matrix
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
AC3D8HexWithSensitivity::getInitialStiff(void)
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
AC3D8HexWithSensitivity::getMass(void)
{
  double r  = 0.0;
  double rw = 0.0;
  double s  = 0.0;
  double sw = 0.0;
  double t  = 0.0;
  double tw = 0.0;
  double weight = 0.0;
  double RHO = 0.0;
  
  short where = 0;
  
  // zero [M] first
  M.Zero();
  
  const Matrix &D = theMaterial[0]->getTangent();
  if(D(0,0) == 0.0){
    opserr << "ERROR: The Kf is zero!\n";
    exit(-1);
  }    
  
  this->computeHH();
  
  for(short GP_c_r = 1; GP_c_r <= r_integration_order; GP_c_r++) {
    r = get_Gauss_p_c(r_integration_order, GP_c_r);
    rw = get_Gauss_p_w(r_integration_order, GP_c_r);
    for(short GP_c_s = 1; GP_c_s <= s_integration_order; GP_c_s++) {
      s = get_Gauss_p_c(s_integration_order, GP_c_s);
      sw = get_Gauss_p_w(s_integration_order, GP_c_s);
      for(short GP_c_t = 1; GP_c_t <= t_integration_order; GP_c_t++) {
        t = get_Gauss_p_c(t_integration_order, GP_c_t);
        tw = get_Gauss_p_w(t_integration_order, GP_c_t);
        
        weight = rw * sw * tw * detJ[where] / D(0,0);
        Matrix &hh = *HH[where];
        M.addMatrix(1.0, hh, weight);
        
        where++;
      }
    }
  }
  
  // lumped mass
  int i, j;
  for(i = 0; i < numDOF; i++) {
    weight = 0.0;
    for(j = 0; j < numDOF; j++) {
      weight += M(i,j);
      M(i,j) = 0.0;
    }
    M(i,i) =weight;
  }
  
  
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
AC3D8HexWithSensitivity::getConsMass(void)
{
  return M;
}

const Matrix &
AC3D8HexWithSensitivity::getDamp(void)
{
  C.Zero();
  if (impVals == 0) {
    return C;
  }
  
  int i, nodes_in_face = 8;
  ID face_nodes(nodes_in_face);
  Matrix Cf(nodes_in_face, nodes_in_face);
  
  for(i = 1; i <= 6; i++) {
    if(impVals[i-1] != 0.0) {
      Cf = get_face_impedance(i);
      localFaceMapping(i, face_nodes);
      
      if(impVals[i-1] != 1.0) {
        Cf = Cf*impVals[i-1];
      }
      
      C.Assemble(Cf, face_nodes, face_nodes);
    }
  }
  
  return C;
}


// Matrix &
// AC3D8HexWithSensitivity::getBMatrix(Matrix dhGlobal)
// {
//   B.Zero();
//   
//   int i, counter = 0;
//   for(i = 0; i < nodes_in_elem; i++) {
//     B(0,counter) = dhGlobal(0,i);
//     B(2,counter) = dhGlobal(1,i);
//     counter++;
//     B(1,counter) = dhGlobal(1,i);
//     B(2,counter) = dhGlobal(0,i);
//     counter++;
//   }
//   
//   return B;
// }

//
// /////////////////////////////////////////////////////////////////////////////
//

// shape function
Matrix 
AC3D8HexWithSensitivity::interp_fun(double r1, double r2, double r3)
{
  Matrix h(1,nodes_in_elem);

  // influence of the node number 8
  h(0,7)=(1.0-r1)*(1.0+r2)*(1.0+r3)*0.125;
  // influence of the node number 7
  h(0,6)=(1.0+r1)*(1.0+r2)*(1.0+r3)*0.125;
  // influence of the node number 6
  h(0,5)=(1.0+r1)*(1.0-r2)*(1.0+r3)*0.125;
  // influence of the node number 5
  h(0,4)=(1.0-r1)*(1.0-r2)*(1.0+r3)*0.125;

  // influence of the node number 4
  h(0,3)=(1.0-r1)*(1.0+r2)*(1.0-r3)*0.125;
  // influence of the node number 3
  h(0,2)=(1.0+r1)*(1.0+r2)*(1.0-r3)*0.125;
  // influence of the node number 2
  h(0,1)=(1.0+r1)*(1.0-r2)*(1.0-r3)*0.125;
  // influence of the node number 1
  h(0,0)=(1.0-r1)*(1.0-r2)*(1.0-r3)*0.125;
  
  return h;
}

Matrix 
AC3D8HexWithSensitivity::diff_interp_fun(double r1, double r2, double r3)
{
  Matrix dh(3,nodes_in_elem);

  // influence of the node number 8
  //dh(0,7)= (1.0-r2)*(1.0-r3)/8.0;
  dh(0,7)=-(1.0+r2)*(1.0+r3)*0.125;
  dh(1,7)= (1.0-r1)*(1.0+r3)*0.125;
  dh(2,7)= (1.0-r1)*(1.0+r2)*0.125;
  // influence of the node number 7
  dh(0,6)= (1.0+r2)*(1.0+r3)*0.125;
  dh(1,6)= (1.0+r1)*(1.0+r3)*0.125;
  dh(2,6)= (1.0+r1)*(1.0+r2)*0.125;
  // influence of the node number 6
  dh(0,5)= (1.0-r2)*(1.0+r3)*0.125;
  dh(1,5)=-(1.0+r1)*(1.0+r3)*0.125;
  dh(2,5)= (1.0+r1)*(1.0-r2)*0.125;
  // influence of the node number 5
  dh(0,4)=-(1.0-r2)*(1.0+r3)*0.125;
  dh(1,4)=-(1.0-r1)*(1.0+r3)*0.125;
  dh(2,4)= (1.0-r1)*(1.0-r2)*0.125;

  // influence of the node number 4
  dh(0,3)=-(1.0+r2)*(1.0-r3)*0.125;
  dh(1,3)= (1.0-r1)*(1.0-r3)*0.125;
  dh(2,3)=-(1.0-r1)*(1.0+r2)*0.125;
  // influence of the node number 3
  dh(0,2)= (1.0+r2)*(1.0-r3)*0.125;
  dh(1,2)= (1.0+r1)*(1.0-r3)*0.125;
  dh(2,2)=-(1.0+r1)*(1.0+r2)*0.125;
  // influence of the node number 2
  dh(0,1)= (1.0-r2)*(1.0-r3)*0.125;
  dh(1,1)=-(1.0+r1)*(1.0-r3)*0.125;
  dh(2,1)=-(1.0+r1)*(1.0-r2)*0.125;
  // influence of the node number 1
  dh(0,0)=-(1.0-r2)*(1.0-r3)*0.125;
  dh(1,0)=-(1.0-r1)*(1.0-r3)*0.125;
  dh(2,0)=-(1.0-r1)*(1.0-r2)*0.125;
  
  return dh;
}

int
AC3D8HexWithSensitivity::computeH(void)
{
  if (H == 0 || DH == 0) {
    H = new Matrix*[numGP];
    DH = new Matrix*[numGP];
    if (H == 0 || DH == 0) {
      opserr << "AC3D8HexWithSensitivity::computeH - out of memory!\n";
      return -3;
    }

    double r = 0.0;
    double s = 0.0;
    double t = 0.0;
    short where = 0;
    
    for(short GP_c_r = 1; GP_c_r <= r_integration_order; GP_c_r++) {
      r = get_Gauss_p_c(r_integration_order, GP_c_r);
      for(short GP_c_s = 1; GP_c_s <= s_integration_order; GP_c_s++) {
        s = get_Gauss_p_c(s_integration_order, GP_c_s);
        for(short GP_c_t = 1; GP_c_t <= t_integration_order; GP_c_t++) {
          t = get_Gauss_p_c(t_integration_order, GP_c_t);
          
          H[where] = new Matrix(1, nodes_in_elem);
          DH[where] = new Matrix(dim, nodes_in_elem);
          if(H[where] == 0 || DH[where] == 0) {
            opserr << "AC3D8HexWithSensitivity::computeH - out of memory!\n";
            return -3;
          }
          
          *H[where] = interp_fun(r, s, t);
          *DH[where] = diff_interp_fun(r, s, t);
          
          where++;
        } // for GP_c_t
      } // for GP_c_s
    } // for GP_c_r
  }
  
  return 0;
}

int
AC3D8HexWithSensitivity::computeHH(void)
{
  if (HH == 0) {
    HH = new Matrix*[numGP];
    if (HH == 0) {
      opserr << "AC3D8HexWithSensitivity::computeHH - out of memory!\n";
      return -3;
    }
    
    // compute H first
    this->computeH();
    
    for(int i = 0; i < numGP; i++) {
      HH[i] = new Matrix(nodes_in_elem, nodes_in_elem);
      if (HH[i] == 0) {
        opserr << "AC3D8HexWithSensitivity::computeHH - out of memory!\n";
        return -3;
      }
      
      HH[i]->addMatrixTransposeProduct(0.0, *H[i], *H[i], 1.0);
    }
  }
  
  
  return 0;
}

int
AC3D8HexWithSensitivity::computeDiff(void)
{
  if (L == 0 || detJ == 0) {
    L = new Matrix*[numGP];
    detJ = new double[numGP];
    if (L == 0 || detJ == 0) {
      opserr << "AC3D8HexWithSensitivity::computeDiff - out of memory!\n";
      return -3;
    }
    
    Matrix Jacobian(3,3);
    
    this->computeH();
    Matrix NC = getNodalCoords();
    
    for(int i = 0; i < numGP; i++) {
      L[i] = new Matrix(3, nodes_in_elem);
      if(L[i] == 0) {
        opserr << "AC3D8HexWithSensitivity::computDiff() - out of memory!\n";
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

double 
AC3D8HexWithSensitivity::Jacobian_det(Matrix Jac)
{
  // for 3*3 matrix, detJ = J11*J22*J33 + J12*J23*J31 + J13*J21*J32
  // - J13*J22*J31 - J12*J21*J33 - J11*J23*J32
  double res = Jac(0,0)*Jac(1,1)*Jac(2,2);
  res += Jac(0,1)*Jac(1,2)*Jac(2,0);
  res += Jac(0,2)*Jac(1,0)*Jac(2,1);
  res -= Jac(0,2)*Jac(1,1)*Jac(2,0);
  res -= Jac(0,1)*Jac(1,0)*Jac(2,2);
  res -= Jac(0,0)*Jac(1,2)*Jac(2,1);
  
  return res;
}

//
// /////////////////////////////////////////////////////////////////////////////
//

Matrix 
AC3D8HexWithSensitivity::interp_fun_face(double r1, double r2)
{
  Matrix h(1,8);
  
  // influence of the node number 8
  h(0,7)=(1.0-r1)*(1.0-r2*r2)*0.5;
  // influence of the node number 7
  h(0,6)=(1.0-r1*r1)*(1.0+r2)*0.5;
  // influence of the node number 6
  h(0,5)=(1.0+r1)*(1.0-r2*r2)*0.5;
  // influence of the node number 5
  h(0,4)=(1.0-r1*r1)*(1.0-r2)*0.5;

  // influence of the node number 4
  h(0,3)=(1.0-r1)*(1.0+r2)*0.25 - (h(0,6)+h(0,7))*0.5;
  // influence of the node number 3
  h(0,2)=(1.0+r1)*(1.0+r2)*0.25 - (h(0,5)+h(0,6))*0.5;
  // influence of the node number 2
  h(0,1)=(1.0+r1)*(1.0-r2)*0.25 - (h(0,4)+h(0,5))*0.5;
  // influence of the node number 1
  h(0,0)=(1.0-r1)*(1.0-r2)*0.25 - (h(0,4)+h(0,7))*0.5;

  return h;
}


Matrix 
AC3D8HexWithSensitivity::diff_interp_fun_face(double r1, double r2)
{
  Matrix dh(2,8);
  
  // influence of the node number 8
  dh(0,7)=-(1.0-r2*r2)*0.5;
  dh(1,7)=-(1.0-r1)*r2;
  // influence of the node number 7
  dh(0,6)=-r1*(1.0+r2);
  dh(1,6)= (1.0-r1*r1)*0.5;
  // influence of the node number 6
  dh(0,5)=-dh(0,7); // (1.0-r2*r2)*0.5;
  dh(1,5)=-(1.0+r1)*r2;
  // influence of the node number 5
  dh(0,4)=-r1*(1.0-r2);
  dh(1,4)=-dh(1,6); //-(1.0-r1*r1)*0.5;

  // influence of the node number 4
  dh(0,3)=-(1.0+r2)*0.25 - (dh(0,6)+dh(0,7))*0.5;
  dh(1,3)= (1.0-r1)*0.25 - (dh(1,6)+dh(1,7))*0.5;
  // influence of the node number 3
  dh(0,2)= (1.0+r2)*0.25 - (dh(0,5)+dh(0,6))*0.5;
  dh(1,2)= (1.0+r1)*0.25 - (dh(1,5)+dh(1,6))*0.5;
  // influence of the node number 2
  dh(0,1)= (1.0-r2)*0.25 - (dh(0,4)+dh(0,5))*0.5;
  dh(1,1)=-(1.0+r1)*0.25 - (dh(1,4)+dh(1,5))*0.5;
  // influence of the node number 1
  dh(0,0)=-(1.0-r2)*0.25 - (dh(0,4)+dh(0,7))*0.5;
  dh(1,0)=-(1.0-r1)*0.25 - (dh(1,4)+dh(1,7))*0.5;
  
  return dh;
}


Matrix 
AC3D8HexWithSensitivity::getFaceNodalCoords(int face_num)
{
  int nodes_in_face = 8;
  Matrix N_coord(nodes_in_face,dim);
  
  if (face_num < 1 || face_num > 6) {
    opserr << "invalid face number!\n";
    return N_coord;
  }

  ID face_nodes(nodes_in_face);
  localFaceMapping(face_num, face_nodes);
  
  int i;
  for(i = 0; i < nodes_in_face; i++) {
    const Vector &ndCrds = theNodes[face_nodes(i)]->getCrds();
    
    N_coord(i,0) = ndCrds(0);
    N_coord(i,1) = ndCrds(1);
    N_coord(i,2) = ndCrds(2);
  }
  
  return N_coord;
}

void 
AC3D8HexWithSensitivity::localFaceMapping(int face_num, ID &face_nodes)
{
  if(face_num == 1) {
    face_nodes(0) =  0; face_nodes(1) =  1; face_nodes(2) =  2; face_nodes(3) =  3;
    face_nodes(4) =  8; face_nodes(5) =  9; face_nodes(6) = 10; face_nodes(7) = 11;
  }
  else if(face_num == 2) {
    face_nodes(0) =  4; face_nodes(1) =  7; face_nodes(2) =  6; face_nodes(3) =  5;
    face_nodes(4) = 15; face_nodes(5) = 14; face_nodes(6) = 13; face_nodes(7) = 12;
  }
  else if(face_num == 3) {
    face_nodes(0) =  0; face_nodes(1) =  4; face_nodes(2) =  5; face_nodes(3) =  1;
    face_nodes(4) = 16; face_nodes(5) = 12; face_nodes(6) = 17; face_nodes(7) =  8;
  }
  else if(face_num == 4) {
    face_nodes(0) =  1; face_nodes(1) =  5; face_nodes(2) =  6; face_nodes(3) =  2;
    face_nodes(4) = 17; face_nodes(5) = 13; face_nodes(6) = 18; face_nodes(7) =  9;
  }
  else if(face_num == 5) {
    face_nodes(0) =  2; face_nodes(1) =  6; face_nodes(2) =  7; face_nodes(3) =  3;
    face_nodes(4) = 18; face_nodes(5) = 14; face_nodes(6) = 19; face_nodes(7) = 10;
  }
  else if(face_num == 6) {
    face_nodes(0) =  3; face_nodes(1) =  7; face_nodes(2) =  4; face_nodes(3) =  0;
    face_nodes(4) = 19; face_nodes(5) = 15; face_nodes(6) = 16; face_nodes(7) = 11;
  }
}


int 
AC3D8HexWithSensitivity::setImpedance(int face, double val)
{
  if (face < 1 || face > 6) {
    printf("AC3D8HexWithSensitivity::setImpedance - invalid face number %d!\n", face);
    return -2;
  }
  
  if (val == 0.0) {
    return 0;
  }
  
  if (impVals == 0) {
    impVals = new double[6];
    if (impVals == 0) {
      printf("AC3D8HexWithSensitivity::setImpedance - out of memory\n");
      return -3;
    }
    impVals[0] = 0.0;
    impVals[1] = 0.0;
    impVals[2] = 0.0;
    impVals[3] = 0.0;
    impVals[4] = 0.0;
    impVals[5] = 0.0;
  }
  
  impVals[face-1] = val;
  
  return 0;
}

Matrix 
AC3D8HexWithSensitivity::get_face_impedance(int face_num)
{
  double r  = 0.0;
  double rw = 0.0;
  double s  = 0.0;
  double sw = 0.0;

  double weight;
  double RHO, Kf, cs;
  double x1,y1,z1,area;
  
  int nodes_in_face = 8;
  Matrix Cf(nodes_in_face, nodes_in_face);
  Matrix Jacobian(2,3);
  Matrix dh(2,nodes_in_face);
  Matrix h(1,nodes_in_face);
  
  Matrix N_coord = getFaceNodalCoords(face_num);
  
  // Get mass density of fluid
  RHO = theMaterial[0]->getRho();
  if(RHO == 0.0){
    opserr << "ERROR: The mass density is zero!\n";
    exit(-1);
  }
  Kf = (theMaterial[0]->getTangent())(0,0);
  cs = sqrt(Kf/RHO);
  
  // zero matrix first
  Cf.Zero();
  
  for(short GP_c_r = 1; GP_c_r <= r_integration_order; GP_c_r++) {
    r = get_Gauss_p_c(r_integration_order, GP_c_r);
    rw = get_Gauss_p_w(r_integration_order, GP_c_r);
    for(short GP_c_s = 1; GP_c_s <= s_integration_order; GP_c_s++) {
      s = get_Gauss_p_c(s_integration_order, GP_c_s);
      sw = get_Gauss_p_w(s_integration_order, GP_c_s);
      
      dh = diff_interp_fun_face(r, s);
      Jacobian = dh*N_coord;
      
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
      
      h = interp_fun_face(r, s);
      weight = rw * sw* area / RHO /cs;
      Cf.addMatrixTransposeProduct(1.0, h, h, weight);
    }
  }
  
  return Cf;
}

//
// /////////////////////////////////////////////////////////////////////////////
//

/*
void 
AC3D8HexWithSensitivity::addHydroStaticLoad(int face, double mag, double top, 
  double bot, int varCrd)
{
  if (face < 1 || face > 6) {
    opserr << "AC3D8HexWithSensitivity::addHydroStaticLoad - invalid face number:";
    opserr << face << "\n";
    return ;
  }
  
  // see if we can quick return
  if (mag == 0.0) 
    return ;
  
  double r, rw;
  double s, sw;
  int nodes_in_face = 8;
  double pCrd, pMag;
  double weight = 0.0;
  double x1, y1, len;
  
  Matrix h(1,nodes_in_face);
  Matrix dh(1,nodes_in_face);
  Matrix Jacobian(1,2);
  Matrix crd(1,2);

  Matrix fs(1,nodes_in_face);
  fs.Zero();
  
  Matrix NC = getFaceNodalCoords(face);

  for(short GP_c_r = 1 ; GP_c_r <= r_integration_order ; GP_c_r++) {
    r = get_Gauss_p_c(r_integration_order, GP_c_r);
    rw = get_Gauss_p_w(r_integration_order, GP_c_r);
    
    dh = diff_interp_fun_face(r);
    Jacobian = dh*NC;
    
    x1 = Jacobian(0,0);
    y1 = Jacobian(0,1);
    len = sqrt(x1*x1 + y1*y1);
    
    h = interp_fun_face(r);
    crd = h*NC;
    
    if (varCrd == 1) { // pressure varies with x coordinates
      pCrd = crd(0,0);
    }
    else { // pressure varies with y coordinates
      pCrd = crd(0,1);
    }
    
    if (top == bot) {
      printf("WARNING zero pressure point is the same as maximum pressure point!\n");
      pMag = -1.0;
    }
    else {
      pMag = (top - pCrd)/(top - bot);
    }
    
    if (pMag < 0.0 || pMag > 1.0) {
      pMag = 0.0;
    }
    else {
      pMag *= mag;
    }
    printf("pMag is %g!\n", pMag);
    
    if (pMag != 0.0) {
      weight = rw*pMag*len;
      fs.addMatrix(1.0, h, weight);
    }
    
  }
  
  // get node index of face in element 
  ID face_nodes(nodes_in_face);
  localFaceMapping(face, face_nodes);
  
  for (int i = 0; i < nodes_in_face; i++) {
    Q(face_nodes(i)) = fs(0,i);
  }
  
  for (int i = 0; i < numDOF; i++) {
    printf("Q(%d) = %g;\n", i+1, Q(i));
  }
  
}


void 
AC3D8HexWithSensitivity::addUniformPressure(int face, double mag)
{
  if (face < 1 || face > 4) {
    opserr << "AC3D8HexWithSensitivity::addUniformPressure - invalid face number:" << face << "\n";
    return ;
  }
  
  // see if we can quick return
  if (mag == 0.0) 
    return ;
  
  double r, rw;
  double weight = 0.0;
  double x1, y1,len;
  int nodes_in_face = 3;
  
  Matrix h(1,nodes_in_face);
  Matrix dh(1,nodes_in_face);
  Matrix Jacobian(1,2);

  Matrix fs(1,nodes_in_face);
  fs.Zero();
  
  Matrix NC = getFaceNodalCoords(face);
  
  for(short GP_c_r = 1 ; GP_c_r <= r_integration_order ; GP_c_r++) {
    r = get_Gauss_p_c(r_integration_order, GP_c_r);
    rw = get_Gauss_p_w(r_integration_order, GP_c_r);
    
    dh = diff_interp_fun_face(r);
    Jacobian = dh*NC;
    
    x1 = Jacobian(0,0);
    y1 = Jacobian(0,1);
    len = sqrt(x1*x1 + y1*y1);
    
    h = interp_fun_face(r);
    
    weight = rw*mag*len;
    fs.addMatrix(1.0, h, weight);
  }
  
  // get node index of face in element 
  ID face_nodes(nodes_in_face);
  localFaceMapping(face, face_nodes);
  
  for (int i = 0; i < nodes_in_face; i++) {
    Q(face_nodes(i)) = fs(0,i);
  }
  
  for (int i = 0; i < numDOF; i++) {
    printf("Q(%d) = %g;\n", i+1, Q(i));
  }
}
*/

//
// /////////////////////////////////////////////////////////////////////////////
//


Vector 
AC3D8HexWithSensitivity::nodal_forces_from_displacement(const Vector &u)
{
  Vector nodalF(numDOF);
  
  double r, s, t;
  double rw, sw, tw;
  double weight = 0.0;
  double det_of_Jacobian = 0.0;

  short where = 0;
  
  Matrix sigma(1,3);
  Matrix NF(1,nodes_in_elem);
  
  Vector epsilon(3);
  Matrix sstrain(3,1);
  Matrix tmp_disp(nodes_in_elem, 1);
  
  int i;
  for(i = 0; i < nodes_in_elem; i++) {
    tmp_disp(i,0) = u(i);
  }
  
  this->computeDiff();
  
  for(short GP_c_r = 1; GP_c_r <= r_integration_order; GP_c_r++) {
    r = get_Gauss_p_c(r_integration_order, GP_c_r);
    rw = get_Gauss_p_w(r_integration_order, GP_c_r);
    for(short GP_c_s = 1; GP_c_s <= s_integration_order; GP_c_s++) {
      s = get_Gauss_p_c(s_integration_order, GP_c_s);
      sw = get_Gauss_p_w(s_integration_order, GP_c_s);
      for(short GP_c_t = 1; GP_c_t <= t_integration_order; GP_c_t++) {
        t = get_Gauss_p_c(t_integration_order, GP_c_t);
        tw = get_Gauss_p_w(t_integration_order, GP_c_t);
        
        det_of_Jacobian = detJ[where];
        Matrix &dhGlobal = *L[where];
        
        sstrain.addMatrixProduct(0.0, dhGlobal, tmp_disp, 1.0);
        epsilon(0) = sstrain(0,0);
        epsilon(1) = sstrain(1,0);
        epsilon(2) = sstrain(2,0);
        
        theMaterial[where]->setTrialStrain(epsilon);
        
        weight = sw * rw * tw * det_of_Jacobian;
        
        // assemble stress tensor
        const Vector &stressvec = theMaterial[where]->getStress();
        sigma(0,0) = stressvec(0);
        sigma(0,1) = stressvec(1);
        sigma(0,2) = stressvec(2);
        
        NF.addMatrixProduct(1.0, sigma, dhGlobal, weight);
        
        where++;
      }
    }
  }
  
  for(i = 1; i <= nodes_in_elem; i++) {
    nodalF(i) = NF(0,i);
  }
  
  return nodalF;
}

//
// /////////////////////////////////////////////////////////////////////////////
//

double 
AC3D8HexWithSensitivity::get_Gauss_p_c(short order, short point_numb)
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
AC3D8HexWithSensitivity::get_Gauss_p_w(short order, short point_numb)
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
AC3D8HexWithSensitivity::setParameter(const char **argv, int argc, Parameter &param)
{


	int numberGauss=8;

    if (strstr(argv[0],"material") != 0) {
		int ok;
		for (int i=0; i<numberGauss; i++) {
			ok = theMaterial[i]->setParameter(&argv[1], argc-1, param);
			if (ok < 0 ) {
				opserr<<"AC3D8HexWithSensitivity::setParameter() can not setParameter for "<<i<<"th Gauss Point\n";
				return -1;
			}
		}
		return ok;
	}
	else{
		opserr<<"AC3D8HexWithSensitivity can not setParameter!"<<endln;
    // otherwise parameter is unknown for the Truss class
	return -1;
	}


}
 int
AC3D8HexWithSensitivity::updateParameter(int parameterID, Information &info)
{

  opserr<<"warnning: AC3D8HexWithSensitivity can not updateParameter!"<<endln;
  return -1;
}

 int
AC3D8HexWithSensitivity::activateParameter(int passedParameterID)
{

	int numberGauss=8;
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
		opserr << "AC3D8HexWithSensitivity::activateParameter() -- unknown parameter " << endln;
	}

	return 0;
}

 const Matrix &
AC3D8HexWithSensitivity::getMassSensitivity(int gradNumber)
{
	mass.Zero();
	return mass;
}

 const Vector &
AC3D8HexWithSensitivity::getResistingForceSensitivity(int gradNumber)
{
  double r  = 0.0;
  double rw = 0.0;
  double s  = 0.0;
  double sw = 0.0;
  double t  = 0.0;
  double tw = 0.0;
  double weight = 0.0;
  double det_of_Jacobian = 0.0;
  static const int nstress = 3 ;
  static Vector stress (nstress);
  short where = 0;
  
  Matrix sigma(1,3);
  Matrix NF(1,nodes_in_elem);
  
  this->computeDiff();
  NF.Zero();
  
  for(short GP_c_r = 1; GP_c_r <= r_integration_order; GP_c_r++) {
    r = get_Gauss_p_c(r_integration_order, GP_c_r);
    rw = get_Gauss_p_w(r_integration_order, GP_c_r);
    for( short GP_c_s = 1; GP_c_s <= s_integration_order; GP_c_s++ ) {
      s = get_Gauss_p_c(s_integration_order, GP_c_s);
      sw = get_Gauss_p_w(s_integration_order, GP_c_s);
      for(short GP_c_t = 1; GP_c_t <= t_integration_order; GP_c_t++) {
        t = get_Gauss_p_c(t_integration_order, GP_c_t);
        tw = get_Gauss_p_w(t_integration_order, GP_c_t);
        
        det_of_Jacobian = detJ[where];
        Matrix &dhGlobal = *L[where];
        
        weight = sw * rw * tw * det_of_Jacobian;
        
        // assemble stress tensor
        const Vector &stressvec = theMaterial[where]->getStressSensitivity(gradNumber,true);
        sigma(0,0) = stressvec(0);
        sigma(0,1) = stressvec(1);
        sigma(0,2) = stressvec(2);
        
        NF.addMatrixProduct(1.0, sigma, dhGlobal, weight);

	//	opserr<<"sigma"<<sigma<<endln;
     //   opserr<<"dhGlobal"<<dhGlobal<<endln;
      //  opserr<<"weight"<<weight<<endln;
      //  opserr<<"NF"<<NF<<endln;

        
        // nodal forces See Zienkievicz part 1 pp 108
        // nodal_forces = nodal_forces + dhGlobal("ib")*stress_at_GP("ab")*weight;
        
        where++;
      }
    }
  }
  
  int i, counter = 0;
  //converting nodalforce matrix to vector
  for (i = 0; i < nodes_in_elem; i++){
    P(i) = NF(0,i);
  }
        
 // opserr<<"P in getResistingForceSensitivity is  "<<P<<endln;
  // P.addVector(1.0, Q, -1.0);
  
  return P;


  //return NF;



 }


int
AC3D8HexWithSensitivity::commitSensitivity(int gradNumber, int numGrads)
{
  double r  = 0.0;
  double rw = 0.0;
  double s  = 0.0;
  double sw = 0.0;
  double t  = 0.0;
  double tw = 0.0;
  double weight = 0.0;
  double det_of_Jacobian = 0.0;
  static const int nstress = 3 ;
  int success;

  static Vector stress (nstress);
  short where = 0;
  
  Matrix sigma(1,3);
  Matrix NF(1,nodes_in_elem);
  
  this->computeDiff();
  NF.Zero();

  static Matrix ul(1,3);

  for(short GP_c_r = 1; GP_c_r <= r_integration_order; GP_c_r++) {
    r = get_Gauss_p_c(r_integration_order, GP_c_r);
    rw = get_Gauss_p_w(r_integration_order, GP_c_r);
    for( short GP_c_s = 1; GP_c_s <= s_integration_order; GP_c_s++ ) {
      s = get_Gauss_p_c(s_integration_order, GP_c_s);
      sw = get_Gauss_p_w(s_integration_order, GP_c_s);
      for(short GP_c_t = 1; GP_c_t <= t_integration_order; GP_c_t++) {
        t = get_Gauss_p_c(t_integration_order, GP_c_t);
        tw = get_Gauss_p_w(t_integration_order, GP_c_t);
        
        det_of_Jacobian = detJ[where];
        Matrix &dhGlobal = *L[where];
        
        weight = sw * rw * tw * det_of_Jacobian;

        where++;

      }
    }
  }
  
 
  int i;

  Vector epsilon(3); 
  Matrix sstrain(3,1);
  
  ul(0,0)= theNodes[0]->getDispSensitivity(1,gradNumber);
  ul(0,1)= theNodes[1]->getDispSensitivity(2,gradNumber);
  ul(0,2)= theNodes[2]->getDispSensitivity(3,gradNumber);
 // opserr<<"ul"<<ul<<endln;

  for(i = 0; i < numGP; i++) {
    const Matrix &dhGlobal = *L[i];
    
    sstrain.addMatrixProduct(0.0, dhGlobal, ul, 1.0);
    epsilon(0) = sstrain(0,0);
    epsilon(1) = sstrain(1,0);
    epsilon(2) = sstrain(2,0);
  
    success = theMaterial[i]->commitSensitivity(epsilon,gradNumber,numGrads ) ;


   } 
  return success;


  //return NF;



 }
