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
                                                                        
// $Revision: 1.6 $
// $Date: 2007-05-03 23:03:26 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/nD/PlateFiberMaterial.cpp,v $

//
// Ed "C++" Love
//
// Generic Plate Fiber Material
//


#include <PlateFiberMaterial.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <elementAPI.h>

//static vector and matrices
Vector  PlateFiberMaterial::stress(5);
Matrix  PlateFiberMaterial::tangent(5,5);

//      0  1  2  3  4  5
// ND: 11 22 33 12 23 31
// PF: 11 22 12 23 31 33

void* OPS_PlateFiberMaterial()
{
    int numdata = OPS_GetNumRemainingInputArgs();
    if (numdata < 2) {
	opserr << "WARNING insufficient arguments\n";
	opserr << "Want: nDMaterial PlateFiber tag? matTag?" << endln;
	return 0;
    }

    int tag[2];
    numdata = 2;
    if (OPS_GetIntInput(&numdata,tag)<0) {
	opserr << "WARNING invalid tags\n";
	return 0;
    }

    NDMaterial *threeDMaterial = OPS_getNDMaterial(tag[1]);
    if (threeDMaterial == 0) {
	opserr << "WARNING nD material does not exist\n";
	opserr << "nD material: " << tag[1];
	opserr << "\nPlateFiber nDMaterial: " << tag[0] << endln;
	return 0;
    }
      
    NDMaterial* mat = new PlateFiberMaterial( tag[0], *threeDMaterial );

    if (mat == 0) {
	opserr << "WARNING: failed to create PlaneStrain material\n";
	return 0;
    }

    return mat;
}

//null constructor
PlateFiberMaterial::PlateFiberMaterial() : 
NDMaterial(0, ND_TAG_PlateFiberMaterial), 
strain(5) 
{ 

}


//full constructor
PlateFiberMaterial::PlateFiberMaterial(   
				   int tag, 
                                   NDMaterial &the3DMaterial) :
NDMaterial(tag, ND_TAG_PlateFiberMaterial),
strain(5)
{
  theMaterial = the3DMaterial.getCopy("ThreeDimensional");

  Tstrain22 = 0.0;
  Cstrain22 = 0.0;
}



//destructor
PlateFiberMaterial::~PlateFiberMaterial() 
{ 
  delete theMaterial;
} 



//make a clone of this material
NDMaterial*
PlateFiberMaterial::getCopy() 
{
  PlateFiberMaterial *clone;   //new instance of this class

  clone = new PlateFiberMaterial(this->getTag(), 
                                   *theMaterial); //make the copy

  clone->Tstrain22 = this->Tstrain22;
  clone->Cstrain22 = this->Cstrain22;

  return clone;
}


//make a clone of this material
NDMaterial* 
PlateFiberMaterial::getCopy(const char *type) 
{
  return this->getCopy();
}


//send back order of strain in vector form
int 
PlateFiberMaterial::getOrder() const
{
  return 5;
}


const char*
PlateFiberMaterial::getType() const 
{
  return "PlateFiber"; 
}



//swap history variables
int 
PlateFiberMaterial::commitState() 
{
  Cstrain22 = Tstrain22;

  return theMaterial->commitState();
}



//revert to last saved state
int 
PlateFiberMaterial::revertToLastCommit()
{
  Tstrain22 = Cstrain22;

  return theMaterial->revertToLastCommit();
}


//revert to start
int
PlateFiberMaterial::revertToStart()
{
  this->Tstrain22 = 0.0;
  this->Cstrain22 = 0.0;

  return theMaterial->revertToStart();
}


//mass per unit volume
double
PlateFiberMaterial::getRho()
{
  return theMaterial->getRho();
}


//receive the strain
int 
PlateFiberMaterial::setTrialStrain(const Vector &strainFromElement)
{
  static const double tolerance = 1.0e-08;

  strain(0) = strainFromElement(0); //11
  strain(1) = strainFromElement(1); //22
  strain(2) = strainFromElement(2); //12
  strain(3) = strainFromElement(3); //23
  strain(4) = strainFromElement(4); //31

  double norm;
  double condensedStress;
  double strainIncrement;
  static Vector threeDstrain(6);
  double dd22;

  int count = 0;
  const int maxCount = 20;
  double norm0;

  //newton loop to solve for out-of-plane strains
  do {

    //set three dimensional strain
    threeDstrain(0) = this->strain(0);
    threeDstrain(1) = this->strain(1);
    threeDstrain(2) = this->Tstrain22;
    threeDstrain(3) = this->strain(2); 
    threeDstrain(4) = this->strain(3);
    threeDstrain(5) = this->strain(4);

    if (theMaterial->setTrialStrain(threeDstrain) < 0) {
      opserr << "PlateFiberMaterial::setTrialStrain - material failed in setTrialStrain() with strain " << threeDstrain;
      return -1;
    }

    //three dimensional stress
    const Vector &threeDstress = theMaterial->getStress();

    //three dimensional tangent 
    const Matrix &threeDtangent = theMaterial->getTangent();

    //NDmaterial strain order          = 11, 22, 33, 12, 23, 31 
    //PlateFiberMaterial strain order =  11, 22, 12, 23, 31, 33 

    condensedStress = threeDstress(2);

    dd22 = threeDtangent(2,2);

    //set norm
    norm = fabs(condensedStress);
    if (count == 0)
      norm0 = norm;

    //condensation 
    strainIncrement = condensedStress/dd22;

    //update out of plane strains
    Tstrain22 -= strainIncrement;

  } while (count++ < maxCount && norm > tolerance);

  return 0;
}


//send back the strain
const Vector& 
PlateFiberMaterial::getStrain()
{
  return strain;
}


//send back the stress 
const Vector&  
PlateFiberMaterial::getStress()
{
  const Vector &threeDstress = theMaterial->getStress();

  stress(0) = threeDstress(0);
  stress(1) = threeDstress(1);
  stress(2) = threeDstress(3);
  stress(3) = threeDstress(4);
  stress(4) = threeDstress(5);

  return stress;
}

const Vector& 
PlateFiberMaterial::getStressSensitivity(int gradIndex,
                                         bool conditional)
{
  const Vector &threeDstress = theMaterial->getStressSensitivity(gradIndex, conditional);

  stress(0) = threeDstress(0);
  stress(1) = threeDstress(1);
  stress(2) = threeDstress(3);
  stress(3) = threeDstress(4);
  stress(4) = threeDstress(5);

  const Matrix &threeDtangent = theMaterial->getTangent();

  static Vector dd12(5);
  dd12(0) = threeDtangent(0,2);
  dd12(1) = threeDtangent(1,2);
  dd12(2) = threeDtangent(3,2);
  dd12(3) = threeDtangent(4,2);
  dd12(4) = threeDtangent(5,2);

  double dd22 = threeDtangent(2,2);

  double sigma2 = threeDstress(2);

  double dd22sigma2 = sigma2/dd22;

  stress.addVector(1.0, dd12, -dd22sigma2);

  return stress;
}

//send back the tangent 
const Matrix&  
PlateFiberMaterial::getTangent()
{
  const Matrix &threeDtangent = theMaterial->getTangent();

  static Matrix dd11(5,5);
  dd11(0,0) = threeDtangent(0,0);
  dd11(1,0) = threeDtangent(1,0);
  dd11(2,0) = threeDtangent(3,0);
  dd11(3,0) = threeDtangent(4,0);
  dd11(4,0) = threeDtangent(5,0);

  dd11(0,1) = threeDtangent(0,1);
  dd11(1,1) = threeDtangent(1,1);
  dd11(2,1) = threeDtangent(3,1);
  dd11(3,1) = threeDtangent(4,1);
  dd11(4,1) = threeDtangent(5,1);

  dd11(0,2) = threeDtangent(0,3);
  dd11(1,2) = threeDtangent(1,3);
  dd11(2,2) = threeDtangent(3,3);
  dd11(3,2) = threeDtangent(4,3);
  dd11(4,2) = threeDtangent(5,3);

  dd11(0,3) = threeDtangent(0,4);
  dd11(1,3) = threeDtangent(1,4);
  dd11(2,3) = threeDtangent(3,4);
  dd11(3,3) = threeDtangent(4,4);
  dd11(4,3) = threeDtangent(5,4);

  dd11(0,4) = threeDtangent(0,5);
  dd11(1,4) = threeDtangent(1,5);
  dd11(2,4) = threeDtangent(3,5);
  dd11(3,4) = threeDtangent(4,5);
  dd11(4,4) = threeDtangent(5,5);

  static Matrix dd12(5,1);
  dd12(0,0) = threeDtangent(0,2);
  dd12(1,0) = threeDtangent(1,2);
  dd12(2,0) = threeDtangent(3,2);
  dd12(3,0) = threeDtangent(4,2);
  dd12(4,0) = threeDtangent(5,2);

  static Matrix dd21(1,5);
  dd21(0,0) = threeDtangent(2,0);
  dd21(0,1) = threeDtangent(2,1);
  dd21(0,2) = threeDtangent(2,3);
  dd21(0,3) = threeDtangent(2,4);
  dd21(0,4) = threeDtangent(2,5);

  double dd22 = threeDtangent(2,2);

  //int Solve(const Vector &V, Vector &res) const;
  //int Solve(const Matrix &M, Matrix &res) const;
  //condensation 
  static Matrix dd22invdd21(1,5);
  //dd22.Solve(dd21, dd22invdd21);
  dd22invdd21.addMatrix(0.0, dd21, 1.0/dd22);

  //this->tangent   = dd11; 
  //this->tangent  -= (dd12*dd22invdd21);
  dd11.addMatrixProduct(1.0, dd12, dd22invdd21, -1.0);
  tangent = dd11;

  return tangent;
}


const Matrix&  
PlateFiberMaterial::getInitialTangent()
{
  const Matrix &threeDtangent = theMaterial->getInitialTangent();

  static Matrix dd11(5,5);
  dd11(0,0) = threeDtangent(0,0);
  dd11(1,0) = threeDtangent(1,0);
  dd11(2,0) = threeDtangent(3,0);
  dd11(3,0) = threeDtangent(4,0);
  dd11(4,0) = threeDtangent(5,0);

  dd11(0,1) = threeDtangent(0,1);
  dd11(1,1) = threeDtangent(1,1);
  dd11(2,1) = threeDtangent(3,1);
  dd11(3,1) = threeDtangent(4,1);
  dd11(4,1) = threeDtangent(5,1);

  dd11(0,2) = threeDtangent(0,3);
  dd11(1,2) = threeDtangent(1,3);
  dd11(2,2) = threeDtangent(3,3);
  dd11(3,2) = threeDtangent(4,3);
  dd11(4,2) = threeDtangent(5,3);

  dd11(0,3) = threeDtangent(0,4);
  dd11(1,3) = threeDtangent(1,4);
  dd11(2,3) = threeDtangent(3,4);
  dd11(3,3) = threeDtangent(4,4);
  dd11(4,3) = threeDtangent(5,4);

  dd11(0,4) = threeDtangent(0,5);
  dd11(1,4) = threeDtangent(1,5);
  dd11(2,4) = threeDtangent(3,5);
  dd11(3,4) = threeDtangent(4,5);
  dd11(4,4) = threeDtangent(5,5);

  static Matrix dd12(5,1);
  dd12(0,0) = threeDtangent(0,2);
  dd12(1,0) = threeDtangent(1,2);
  dd12(2,0) = threeDtangent(3,2);
  dd12(3,0) = threeDtangent(4,2);
  dd12(4,0) = threeDtangent(5,2);

  static Matrix dd21(1,5);
  dd21(0,0) = threeDtangent(2,0);
  dd21(0,1) = threeDtangent(2,1);
  dd21(0,2) = threeDtangent(2,3);
  dd21(0,3) = threeDtangent(2,4);
  dd21(0,4) = threeDtangent(2,5);

  double dd22 = threeDtangent(2,2);

  //int Solve(const Vector &V, Vector &res) const;
  //int Solve(const Matrix &M, Matrix &res) const;
  //condensation 
  static Matrix dd22invdd21(1,5);
  //dd22.Solve(dd21, dd22invdd21);
  dd22invdd21.addMatrix(0.0, dd21, 1.0/dd22);

  //this->tangent   = dd11; 
  //this->tangent  -= (dd12*dd22invdd21);
  dd11.addMatrixProduct(1.0, dd12, dd22invdd21, -1.0);
  tangent = dd11;

  return tangent;
}

//print out data
void  
PlateFiberMaterial::Print(OPS_Stream &s, int flag)
{
  s << "General Plate Fiber Material \n";
  s << " Tag: " << this->getTag() << "\n"; 
  s << "using the 3D material : \n";

  theMaterial->Print(s, flag);

  return;
}


int 
PlateFiberMaterial::sendSelf(int commitTag, Channel &theChannel) 
{
  int res = 0;

  // put tag and assocaited materials class and database tags into an id and send it
  static ID idData(3);
  idData(0) = this->getTag();
  idData(1) = theMaterial->getClassTag();
  int matDbTag = theMaterial->getDbTag();
  if (matDbTag == 0) {
    matDbTag = theChannel.getDbTag();
    theMaterial->setDbTag(matDbTag);
  }
  idData(2) = matDbTag;
  
  res = theChannel.sendID(this->getDbTag(), commitTag, idData);
  if (res < 0) {
    opserr << "PlateFiberMaterial::sendSelf() - failed to send id data\n";
    return res;
  }

  // put the strains in a vector and send it
  static Vector vecData(1);
  vecData(0) = Cstrain22;

  res = theChannel.sendVector(this->getDbTag(), commitTag, vecData);
  if (res < 0) {
    opserr << "PlateFiberMaterial::sendSelf() - failed to send vector data\n";
    return res;
  }

  // now send the materials data
  res = theMaterial->sendSelf(commitTag, theChannel);
  if (res < 0) 
    opserr << "PlateFiberMaterial::sendSelf() - failed to send vector material\n";

  return res;
}

int 
PlateFiberMaterial::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  int res = 0;

  // recv an id containg the tag and associated materials class and db tags
  static ID idData(3);
  res = theChannel.recvID(this->getDbTag(), commitTag, idData);
  if (res < 0) {
    opserr << "PlateFiberMaterial::sendSelf() - failed to send id data\n";
    return res;
  }

  this->setTag(idData(0));
  int matClassTag = idData(1);

  // if the associated material has not yet been created or is of the wrong type
  // create a new material for recvSelf later
  if (theMaterial == 0 || theMaterial->getClassTag() != matClassTag) {
    if (theMaterial != 0)
      delete theMaterial;
    theMaterial = theBroker.getNewNDMaterial(matClassTag);
    if (theMaterial == 0) {
      opserr << "PlateFiberMaterial::recvSelf() - failed to get a material of type: " << matClassTag << endln;
      return -1;
    }
  }
  theMaterial->setDbTag(idData(2));

  // recv a vector containing strains and set the strains
  static Vector vecData(1);
  res = theChannel.recvVector(this->getDbTag(), commitTag, vecData);
  if (res < 0) {
    opserr << "PlateFiberMaterial::sendSelf() - failed to send vector data\n";
    return res;
  }

  Cstrain22 = vecData(0);
  Tstrain22 = Cstrain22;

  // now receive the assocaited materials data
  res = theMaterial->recvSelf(commitTag, theChannel, theBroker);
  if (res < 0) 
    opserr << "PlateFiberMaterial::sendSelf() - failed to send vector material\n";
  
  return res;
}
 
int
PlateFiberMaterial::setParameter(const char **argv, int argc,
				 Parameter &param)
{
  return theMaterial->setParameter(argv, argc, param);
}
