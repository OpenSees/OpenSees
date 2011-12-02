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
                                                                        
// $Revision: 1.5 $
// $Date: 2003-02-14 23:01:26 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/nD/PlateFiberMaterial.cpp,v $

//
// Ed "C++" Love
//
// Generic Plate Fiber Material
//


#include <PlateFiberMaterial.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>

//static vector and matrices
Vector  PlateFiberMaterial::stress(5);
Matrix  PlateFiberMaterial::tangent(5,5);

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

  this->strain(0) = strainFromElement(0); //11
  this->strain(1) = strainFromElement(1); //22
  this->strain(2) = strainFromElement(2); //12
  this->strain(3) = strainFromElement(3); //23
  this->strain(4) = strainFromElement(4); //31

  double norm;
  static Vector outOfPlaneStress(1);
  static Vector strainIncrement(1);
  static Vector threeDstress(6);
  static Vector threeDstrain(6);
  static Matrix threeDtangent(6,6);
  static Vector threeDstressCopy(6); 

  static Matrix threeDtangentCopy(6,6);
  static Matrix dd22(1,1);

  int i, j;
  int ii, jj;

  int count = 0;
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
    threeDstress = theMaterial->getStress();

    //three dimensional tangent 
    threeDtangent = theMaterial->getTangent();

    //NDmaterial strain order          = 11, 22, 33, 12, 23, 31 
    //PlateFiberMaterial strain order =  11, 22, 12, 23, 31, 33 

    //swap matrix indices to sort out-of-plane components 
    for (i=0; i<6; i++) {

      ii = this->indexMap(i);

      threeDstressCopy(ii) = threeDstress(i);

      for (j=0; j<6; j++) {

	jj = this->indexMap(j);

	threeDtangentCopy(ii,jj) = threeDtangent(i,j);

      }//end for j

    }//end for i


    //partitioned stresses and tangent
    
    outOfPlaneStress(0) = threeDstress(2);

    dd22(0,0) = threeDtangentCopy(5,5);

    //set norm
    norm = outOfPlaneStress.Norm();

    //int Solve(const Vector &V, Vector &res) const;
    //int Solve(const Matrix &M, Matrix &res) const;
    //condensation 
    dd22.Solve(outOfPlaneStress, strainIncrement);

    //update out of plane strains
    this->Tstrain22 -= strainIncrement(0);

    count++;
  } while (norm > tolerance && count < 10);

  return 0;
}


//send back the strain
const Vector& 
PlateFiberMaterial::getStrain()
{
  return this->strain;
}


//send back the stress 
const Vector&  
PlateFiberMaterial::getStress()
{
  const Vector &threeDstress = theMaterial->getStress();
  static Vector threeDstressCopy(6);

  //swap matrix indices to sort out-of-plane components 
  int i, ii;
  for (i=0; i<6; i++) {
    
    ii = this->indexMap(i);
    
    threeDstressCopy(ii) = threeDstress(i);
    
  }//end for i

  for (i=0; i<5; i++) 
    this->stress(i)     = threeDstressCopy(i);

  return this->stress;
}

//send back the tangent 
const Matrix&  
PlateFiberMaterial::getTangent()
{
  static Matrix dd11(5,5);
  static Matrix dd12(5,1);
  static Matrix dd21(1,5);
  static Matrix dd22(1,1);
  static Matrix dd22invdd21(1,5);

  static Matrix threeDtangentCopy(6,6);
  const Matrix &threeDtangent = theMaterial->getTangent();

  //swap matrix indices to sort out-of-plane components 
  int i,j, ii, jj;
  for (i=0; i<6; i++) {

    ii = this->indexMap(i);

    for (j=0; j<6; j++) {

      jj = this->indexMap(j);

      threeDtangentCopy(ii,jj) = threeDtangent(i,j);

    }//end for j

  }//end for i
  
  dd22(0,0) = threeDtangentCopy(5,5);

  for (i=0; i<5; i++) {

    dd12(i,0) = threeDtangentCopy(i,5);
    dd21(0,i) = threeDtangentCopy(5,i);
    
    for (int j=0; j<5; j++) 
      dd11(i,j) = threeDtangentCopy(i,j);
    
  }//end for i
    
  //int Solve(const Vector &V, Vector &res) const;
  //int Solve(const Matrix &M, Matrix &res) const;
  //condensation 
  dd22.Solve(dd21, dd22invdd21);
  this->tangent   = dd11; 
  this->tangent  -= (dd12*dd22invdd21);

  return this->tangent;
}


const Matrix&  
PlateFiberMaterial::getInitialTangent()
{
  opserr << "PlateFiberMaterial::getInitialTangent() - not yet implemented\n";
  return this->getTangent();
}


int 
PlateFiberMaterial::indexMap(int i)
{
  int ii;

  switch (i+1) { //add 1 for standard vector indices

    case 1 :
      ii = 1; 
      break;
 
    case 2 :
      ii = 2;
      break;

    case 3 :
      ii = 6;
      break;

    case 4 :
      ii = 3;
      break;

    case 5 :
      ii = 4;
      break;

    case 6 :
      ii = 5;
      break;

    default :
      ii = 1;
      break;

  } //end switch

  ii--;

  return ii;
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
  res = theChannel.sendID(this->getDbTag(), commitTag, idData);
  if (res < 0) {
    opserr << "BeamFiberMaterial::sendSelf() - failed to send id data\n";
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
      opserr << "BeamFiberMaterial::recvSelf() - failed to get a material of type: " << matClassTag << endln;
      return -1;
    }
  }
  theMaterial->setDbTag(idData(2));

  // recv a vector containing strains and set the strains
  static Vector vecData(1);
  res = theChannel.recvVector(this->getDbTag(), commitTag, vecData);
  if (res < 0) {
    opserr << "BeamFiberMaterial::sendSelf() - failed to send vector data\n";
    return res;
  }

  Cstrain22 = vecData(0);
  Tstrain22 = Cstrain22;

  // now receive the assocaited materials data
  res = theMaterial->recvSelf(commitTag, theChannel, theBroker);
  if (res < 0) 
    opserr << "BeamFiberMaterial::sendSelf() - failed to send vector material\n";
  
  return res;
}
 


