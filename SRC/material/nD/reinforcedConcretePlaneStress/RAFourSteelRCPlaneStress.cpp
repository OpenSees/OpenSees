// Written: JZhong
// Created: 2004.11
//
// Description: This file contains the class definition for
// RAFourSteelRCPlaneStress
// For Detailed explanation of the model, please refer to the book
// entitled "Unified Theory of Concrete Structures,"
// by Thomas T.C. Hsu and Y.L. Mo, John Wiley & Sons, April 2010.

#include "RAFourSteelRCPlaneStress.h"
#include <UniaxialMaterial.h>
#include <Vector.h>
#include <Matrix.h>
#include <math.h>
#include <float.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <stdlib.h>

#include <DummyStream.h>
#include <MaterialResponse.h>
#include <elementAPI.h>
#define OPS_Export 

static int numRAFourSteelRCPPlaneStressMaterials = 0;

OPS_Export void * OPS_ADD_RUNTIME_VPV(OPS_RAFourSteelRCPlaneStressMaterial)
{
  if (numRAFourSteelRCPPlaneStressMaterials == 0) {
    numRAFourSteelRCPPlaneStressMaterials++;
    opserr << "RAFourSteelRCPPlaneStress unaxial material - Written by J.Zhong, Thomas T.C. Hsu and Y.L. Mo - Copyright@2009\n";
  }

  // Pointer to a uniaxial material that will be returned
  NDMaterial *theMaterial = 0;

  int numRemainingArgs = OPS_GetNumRemainingInputArgs();

  if (numRemainingArgs < 20) {
       opserr << "Want: nDMaterial RAFourSteelRCPlaneStress matTag? rho? UniaxiaMatTag1? UniaxiaMatTag2? UniaxiaMatTag3? UniaxiaMatTag4? UniaxiaMatTag5? UniaxiaMatTag6? angle1? angle2? angle3? angle4? rou1? rou2? rou3? rou4? fpc? fy? E0? epsc0?\n";
    return 0;	
  }

  int tag;
  double rho;
  int    iData[6];
  double dData[12];
  int numData = 0;

  numData = 1;
  if (OPS_GetInt(&numData, &tag) != 0) {
    opserr << "WARNING invalid uniaxialMaterial RAFourSteelRCPPlaneStress tag" << endln;
    return 0;
  }

  numData = 1;
  if (OPS_GetDouble(&numData, &rho) != 0) {
    opserr << "Invalid Arg rho: uniaxialMaterial RAFourSteelRCPPlaneStress tag: " << tag << endln;
    return 0;	
  }

  numData = 6;
  if (OPS_GetInt(&numData, iData) != 0) {
    opserr << "WARNING invalid uniaxialMaterial RAFourSteelRCPPlaneStress tag: " << tag << endln;
    return 0;
  }

  numData = 12;
  if (OPS_GetDouble(&numData, dData) != 0) {
    opserr << "WARNING invalid data RAFourSteelRCPPlaneStress tag: " << tag << endln;
    return 0;
  }

  UniaxialMaterial *theUniaxialMaterial1 = OPS_GetUniaxialMaterial(iData[0]);
    
  if (theUniaxialMaterial1 == 0) {
    opserr << "WARNING material not found\n";
    opserr << "Material: " << iData[0];
    opserr << "\nRAFourSteelRCPPlaneStress tag: " << tag << endln;
    return 0;
  }
  
  UniaxialMaterial *theUniaxialMaterial2 = OPS_GetUniaxialMaterial(iData[1]);

  if (theUniaxialMaterial2 == 0) {
    opserr << "WARNING material not found\n";
    opserr << "Material: " << iData[1];
    opserr << "\nRAFourSteelRCPPlaneStress tag: " << tag << endln;
    return 0;
  }
  
  UniaxialMaterial *theUniaxialMaterial3 = OPS_GetUniaxialMaterial(iData[2]);
  if (theUniaxialMaterial3 == 0) {
    opserr << "WARNING material not found\n";
    opserr << "Material: " << iData[2];
    opserr << "\nRAFourSteelRCPPlaneStress tag: " << tag << endln;
    return 0;
  }

  UniaxialMaterial *theUniaxialMaterial4 = OPS_GetUniaxialMaterial(iData[3]);  
  if (theUniaxialMaterial4 == 0) {
    opserr << "WARNING material not found\n";
    opserr << "Material: " << iData[3];
    opserr << "\nRAFourSteelRCPPlaneStress tag: " << tag << endln;
    return 0;
  }

  UniaxialMaterial *theUniaxialMaterial5 = OPS_GetUniaxialMaterial(iData[4]);  
  if (theUniaxialMaterial5 == 0) {
    opserr << "WARNING material not found\n";
    opserr << "Material: " << iData[4];
    opserr << "\nRAFourSteelRCPPlaneStress tag: " << tag << endln;
    return 0;
  }


  UniaxialMaterial *theUniaxialMaterial6 = OPS_GetUniaxialMaterial(iData[5]);  
  if (theUniaxialMaterial6 == 0) {
    opserr << "WARNING material not found\n";
    opserr << "Material: " << iData[5];
    opserr << "\nRAFourSteelRCPPlaneStress tag: " << tag << endln;
    return 0;
  }

  //now create the RAFourSteelRCPPlaneStress
  theMaterial = new RAFourSteelRCPlaneStress (tag, 
					      rho,
					      theUniaxialMaterial1, 
					      theUniaxialMaterial2, 
					      theUniaxialMaterial3, 
					      theUniaxialMaterial4,
					      theUniaxialMaterial5,
					      theUniaxialMaterial6,
					      dData[0],
					      dData[1],
					      dData[2],
					      dData[3],
					      dData[4],
					      dData[5],
					      dData[6],
					      dData[7],
					      dData[8],
					      dData[9],
					      dData[10],
					      dData[11]);
  
  if (theMaterial == 0) {
    opserr << "WARNING ran out of memory creating material\n";
    opserr << "RAFourSteelRCPPlaneStress tag: " << tag << endln;
    return 0;
  }

  return theMaterial;
}
 
RAFourSteelRCPlaneStress::RAFourSteelRCPlaneStress (int      tag, 
						    double   RHO,
						    UniaxialMaterial *s1,
						    UniaxialMaterial *s2,
						    UniaxialMaterial *s3,
						    UniaxialMaterial *s4,
						    UniaxialMaterial *c1,
						    UniaxialMaterial *c2,
						    double   ANGLE1,
						    double   ANGLE2,
						    double   ANGLE3,
						    double   ANGLE4,
						    double   ROU1,
						    double   ROU2,
						    double   ROU3,
						    double   ROU4,
						    double   FPC,
						    double   FY,
						    double   E,
						    double   EPSC0) :
  NDMaterial(tag, ND_TAG_RAFourSteelRCPlaneStress), rho(RHO), 
  angle1(ANGLE1), angle2(ANGLE2), angle3(ANGLE3), angle4(ANGLE4), 
  rou1(ROU1), rou2(ROU2), rou3(ROU3), rou4(ROU4),
  fpc(FPC), fy(FY), E0(E), epsc0(EPSC0), strain_vec(3), stress_vec(3),tangent_matrix(3,3)
{
  steelStatus = 0;
  dirStatus = 0;
  G12 = 0;
  citaStrain = 10;
  citaStress = 10;
  
  TOneReverseStatus = 0;         
  TOneNowMaxComStrain = 0.0;
  TOneLastMaxComStrain = 0.0;
  
  TTwoReverseStatus = 0;         
  TTwoNowMaxComStrain = 0.0;
  TTwoLastMaxComStrain = 0.0;
  
  COneReverseStatus = 0;         
  COneNowMaxComStrain = 0.0;
  COneLastMaxComStrain = 0.0;
  
  CTwoReverseStatus = 0;         
  CTwoNowMaxComStrain = 0.0;
  CTwoLastMaxComStrain = 0.0;
  
  lastStress[0] = 0.0;  // add at 7.28/04
  lastStress[1] = 0.0;  // add at 7.28
  lastStress[2] = 0.0;  // add at 7.28
  
  if ( fpc < 0.0 ) { fpc = -fpc; } // set fpc > 0
  
  theMaterial = 0;
  
  // Allocate pointers to Steel and concrete
  theMaterial = new UniaxialMaterial *[6];
  
  if ( theMaterial == 0 ) {
    opserr << " RAFourSteelRCPlaneStress::RAFourSteelRCPlaneStress - failed to allocate material array\n";
    exit(-1);
  }
  
  // Get the copy for theSteel1
  theMaterial[0] = s1->getCopy();
  // Check allocation    
  if ( theMaterial[0] == 0 ) {
    opserr << " RAFourSteelRCPlaneStress::RAFourSteelRCPlaneStress - failed to get a copy for steel1\n";
    exit(-1);
  }
  

  // Get the copy for theSteel2
  theMaterial[1] = s2->getCopy();	
  // Check allocation    
  if ( theMaterial[1] == 0 ) {
    opserr << " RAFourSteelRCPlaneStress::RAFourSteelRCPlaneStress - failed to get a copy for steel2\n";
    exit(-1);
  }
  
  // Get the copy for theSteel3
  theMaterial[2] = s3->getCopy();	
  // Check allocation    
  if ( theMaterial[2] == 0 ) {
    opserr << " RAFourSteelRCPlaneStress::RAFourSteelRCPlaneStress - failed to get a copy for steel3\n";
    exit(-1);
  }
  
  // Get the copy for theSteel4
  theMaterial[3] = s4->getCopy();	
  // Check allocation    
  if ( theMaterial[3] == 0 ) {
    opserr << " RAFourSteelRCPlaneStress::RAFourSteelRCPlaneStress - failed to get a copy for steel4\n";
    exit(-1);
  }
  
  
  
  // Get the copy for theConcrete1
  theMaterial[4] = c1->getCopy();	
  // Check allocation    
  if ( theMaterial[4] == 0 ) {
    opserr << " RAFourSteelRCPlaneStress::RAFourSteelRCPlaneStress - failed to get a copy for concrete1\n";
    exit(-1);
  }
  
  // Get the copy for theConcrete2
  theMaterial[5] = c2->getCopy();	
  // Check allocation    
  if ( theMaterial[5] == 0 ) {
    opserr << " RAFourSteelRCPlaneStress::RAFourSteelRCPlaneStress - failed to get a copy for concrete2\n";
    exit(-1);
  }
  
  /* FMK */
  theResponses = new Response *[8];  
  
  if ( theResponses == 0) {
    opserr << " ReinforcedConcretePlaneStress::ReinforcedConcretePlaneStress - failed allocate responses  array\n";
    exit(-1);
  }
  
  OPS_Stream *theDummyStream = new DummyStream();
  
  const char **argv = new const char *[1];
  argv[0] = "getCommittedStrain";
  theResponses[0] = theMaterial[0]->setResponse(argv, 1, *theDummyStream);
  theResponses[1] = theMaterial[1]->setResponse(argv, 1, *theDummyStream);
  theResponses[2] = theMaterial[2]->setResponse(argv, 1, *theDummyStream);
  theResponses[3] = theMaterial[3]->setResponse(argv, 1, *theDummyStream);
  
  argv[0] = "setWallVar";
  theResponses[4] = theMaterial[4]->setResponse(argv, 1, *theDummyStream);
  theResponses[5] = theMaterial[5]->setResponse(argv, 1, *theDummyStream);
  
  argv[0] = "getPD";
  theResponses[6] = theMaterial[4]->setResponse(argv, 1, *theDummyStream);
  theResponses[7] = theMaterial[5]->setResponse(argv, 1, *theDummyStream);
  
  if ((theResponses[0] == 0) || (theResponses[1] == 0) ||
      (theResponses[2] == 0) || (theResponses[3] == 0) ||
      (theResponses[4] == 0) || (theResponses[5] == 0) ||
      (theResponses[6] == 0) || (theResponses[7] == 0)) {
    opserr << " ReinforcedConcretePlaneStress::ReinforcedConcretePlaneStress - failed to set appropriate concrete material\n";
    exit(-1);
  }
  
  // end FMK

  this->revertToStart();
}

RAFourSteelRCPlaneStress::RAFourSteelRCPlaneStress()
:NDMaterial(0, ND_TAG_RAFourSteelRCPlaneStress), strain_vec(3),
stress_vec(3),tangent_matrix(3,3)
{
  theMaterial = 0;
  theResponses = 0;
}


RAFourSteelRCPlaneStress::~RAFourSteelRCPlaneStress()
{
  // Delete the pointers
  if (theMaterial != 0) {
    for (int i=0; i<6; i++)
      {
	if (theMaterial[i])
	  delete theMaterial[i];
      }
    delete [] theMaterial;
  }
  if (theResponses != 0) {
    for (int j=0; j<8; j++)
      {
	if (theResponses[j] != 0)
	  delete theResponses[j];
      }
    delete [] theResponses;
  }    
}
								  

double RAFourSteelRCPlaneStress::getRho(void)
{
	return rho;
}


int RAFourSteelRCPlaneStress::setTrialStrain(const Vector &v)
{

	// Set values for strain_vec
	strain_vec = v;


	// Set initial values for Tstress
	Tstress[0] = 0.0;
	Tstress[1] = 0.0;
	Tstress[2] = 0.0;


    TOneReverseStatus = COneReverseStatus;         
	TOneNowMaxComStrain = COneNowMaxComStrain;
	TOneLastMaxComStrain = COneLastMaxComStrain;

	TTwoReverseStatus = CTwoReverseStatus;         
	TTwoNowMaxComStrain = CTwoNowMaxComStrain;
	TTwoLastMaxComStrain = CTwoLastMaxComStrain;



	determineTrialStress();
	
	return 0;
}

///*
int RAFourSteelRCPlaneStress::setTrialStrain(const Vector &v, const Vector &r)
{
   	// Set values for strain_vec
	strain_vec(0) = v(0);
	strain_vec(1) = v(1);
	strain_vec(2) = v(2);

	// Set initial values for Tstress
	Tstress[0] = 0.0;
	Tstress[1] = 0.0;
	Tstress[2] = 0.0;


    TOneReverseStatus = COneReverseStatus;         
	TOneNowMaxComStrain = COneNowMaxComStrain;
	TOneLastMaxComStrain = COneLastMaxComStrain;

	TTwoReverseStatus = CTwoReverseStatus;         
	TTwoNowMaxComStrain = CTwoNowMaxComStrain;
	TTwoLastMaxComStrain = CTwoLastMaxComStrain;


	determineTrialStress();
	
	return 0;
}


int RAFourSteelRCPlaneStress::setTrialStrainIncr(const Vector &v)
{
	// Set values for strain_vec
	strain_vec(0) = v(0);
	strain_vec(1) = v(1);
	strain_vec(2) = v(2);

	// Set initial values for Tstress
	Tstress[0] = 0.0;
	Tstress[1] = 0.0;
	Tstress[2] = 0.0;

    TOneReverseStatus = COneReverseStatus;         
	TOneNowMaxComStrain = COneNowMaxComStrain;
	TOneLastMaxComStrain = COneLastMaxComStrain;

	TTwoReverseStatus = CTwoReverseStatus;         
	TTwoNowMaxComStrain = CTwoNowMaxComStrain;
	TTwoLastMaxComStrain = CTwoLastMaxComStrain;

	determineTrialStress();
		
	return 0;
}


int RAFourSteelRCPlaneStress::setTrialStrainIncr(const Vector &v, const Vector &r)
{
	// Set values for strain_vec
	strain_vec(0) = v(0);
	strain_vec(1) = v(1);
	strain_vec(2) = v(2);

	// Set initial values for Tstress
	Tstress[0] = 0.0;
	Tstress[1] = 0.0;
	Tstress[2] = 0.0;

    TOneReverseStatus = COneReverseStatus;         
	TOneNowMaxComStrain = COneNowMaxComStrain;
	TOneLastMaxComStrain = COneLastMaxComStrain;

	TTwoReverseStatus = CTwoReverseStatus;         
	TTwoNowMaxComStrain = CTwoNowMaxComStrain;
	TTwoLastMaxComStrain = CTwoLastMaxComStrain;

	determineTrialStress();
		
	return 0;
}
//*/

const Matrix& RAFourSteelRCPlaneStress::getTangent(void)
{
	return tangent_matrix;
}

const Vector& RAFourSteelRCPlaneStress::getStress(void)
{

	return stress_vec;
}


const Vector& RAFourSteelRCPlaneStress :: getStrain()
{
	return strain_vec;
}
    

const Vector& RAFourSteelRCPlaneStress::getCommittedStress(void)
{
    return stress_vec;
}


const Vector& RAFourSteelRCPlaneStress::getCommittedStrain(void)
{
    return strain_vec;
}
    

int RAFourSteelRCPlaneStress::commitState(void)
{
    for (int i=0; i<6; i++)
	{
		theMaterial[i]->commitState();
	}

    COneReverseStatus = TOneReverseStatus;         
	COneNowMaxComStrain = TOneNowMaxComStrain;
	COneLastMaxComStrain = TOneLastMaxComStrain;

	CTwoReverseStatus = TTwoReverseStatus;         
	CTwoNowMaxComStrain = TTwoNowMaxComStrain;
	CTwoLastMaxComStrain = TTwoLastMaxComStrain;

	lastStress[0] = stress_vec(0);
	lastStress[1] = stress_vec(1);
	lastStress[2] = stress_vec(2);

	return 0;
}


int RAFourSteelRCPlaneStress::revertToLastCommit(void)
{

    for (int i=0; i<6; i++)
	{
		theMaterial[i]->revertToLastCommit();
	}

    TOneReverseStatus = COneReverseStatus;         
	TOneNowMaxComStrain = COneNowMaxComStrain;
	TOneLastMaxComStrain = COneLastMaxComStrain;

	TTwoReverseStatus = CTwoReverseStatus;         
	TTwoNowMaxComStrain = CTwoNowMaxComStrain;
	TTwoLastMaxComStrain = CTwoLastMaxComStrain;

	return 0;
}


int RAFourSteelRCPlaneStress::revertToStart(void)
{
  int i;
  for (i=0; i<6; i++)
    {
      theMaterial[i]->revertToStart();
    }
  
  for (i=0; i<3; i++)
    {
      Tstress[i] = 0.0;
    }
  strain_vec.Zero();
  stress_vec.Zero();
  
    steelStatus = 0;
    dirStatus = 0;
    G12 = 0;
    
    TOneReverseStatus = 0;         
    TOneNowMaxComStrain = 0.0;
    TOneLastMaxComStrain = 0.0;
    
    TTwoReverseStatus = 0;         
    TTwoNowMaxComStrain = 0.0;
    TTwoLastMaxComStrain = 0.0;
	
    COneReverseStatus = 0;         
    COneNowMaxComStrain = 0.0;
    COneLastMaxComStrain = 0.0;
    
    CTwoReverseStatus = 0;         
    CTwoNowMaxComStrain = 0.0;
    CTwoLastMaxComStrain = 0.0;

    return 0;
}

NDMaterial* RAFourSteelRCPlaneStress::getCopy(void)
{
  /*
    RAFourSteelRCPlaneStress *clone;
    clone = new RAFourSteelRCPlaneStress();
    *clone = *this;
    return clone;
  //*/
  
  RAFourSteelRCPlaneStress* theCopy =
    new RAFourSteelRCPlaneStress( this->getTag(), 
				  rho, theMaterial[0], theMaterial[1], theMaterial[2], theMaterial[3], theMaterial[4], theMaterial[5],
				  angle1, angle2, angle3, angle4, rou1, rou2, rou3, rou4, fpc, fy, E0, epsc0 );
  //theCopy->strain_vec = strain_vec;
  return theCopy;
}


NDMaterial* RAFourSteelRCPlaneStress::getCopy(const char *type)
{
  /*
    RAFourSteelRCPlaneStress *clone;
    clone = new RAFourSteelRCPlaneStress();
    *clone = *this;
    return clone;
  //*/
  RAFourSteelRCPlaneStress* theModel =
    new RAFourSteelRCPlaneStress( this->getTag(), 
				  rho, theMaterial[0], theMaterial[1], theMaterial[2], theMaterial[3], theMaterial[4], theMaterial[5],
				  angle1, angle2, angle3, angle4, rou1, rou2, rou3, rou4, fpc, fy, E0, epsc0 );
  //theModel->strain_vec = strain_vec;
  //theModel->stress_vec = stress_vec;
  return theModel;
}


void RAFourSteelRCPlaneStress::Print(OPS_Stream &s, int flag )
{
  s << "\n\tRAFourSteelRCPlaneStress, material id: " << this->getTag() << endln;
  /*
    s << "\tRho: " << rho << endln;
    s << "\tangle1: " << angle1 << endln;
    s << "\tangle2: " << angle2 << endln;
    s << "\tangle3: " << angle1 << endln;
    s << "\tangle4: " << angle2 << endln;
    s << "\trou1: " << rou1 << endln;
    s << "\trou2: " << rou2 << endln;
    s << "\trou3: " << rou1 << endln;
    s << "\trou4: " << rou2 << endln;
    s << "\tfpc: " << fpc << endln;
    s << "\tfy: " << fy << endln;
    s << "\tE0: " << E0 << endln;
  //*/
  
  //s << " x = " << xxx << endln;
  //s << " k = " << kkk << endln;
  
  //s <<  "Principal Strain: citaStrain = "<< citaStrain/3.14159*180.0 << endln;
  //s <<  "Principal Stress: citaStress = "<< citaStress/3.14159*180.0 << endln;
  //s << " v12 = " << miu12 << " v21 = " << miu21 << endln;	
  //s << " steelStatus " << steelStatus << endln;
  //s << " Damage DOne = " << DDOne << endln;
  //s << " Damage DTwo = " << DDTwo << endln;
  //s << " dirStatus " << dirStatus << endln;	
  //s << " G12 = " << G12 << endln;
  
  //s << " tt1 = " << tt1 << endln;
  //s << " tt2 = " << tt2 << endln;
  
  //int i;
  //int i, j;
  
   // s << "\tThe tangent stiffness is:" << endln;
   //for ( i=0; i<3; i++)
   //{
   //		for ( j=0; j<3; j++)
	//	{
	//		s << "  " << tangent_matrix(i,j);
	//	}
	//	s<< endln;
	//}


    //s << "\tThe strain is:" << endln;
	//for ( i=0; i<3; i++)
	//{
	//	s << "  " << strain_vec(i);		
	//	s<< endln;
	//}

	//s << "\tThe stress is:" << endln;
	//for ( i=0; i<3; i++)
	//{
	//	s << "  " << stress_vec(i);		
	//	s<< endln;
	//}


	s << "\t call the material print() function : "<< endln;
	
	s << "\t the steel 1 information is : " << endln;
	theMaterial[0]->Print(s,flag);
	s << "\t the steel 2 information is : " << endln;
	theMaterial[1]->Print(s,flag);
	s << "\t the steel 3 information is : " << endln;
	theMaterial[2]->Print(s,flag);
	s << "\t the steel 4 information is : " << endln;
	theMaterial[3]->Print(s,flag);
	s << "\t the concrete 1 information is : " << endln;
	theMaterial[4]->Print(s,flag);
	s << "\t the concrete 2 information is : " << endln;
	theMaterial[5]->Print(s,flag);


	//s << "\tStrain and stress of the uniaxial materials:"<<endln;
	//for ( i=0; i<6; i++)
	//{
	//	s<< "Uniaxial Material "<<i+1<<" :"<<theMaterial[i]->getStrain()<<"   "<<theMaterial[i]->getStress()<< endln;
	//}

}
 

int RAFourSteelRCPlaneStress::sendSelf(int commitTag, Channel &theChannel)
{
	int res = 0;
    
	int dataTag = this->getDbTag();

	// Packs its data into a Vector and sends this to theChannel
	static Vector data(13);
	data(0) = this->getTag();
	data(1) = rho;
	data(2) = angle1;
	data(3) = angle2;
	data(4) = angle3;
	data(5) = angle4;
	data(6) = rou1;
	data(7) = rou2;
	data(8) = rou3;
	data(9) = rou4;
	data(10) = fpc;
	data(11) = fy;
	data(12) = E0;
  
	res += theChannel.sendVector(dataTag, commitTag, data);
    if (res < 0) {
      opserr << "WARNING RAFourSteelRCPlaneStress::sendSelf() - " << this->getTag() << " failed to send Vector\n";
      return res;
	}	      

	
	// Now sends the IDs of its materials
    int matDbTag;
 
    static ID idData(12);

	// NOTE: to ensure that the material has a database
    // tag if sending to a database channel.

    int i;
    for (i=0; i<6; i++)
      {
	idData(i) = theMaterial[i]->getClassTag();
	matDbTag = theMaterial[i]->getDbTag();
	if (matDbTag == 0) {
	  matDbTag = theChannel.getDbTag();
	  if (matDbTag != 0)
	    theMaterial[i]->setDbTag(matDbTag);
	}
	idData(i+6) = matDbTag;
      }
    
    res += theChannel.sendID(dataTag, commitTag, idData);
    if (res < 0) {
      opserr << "WARNING RAFourSteelRCPlaneStress::sendSelf() - " << this->getTag() << " failed to send ID\n";
      return res;
    }
    
    // Finally, quad asks its material objects to send themselves
    for (i = 0; i < 6; i++) {
      res += theMaterial[i]->sendSelf(commitTag, theChannel);
      if (res < 0) {
	opserr << "RAFourSteelRCPlaneStress::sendSelf() - " << this->getTag() << " failed to send its Material\n";
	return res;
      }
    }	
    
    return res;
}
 

int RAFourSteelRCPlaneStress::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  int res = 0;
  
  int dataTag = this->getDbTag();

  // Quad creates a Vector, receives the Vector and then sets the 
  // internal data with the data in the Vector
  static Vector data(13);
  res += theChannel.recvVector(dataTag, commitTag, data);
  if (res < 0) {
    opserr << "WARNING RAFourSteelRCPlaneStress::recvSelf() - failed to receive Vector\n";
    return res;
  }
  
  this->setTag((int)data(0));
  rho = data(1);
  angle1 = data(2);
  angle2 = data(3);
  angle3 = data(4);
  angle4 = data(5);
  rou1 = data(6);
  rou2 = data(7);
  rou3 = data(8);
  rou4 = data(9);
  fpc = data(10);
  fy = data(11);
  E0 = data(12);

  static ID idData(12);
  
  // now receives the tags of its materials
  res += theChannel.recvID(dataTag, commitTag, idData);
  if (res < 0) {
    opserr << "WARNING RAFourSteelRCPlaneStress::recvSelf() - " << this->getTag() << " failed to receive ID\n";
    return res;
  }


  if (theMaterial == 0) {
    // Allocate new materials
    theMaterial = new UniaxialMaterial *[6];
    if (theMaterial == 0) {
      opserr << "RAFourSteelRCPlaneStress::recvSelf() - Could not allocate UniaxialMaterial* array\n";
      return -1;
    }
    for (int i = 0; i < 6; i++) {
      int matClassTag = idData(i);
      int matDbTag = idData(i+6);
      // Allocate new material with the sent class tag
      theMaterial[i] = theBroker.getNewUniaxialMaterial(matClassTag);
      if (theMaterial[i] == 0) {
	     opserr << "RAFourSteelRCPlaneStress::recvSelf() - Broker could not create NDMaterial of class type " << matClassTag << endln;
	     return -1;
      }
      // Now receive materials into the newly allocated space
      theMaterial[i]->setDbTag(matDbTag);
      res += theMaterial[i]->recvSelf(commitTag, theChannel, theBroker);
      if (res < 0) {
         opserr << "RAFourSteelRCPlaneStress::recvSelf() - material " << i << "failed to recv itself\n";
	     return res;
      }
    }
  }

  // materials exist , ensure materials of correct type and recvSelf on them
  else {
    for (int i = 0; i < 6; i++) {
      int matClassTag = idData(i);
      int matDbTag = idData(i+6);
      // Check that material is of the right type; if not,
      // delete it and create a new one of the right type
      if (theMaterial[i]->getClassTag() != matClassTag) {
	     delete theMaterial[i];
	     theMaterial[i] = theBroker.getNewUniaxialMaterial(matClassTag);
	     if (theMaterial[i] == 0) {
            opserr << "RAFourSteelRCPlaneStress::recvSelf() - material " << i << "failed to create\n";
	        return -1;
		 }
	  }
      // Receive the material
      theMaterial[i]->setDbTag(matDbTag);
      res += theMaterial[i]->recvSelf(commitTag, theChannel, theBroker);
      if (res < 0) {
           opserr << "RAFourSteelRCPlaneStress::recvSelf() - material " << i << "failed to recv itself\n";
	       return res;
	  }
    }
  }

  
  return res;

}


int RAFourSteelRCPlaneStress::determineTrialStress(void)
{   
	// Define i, j, k
	int i, j, k;

    // Definition of the variables and matrix
	//Transformation matrix
	double TOne[3][3];     // T(citaOne)
	double TMOne[3][3];    // T(minus citaOne)

	double TMSteelOne[3][3];     // T(minus citaSteelOne)
	double TSteelOne_One[3][3];  // T(citaSteelOne minus citaOne)

	double TMSteelTwo[3][3];     // T(minus citaSteelTwo)
	double TSteelTwo_One[3][3];  // T(citaSteelTwo minus citaOne)

	double TMSteelThr[3][3];     // T(minus citaSteelThree)
	double TSteelThr_One[3][3];  // T(citaSteelThree minus citaOne)

	double TMSteelFor[3][3];     // T(minus citaSteelFour)
	double TSteelFor_One[3][3];  // T(citaSteelFour minus citaOne)

	double V[3][3];        //Matrix considering Hsu/Zhu ratios

	//Strain and stress
	double Tstrain[3];     //epslonx,epslony,0.5*gammaxy	
	double tempStrain[3]; //temp values of strains
	double stressS1, stressS2,stressS3, stressS4; // stress of steel layers in 1, 2, 3 and 4 directions

	//stiffness of element
	double D[3][3];      //tangent stiffness matrix
	double DC[3][3];     //concrete part

	double DSOne[3][3];    //steel part in 1 direction
	double DSTwo[3][3];    //steel part in 2 direction
	double DSThr[3][3];    //steel part in 3 direction
	double DSFor[3][3];    //steel part in 4 direction

	double DC_bar[3][3];   //partial differentiation matrix of concrete part Eq.(49)	
	double tempD[3][3];  //temp matrix for data transfer

	double pi = 3.1415926535;
	double epsy = fy/E0;
	double fcr = 0.31*sqrt(fpc);
	double rou=0.0;   
	if ( rou < rou1 )	{	rou = rou1;	} 
	if ( rou < rou2 )	{	rou = rou2;	} 
	if ( rou < rou3 )	{	rou = rou3;	} 
	if ( rou < rou4 )	{	rou = rou4;	} 	

	if ( rou < 0.0025 ) rou = 0.0025;
	double B = pow((fcr/fy),1.5)/rou;
	double fn = fy * (0.91-2.0*B) / (0.98-0.25*B);

	double citaS1 = angle1;
	double citaS2 = angle2;
	double citaS3 = angle3;
	double citaS4 = angle4;


    //Set values for matrix TMSteelOne, TMSteelTwo, TMSteelThree and TMSteelFour
	TMSteelOne[0][0] = pow(cos(citaS1),2);
	TMSteelOne[0][1] = pow(sin(citaS1),2);
	TMSteelOne[0][2] = -2.0*cos(citaS1)*sin(citaS1);
	TMSteelOne[1][0] = pow(sin(citaS1),2);
	TMSteelOne[1][1] = pow(cos(citaS1),2);
	TMSteelOne[1][2] = 2.0*cos(citaS1)*sin(citaS1);
	TMSteelOne[2][0] = cos(citaS1)*sin(citaS1);
	TMSteelOne[2][1] = -cos(citaS1)*sin(citaS1);
	TMSteelOne[2][2] = pow(cos(citaS1),2)-pow(sin(citaS1),2);

	TMSteelTwo[0][0] = pow(cos(citaS2),2);
	TMSteelTwo[0][1] = pow(sin(citaS2),2);
	TMSteelTwo[0][2] = -2.0*cos(citaS2)*sin(citaS2);
	TMSteelTwo[1][0] = pow(sin(citaS2),2);
	TMSteelTwo[1][1] = pow(cos(citaS2),2);
	TMSteelTwo[1][2] = 2.0*cos(citaS2)*sin(citaS2);
	TMSteelTwo[2][0] = cos(citaS2)*sin(citaS2);
	TMSteelTwo[2][1] = -cos(citaS2)*sin(citaS2);
	TMSteelTwo[2][2] = pow(cos(citaS2),2)-pow(sin(citaS2),2);

	TMSteelThr[0][0] = pow(cos(citaS3),2);
	TMSteelThr[0][1] = pow(sin(citaS3),2);
	TMSteelThr[0][2] = -2.0*cos(citaS3)*sin(citaS3);
	TMSteelThr[1][0] = pow(sin(citaS3),2);
	TMSteelThr[1][1] = pow(cos(citaS3),2);
	TMSteelThr[1][2] = 2.0*cos(citaS3)*sin(citaS3);
	TMSteelThr[2][0] = cos(citaS3)*sin(citaS3);
	TMSteelThr[2][1] = -cos(citaS3)*sin(citaS3);
	TMSteelThr[2][2] = pow(cos(citaS3),2)-pow(sin(citaS3),2);

	TMSteelFor[0][0] = pow(cos(citaS4),2);
	TMSteelFor[0][1] = pow(sin(citaS4),2);
	TMSteelFor[0][2] = -2.0*cos(citaS4)*sin(citaS4);
	TMSteelFor[1][0] = pow(sin(citaS4),2);
	TMSteelFor[1][1] = pow(cos(citaS4),2);
	TMSteelFor[1][2] = 2.0*cos(citaS4)*sin(citaS4);
	TMSteelFor[2][0] = cos(citaS4)*sin(citaS4);
	TMSteelFor[2][1] = -cos(citaS4)*sin(citaS4);
	TMSteelFor[2][2] = pow(cos(citaS4),2)-pow(sin(citaS4),2);
    
    // Get strain values from strain of element
	Tstrain[0] = strain_vec(0);
	Tstrain[1] = strain_vec(1);
	Tstrain[2] = 0.5*strain_vec(2);

    // Get citaR based on Tstrain
	double citaR;
	double temp_citaR;
	double eps = 1e-12;


    
    if ( fabs(Tstrain[0]-Tstrain[1]) < 1e-7 )
	{
	   citaR = 0.25*pi;	
	}
	else // Tstrain[0] != Tstrain[1]
	{
	    temp_citaR = 0.5 * atan(fabs(2.0*1000000.0*Tstrain[2]/(1000000.0*Tstrain[0]-1000000.0*Tstrain[1]))); 
		if ( fabs(Tstrain[2]) < 1e-7 ) 
		{
			citaR = 0;
		} else if ( (Tstrain[0] > Tstrain[1]) && ( Tstrain[2] > 0) )
		{			
			citaR = temp_citaR;
		}
		else if ( (Tstrain[0] > Tstrain[1]) && ( Tstrain[2] < 0) )
		{
			citaR = pi - temp_citaR;
		}
		else if ( (Tstrain[0] < Tstrain[1]) && ( Tstrain[2] > 0) )
		{
			citaR = 0.5*pi - temp_citaR;
		}
		else if ( (Tstrain[0] < Tstrain[1]) && ( Tstrain[2] < 0) )
		{
			citaR = 0.5*pi + temp_citaR;
		}
		else
		{
			opserr << "RAFourSteelRCPlaneStress::determineTrialStress: Failure to calculate citaR\n";
			opserr << " Tstrain[0] = " << Tstrain[0] << endln;
			opserr << " Tstrain[1] = " << Tstrain[1] << endln;
			opserr << " Tstrain[2] = " << Tstrain[2] << endln;
		}
	}


    while (  (citaR - 0.5*pi) > 1e-8 ) {
		citaR = citaR-0.5*pi;
		dirStatus = 1;
	}

	citaStrain = citaR;
	

    // Set values for transformation matrix TOne[3][3],TMOne[3][3],
	// TSteelOne_One[3][3], TSteelTwo_One[3][3], TSteelThr_One[3][3], TSteelFor_One[3][3],
	TOne[0][0] = pow(cos(citaR),2);
	TOne[0][1] = pow(sin(citaR),2);
	TOne[0][2] = 2.0*cos(citaR)*sin(citaR);
	TOne[1][0] = pow(sin(citaR),2);
	TOne[1][1] = pow(cos(citaR),2);
	TOne[1][2] = -2.0*cos(citaR)*sin(citaR);
	TOne[2][0] = -cos(citaR)*sin(citaR);
	TOne[2][1] = cos(citaR)*sin(citaR);
	TOne[2][2] = pow(cos(citaR),2)-pow(sin(citaR),2);

    TMOne[0][0] = pow(cos(citaR),2);
	TMOne[0][1] = pow(sin(citaR),2);
	TMOne[0][2] = -2.0*cos(citaR)*sin(citaR);
	TMOne[1][0] = pow(sin(citaR),2);
	TMOne[1][1] = pow(cos(citaR),2);
	TMOne[1][2] = 2.0*cos(citaR)*sin(citaR);
	TMOne[2][0] = cos(citaR)*sin(citaR);
	TMOne[2][1] = -cos(citaR)*sin(citaR);
	TMOne[2][2] = pow(cos(citaR),2)-pow(sin(citaR),2);

    TSteelOne_One[0][0] = pow(cos(citaS1-citaR),2);
	TSteelOne_One[0][1] = pow(sin(citaS1-citaR),2);
	TSteelOne_One[0][2] = 2.0*cos(citaS1-citaR)*sin(citaS1-citaR);
	TSteelOne_One[1][0] = pow(sin(citaS1-citaR),2);
	TSteelOne_One[1][1] = pow(cos(citaS1-citaR),2);
	TSteelOne_One[1][2] = -2.0*cos(citaS1-citaR)*sin(citaS1-citaR);
	TSteelOne_One[2][0] = -cos(citaS1-citaR)*sin(citaS1-citaR);
	TSteelOne_One[2][1] = cos(citaS1-citaR)*sin(citaS1-citaR);
	TSteelOne_One[2][2] = pow(cos(citaS1-citaR),2)-pow(sin(citaS1-citaR),2);


    TSteelTwo_One[0][0] = pow(cos(citaS2-citaR),2);
	TSteelTwo_One[0][1] = pow(sin(citaS2-citaR),2);
	TSteelTwo_One[0][2] = 2.0*cos(citaS2-citaR)*sin(citaS2-citaR);
	TSteelTwo_One[1][0] = pow(sin(citaS2-citaR),2);
	TSteelTwo_One[1][1] = pow(cos(citaS2-citaR),2);
	TSteelTwo_One[1][2] = -2.0*cos(citaS2-citaR)*sin(citaS2-citaR);
	TSteelTwo_One[2][0] = -cos(citaS2-citaR)*sin(citaS2-citaR);
	TSteelTwo_One[2][1] = cos(citaS2-citaR)*sin(citaS2-citaR);
	TSteelTwo_One[2][2] = pow(cos(citaS2-citaR),2)-pow(sin(citaS2-citaR),2);

	TSteelThr_One[0][0] = pow(cos(citaS3-citaR),2);
	TSteelThr_One[0][1] = pow(sin(citaS3-citaR),2);
	TSteelThr_One[0][2] = 2.0*cos(citaS3-citaR)*sin(citaS3-citaR);
	TSteelThr_One[1][0] = pow(sin(citaS3-citaR),2);
	TSteelThr_One[1][1] = pow(cos(citaS3-citaR),2);
	TSteelThr_One[1][2] = -2.0*cos(citaS3-citaR)*sin(citaS3-citaR);
	TSteelThr_One[2][0] = -cos(citaS3-citaR)*sin(citaS3-citaR);
	TSteelThr_One[2][1] = cos(citaS3-citaR)*sin(citaS3-citaR);
	TSteelThr_One[2][2] = pow(cos(citaS3-citaR),2)-pow(sin(citaS3-citaR),2);

	TSteelFor_One[0][0] = pow(cos(citaS4-citaR),2);
	TSteelFor_One[0][1] = pow(sin(citaS4-citaR),2);
	TSteelFor_One[0][2] = 2.0*cos(citaS4-citaR)*sin(citaS4-citaR);
	TSteelFor_One[1][0] = pow(sin(citaS4-citaR),2);
	TSteelFor_One[1][1] = pow(cos(citaS4-citaR),2);
	TSteelFor_One[1][2] = -2.0*cos(citaS4-citaR)*sin(citaS4-citaR);
	TSteelFor_One[2][0] = -cos(citaS4-citaR)*sin(citaS4-citaR);
	TSteelFor_One[2][1] = cos(citaS4-citaR)*sin(citaS4-citaR);
	TSteelFor_One[2][2] = pow(cos(citaS4-citaR),2)-pow(sin(citaS4-citaR),2);


	//calculate tempStrain: epslon1,epslon2, 0.5*gamma12
    for (i=0;i<3;i++)
	{
		tempStrain[i] = 0.0;
		for (k=0;k<3;k++)
			tempStrain[i] += TOne[i][k]*Tstrain[k];	
	}

    double epslonOne, epslonTwo, halfGammaOneTwo;
	epslonOne = tempStrain[0];
	epslonTwo = tempStrain[1];
	halfGammaOneTwo = tempStrain[2];

    double epsn = epsy*(0.91-2.0*B)/(0.98-0.25*B); 

    // FMK
    //double CS1 = theMaterial[0]->getCommittedStrain();
    //double CS2 = theMaterial[1]->getCommittedStrain();
    //double CS3 = theMaterial[2]->getCommittedStrain();
    //double CS4 = theMaterial[3]->getCommittedStrain();

    theResponses[0]->getResponse();
    theResponses[1]->getResponse();
    theResponses[2]->getResponse();
    theResponses[3]->getResponse();
    Information &theInfo1 = theResponses[0]->getInformation();
    Information &theInfo2 = theResponses[1]->getInformation();
    Information &theInfo3 = theResponses[2]->getInformation();
    Information &theInfo4 = theResponses[3]->getInformation();
    double CS1 = theInfo1.theDouble;
    double CS2 = theInfo2.theDouble;
    double CS3 = theInfo3.theDouble;
    double CS4 = theInfo4.theDouble;
    
    // END FMK 
	
    if (( CS1 > epsn ) || ( CS2 > epsn ) || ( CS3 > epsn ) || ( CS4 > epsn ) ) {
      steelStatus = 1;
    }



	
    //set v12 and v21 obtained from strain
	double strainS1, strainS2, strainS3, strainS4; //Biaxial strain of steel in 1, 2, 3 and 4
	double strainSF=0.0;           //larger one of strainSL, strainST

	strainS1= pow(cos(citaS1),2)*Tstrain[0] + pow(sin(citaS1),2)*Tstrain[1] + 2.0*sin(citaS1)*cos(citaS1)*Tstrain[2];
	strainS2= pow(cos(citaS2),2)*Tstrain[0] + pow(sin(citaS2),2)*Tstrain[1] + 2.0*sin(citaS2)*cos(citaS2)*Tstrain[2];
	strainS3= pow(cos(citaS3),2)*Tstrain[0] + pow(sin(citaS3),2)*Tstrain[1] + 2.0*sin(citaS3)*cos(citaS3)*Tstrain[2];
	strainS4= pow(cos(citaS4),2)*Tstrain[0] + pow(sin(citaS4),2)*Tstrain[1] + 2.0*sin(citaS4)*cos(citaS4)*Tstrain[2];
 
	
	if (strainSF < strainS1)	{	strainSF = strainS1; }
	if (strainSF < strainS2)	{	strainSF = strainS2; }
	if (strainSF < strainS3)	{	strainSF = strainS3; }
	if (strainSF < strainS4)	{	strainSF = strainS4; }

	double v12 = 0.2;
	double v21 = 0.2; //initial values for Hsu/Zhu ratios
    
   
	if ( (epslonOne > 0.0) && (epslonTwo > 0.0) )
	{
		v12 = 0.0;
		v21 = 0.0;
	}
	else if ( (epslonOne > 0.0) && (epslonTwo <= 0.0) )
	{
		v21 = 0.0;
		if (strainSF > 0.002 )
		{
			//v12 = 0.5;
			v12 = 1.0;
			//v12 = 1.9;
		}
		else if (strainSF < 0)
		{
			v12 = 0.2;
		}
		else
		{
			//v12 = 0.2 + 150*strainSF;
			v12 = 0.2 + 400.0*strainSF;
			//v12 = 0.2 + 850.0*strainSF;
		}
		
		if (steelStatus == 1)
			v12 = 1.0;
			//v12 = 1.9;
			//v12 = 0.5;

	}
	else if ( (epslonOne <= 0.0) && (epslonTwo > 0.0) )
	{
		v12 = 0.0;
        if (strainSF > 0.002)
		{
			//v21 = 0.5;
			v21 = 1.0;
			//v21 = 1.9;
		}
		else if (strainSF < 0)
		{
			v21 = 0.2;
		}
		else
		{
			//v21 = 0.2 + 150*strainSF;
			v21 = 0.2 + 400.0*strainSF;
			//v21 = 0.2 + 850.0*strainSF;
		}
		
		if (steelStatus == 1)
			v21 = 1.0;
			//v21 = 1.9;
			//v21 = 0.5;
	}
    else if ( (epslonOne <= 0.0) && (epslonTwo <= 0.0) )
	{
		if (strainSF > 0.002)
		{
			v21 = 0.95; //0.95
			v12 = 0.95; //0.95
			//v21 = 1.9;
			//v12 = 1.9;
			//v21 = 0.5;
			//v12 = 0.5;
		}
		else if (strainSF < 0)
		{
			v21 = 0.2;
			v12 = 0.2;
		}
		else
		{
			//v21 = 0.2 + 850.0*strainSF;
			//v12 = 0.2 + 850.0*strainSF;
			v21 = 0.2 + 375.0*strainSF; //375
			v12 = 0.2 + 375.0*strainSF; //375
			//v21 = 0.2 + 150.0*strainSF;
			//v12 = 0.2 + 150.0*strainSF;
		}
		if (steelStatus == 1)
		{
			//v21 = 1.9;
			//v12 = 1.9;
			v21 = 0.95; //0.95
			v12 = 0.95; //0.95
			//v21 = 0.5;
			//v12 = 0.5;
		}
	}
	else
	{
    	opserr << "RAFourSteelRCPlaneStress::determineTrialStress: failure to get Hsu/Zhu ratio!\n";
	}

  
    miu12 = v12;
	miu21 = v21;


	//set values of matrix V[3][3]
	if ( v12*v21==1.0 )
	{
		opserr << "RAFourSteelRCPlaneStress::determineTrialStress: failure to get matrix [V]!\n";
		opserr << "v12= "<<v12 << endln;
		opserr << "v21= "<<v21 << endln;
		V[0][0] = 1.0;
		V[0][1] = 0.0;
		V[0][2] = 0.0;
		V[1][0] = 0.0;
		V[1][1] = 1.0;
		V[1][2] = 0.0;
		V[2][0] = 0.0;
		V[2][1] = 0.0;
		V[2][2] = 1.0;

	}
	else
	{
		V[0][0] = 1.0/(1.0-v12*v21);
		V[0][1] = v12/(1.0-v12*v21);
		V[0][2] = 0.0;
		V[1][0] = v21/(1.0-v12*v21);
		V[1][1] = 1.0/(1.0-v12*v21);
		V[1][2] = 0.0;
		V[2][0] = 0.0;
		V[2][1] = 0.0;
		V[2][2] = 1.0;
	}

     

    //********** get [DC]**************
    //set temp[DC]=[V]*[TOne]
	for (i=0;i<3;i++)
		for (j=0;j<3;j++)
		{
			tempD[i][j] = 0.0;
			for (k=0;k<3;k++)
				tempD[i][j] += V[i][k]*TOne[k][j];			
		}
    for (i=0;i<3;i++)
		for (j=0;j<3;j++)
		{
			DC[i][j] = tempD[i][j];			
		}
    
    //calculate epslon1_bar, epslon2_bar, 0.5*gamma12
    tempStrain[0] = V[0][0]*tempStrain[0] + V[0][1]*tempStrain[1]; //epslon1_bar
	tempStrain[1] = V[1][0]*tempStrain[0] + V[1][1]*tempStrain[1]; //epslon2_bar


	//get stiffness of unixail strain of concrete in 12 direction
    double cigmaOneC; //stress of concrete in 12 direction
	double cigmaTwoC;
	double GOneTwoC; //shear modulus of concrete in 12 direction

	double ita = 1.0; // set ita to 1.0
	
    
    // get Damage factor: DOne, DTwo
	if ( tempStrain[0] < 0 )
	{
		TOneReverseStatus = 0;
		if ( TOneNowMaxComStrain > tempStrain[0] )
            TOneNowMaxComStrain = tempStrain[0];
	} 
	else // tempStrain[0] > 0
	{
	  if ( TOneReverseStatus == 0 ) // first reverse from compressive strain
	  {
	     TOneReverseStatus = 1;
		 TOneLastMaxComStrain = COneNowMaxComStrain;
		 TOneNowMaxComStrain = 0.0;
	  }
	}

	if ( tempStrain[1] < 0 )
	{
		TTwoReverseStatus = 0;
		if ( TTwoNowMaxComStrain > tempStrain[1] )
            TTwoNowMaxComStrain = tempStrain[1];
	} 
	else // tempStrain[1] > 0
	{
	  if ( TTwoReverseStatus == 0 ) // first reverse from compressive strain
	  {
	     TTwoReverseStatus = 1;
		 TTwoLastMaxComStrain = CTwoNowMaxComStrain;
		 TTwoNowMaxComStrain = 0.0;
	  }
	}

	tt1 = tempStrain[0];
    tt2 = tempStrain[1];

	double DOne = 1.0;
	double DTwo = 1.0;
	if ( tempStrain[0] < 0 ) 
	{
		DOne = 1 - fabs(0.4*TTwoLastMaxComStrain/epsc0);
		if ( DOne < 0.2 )	DOne = 0.2;
	}
	if ( tempStrain[1] < 0 ) 
	{
		DTwo = 1 - fabs(0.4*TOneLastMaxComStrain/epsc0);
		if ( DTwo < 0.2 )	DTwo = 0.2;
	}
  
    DOne = 1.0;
	DTwo = 1.0;
    
	
	DDOne = DOne;  // assign values for screen output
	DDTwo = DTwo;

    //for xx, kk, m, xx=n, kk=delta, keci
	double xx,kk;
	if ( ((fabs(lastStress[0])+fabs(lastStress[0])+fabs(lastStress[0]))==0.0) || (fabs(lastStress[0])==0.0) ) {
		xx = 2.0;
		kk = 1.0;
	} else {
        if ( ( lastStress[0]<0.0 ) && ( lastStress[1]<0.0 ) )
		{
			if ( lastStress[2] == 0.0 ) 
			{
				xx = 2.0;
				kk = 1.0;
				//kk = 0;
			} else {
				
				   
	
				double keci = sqrt( (fabs(lastStress[0])/fabs(lastStress[2])+1) * (fabs(lastStress[1])/fabs(lastStress[2])+1) );
				
				xx = 2/pow(keci,3.0);
                if ( xx < 0.6 ) xx = 0.6;

				///*
				double a, b;
                if ( keci < 1.5 ) {
					a = 0;
					b = 1.0;
				} else {
					a = ( 1.3 - 1.0 ) / (1.9 - 1.5);
					b = 1.0 - a*1.5;
				}
				
				kk = a*keci + b;
                //*/
				//kk = 0.105*keci+1;


			}
		} 
		else if ( ( lastStress[0]>0.0 ) && ( lastStress[1]>0.0 ) )
		{
			kk = 1.0;		
			xx = 2.0;
		}
		else // under tension and compression
		{
			kk = 1.0;
  	        
			double keciN;
			if ( lastStress[0]<0 ) {
			keciN= sqrt(fabs(lastStress[0])/fabs(lastStress[2])+1);
			} else {
            keciN= sqrt(fabs(lastStress[1])/fabs(lastStress[2])+1);
			}
			xx = 2/pow(keciN,3.0);
            if ( xx < 0.6 ) xx = 0.6;

		}
	}
    
	kk = 1.0;  // for no axial load cases
	xx = 2.0;  // for no axial load cases
   
	xxx = xx;
	kkk = kk;
 

	//for ConcreteZ01, ConcreteZ02
	// FMK
	//theMaterial[4]->setTrialStrain(xx, kk, DOne,ita,tempStrain[1],tempStrain[0]); 
	//theMaterial[5]->setTrialStrain(xx, kk, DTwo,ita,tempStrain[0],tempStrain[1]); 
	Information &theInfoC02 = theResponses[4]->getInformation();
	Information &theInfoC03 = theResponses[5]->getInformation();
    
	static Vector theData(5);
	theData(0) = xx;
	theData(1) = kk;
	theData(2) = DOne;
	theData(3) = ita;
	theData(4) = tempStrain[1];
	theInfoC02.setVector(theData);
	theResponses[4]->getResponse();
	theData(2) = DTwo;
	theData(4) = tempStrain[0];
	theInfoC03.setVector(theData);
	theResponses[5]->getResponse();

	theMaterial[4]->setTrialStrain(tempStrain[0]); 
	theMaterial[5]->setTrialStrain(tempStrain[1]); 
	
	// end FMK

	cigmaOneC = theMaterial[4]->getStress();
	cigmaTwoC = theMaterial[5]->getStress();

    // set GOneTwoC = 1.0;
	if (epslonOne == epslonTwo)
	{
		GOneTwoC = 20000.0;
	}
	else
	{
		// change at 03/19/04
		GOneTwoC = fabs((cigmaOneC-cigmaTwoC)/(epslonOne-epslonTwo));
		if (GOneTwoC > 20000.0) // >EGc
		GOneTwoC = 20000.0;
	}
    
	GOneTwoC = 0.01; // For RA Model

    G12 = GOneTwoC;


	DC_bar[0][0] = theMaterial[4]->getTangent();
	//FMK DC_bar[0][1] = theMaterial[4]->getPD();
	theResponses[6]->getResponse();
	Information &theInfoC1 = theResponses[6]->getInformation();
	DC_bar[0][1] = theInfoC1.theDouble;
	// end FMK

	//DC_bar[0][1] = 0.0;
	DC_bar[0][2] = 0.0;

	// FMK DC_bar[1][0] = theMaterial[5]->getPD();
	theResponses[7]->getResponse();
	Information &theInfoC2 = theResponses[7]->getInformation();
	DC_bar[1][0] = theInfoC2.theDouble;
	// end FMK

	//DC_bar[1][0] = 0.0;
	DC_bar[1][1] = theMaterial[5]->getTangent();
	DC_bar[1][2] = 0.0;

	DC_bar[2][0] = 0.0;
	DC_bar[2][1] = 0.0;	
	DC_bar[2][2] = GOneTwoC;

	//before [DC]=[v]*[TOne], now update [DC]=[Dc_bar]*[V]*[TOne]
	for (i=0;i<3;i++)
		for (j=0;j<3;j++)
		{
			tempD[i][j] = 0.0;
			for (k=0;k<3;k++)
				tempD[i][j] += DC_bar[i][k]*DC[k][j];
		}
     for (i=0;i<3;i++)
		for (j=0;j<3;j++)
		{
			DC[i][j] = tempD[i][j];			
		}

	//update [DC]=[TMOne]*[Dc_bar]*[V]*[TOne]
	for (i=0;i<3;i++)
		for (j=0;j<3;j++)
		{
			tempD[i][j] = 0.0;
			for (k=0;k<3;k++)
				tempD[i][j] += TMOne[i][k]*DC[k][j];
		}
     for (i=0;i<3;i++)
		for (j=0;j<3;j++)
		{
			DC[i][j] = tempD[i][j];			
		}



    //***************** get [DSOne] ******************
    //get [DSOne]=[V][TOne]
    for (i=0;i<3;i++)
		for (j=0;j<3;j++)
		{
			tempD[i][j] = 0.0;
			for (k=0;k<3;k++)
				tempD[i][j] += V[i][k]*TOne[k][j];
		}   
     for (i=0;i<3;i++)
		for (j=0;j<3;j++)
		{
			DSOne[i][j] = tempD[i][j];			
		}

	//get [DSOne]=[TSteelOneMOne][V][TOne]
     for (i=0;i<3;i++)
		for (j=0;j<3;j++)
		{
			tempD[i][j] = 0.0;
			for (k=0;k<3;k++)
				tempD[i][j] += TSteelOne_One[i][k]*DSOne[k][j];
		}   
     for (i=0;i<3;i++)
		for (j=0;j<3;j++)
		{
			DSOne[i][j] = tempD[i][j];			
		}
    
     double strainS1_b; //uniaxial strain of steel in 1 direction
	 strainS1_b = DSOne[0][0]*Tstrain[0] + DSOne[0][1]*Tstrain[1] + DSOne[0][2]*Tstrain[2];
     

     //get stiffness
	 double tangentS1;
	 theMaterial[0]->setTrialStrain( strainS1_b );
	 tangentS1 = theMaterial[0]->getTangent();
	 stressS1 = theMaterial[0]->getStress();

	 for (j=0;j<3;j++)
	 {
		 DSOne[0][j] = rou1*tangentS1*DSOne[0][j];
		 DSOne[1][j] = 0.0;
		 DSOne[2][j] = 0.0;
	 }

     //get [DSOne]=[TMSOne][DsOne][TSteelOne_One][V][TOne]
     for (i=0;i<3;i++)
		for (j=0;j<3;j++)
		{
			tempD[i][j] = 0.0;
			for (k=0;k<3;k++)
				tempD[i][j] += TMSteelOne[i][k]*DSOne[k][j];
		}   
     for (i=0;i<3;i++)
		for (j=0;j<3;j++)
		{
			DSOne[i][j] = tempD[i][j];			
		}

 
  
    //***************** get [DSTwo] ******************
    //get [DSTwo]=[V][TOne]
    for (i=0;i<3;i++)
		for (j=0;j<3;j++)
		{
			tempD[i][j] = 0.0;
			for (k=0;k<3;k++)
				tempD[i][j] += V[i][k]*TOne[k][j];
		}   
     for (i=0;i<3;i++)
		for (j=0;j<3;j++)
		{
			DSTwo[i][j] = tempD[i][j];			
		}

	//get [DSTwo]=[TSteelTwoMOne][V][TOne]
     for (i=0;i<3;i++)
		for (j=0;j<3;j++)
		{
			tempD[i][j] = 0.0;
			for (k=0;k<3;k++)
				tempD[i][j] += TSteelTwo_One[i][k]*DSTwo[k][j];
		}   
     for (i=0;i<3;i++)
		for (j=0;j<3;j++)
		{
			DSTwo[i][j] = tempD[i][j];			
		}
    
     double strainS2_b; //uniaxial strain of steel in 2 direction
	 strainS2_b = DSTwo[0][0]*Tstrain[0] + DSTwo[0][1]*Tstrain[1] + DSTwo[0][2]*Tstrain[2];
     

     //get stiffness
	 double tangentS2;
	 theMaterial[1]->setTrialStrain( strainS2_b );
	 tangentS2 = theMaterial[1]->getTangent();
	 stressS2 = theMaterial[1]->getStress();

	 for (j=0;j<3;j++)
	 {
		 DSTwo[0][j] = rou2*tangentS2*DSTwo[0][j];
		 DSTwo[1][j] = 0.0;
		 DSTwo[2][j] = 0.0;
	 }

     //get [DSTwo]=[TMSTwo][DsTwo][TSteelTwo_One][V][TOne]
     for (i=0;i<3;i++)
		for (j=0;j<3;j++)
		{
			tempD[i][j] = 0.0;
			for (k=0;k<3;k++)
				tempD[i][j] += TMSteelTwo[i][k]*DSTwo[k][j];
		}   
     for (i=0;i<3;i++)
		for (j=0;j<3;j++)
		{
			DSTwo[i][j] = tempD[i][j];			
		}




    //***************** get [DSThr] ******************
    //get [DSThr]=[V][TOne]
    for (i=0;i<3;i++)
		for (j=0;j<3;j++)
		{
			tempD[i][j] = 0.0;
			for (k=0;k<3;k++)
				tempD[i][j] += V[i][k]*TOne[k][j];
		}   
     for (i=0;i<3;i++)
		for (j=0;j<3;j++)
		{
			DSThr[i][j] = tempD[i][j];			
		}

	//get [DSThr]=[TSteelThrMOne][V][TOne]
     for (i=0;i<3;i++)
		for (j=0;j<3;j++)
		{
			tempD[i][j] = 0.0;
			for (k=0;k<3;k++)
				tempD[i][j] += TSteelThr_One[i][k]*DSThr[k][j];
		}   
     for (i=0;i<3;i++)
		for (j=0;j<3;j++)
		{
			DSThr[i][j] = tempD[i][j];			
		}
    
     double strainS3_b; //uniaxial strain of steel in 3 direction
	 strainS3_b = DSThr[0][0]*Tstrain[0] + DSThr[0][1]*Tstrain[1] + DSThr[0][2]*Tstrain[2];
     

     //get stiffness
	 double tangentS3;
	 theMaterial[2]->setTrialStrain( strainS3_b );
	 tangentS3 = theMaterial[2]->getTangent();
	 stressS3 = theMaterial[2]->getStress();

	 for (j=0;j<3;j++)
	 {
		 DSThr[0][j] = rou3*tangentS3*DSThr[0][j];
		 DSThr[1][j] = 0.0;
		 DSThr[2][j] = 0.0;
	 }

     //get [DSThr]=[TMSThr][DsThr][TSteelThr_One][V][TOne]
     for (i=0;i<3;i++)
		for (j=0;j<3;j++)
		{
			tempD[i][j] = 0.0;
			for (k=0;k<3;k++)
				tempD[i][j] += TMSteelThr[i][k]*DSThr[k][j];
		}   
     for (i=0;i<3;i++)
		for (j=0;j<3;j++)
		{
			DSThr[i][j] = tempD[i][j];			
		}





    //***************** get [DSFor] ******************
    //get [DSFor]=[V][TOne]
    for (i=0;i<3;i++)
		for (j=0;j<3;j++)
		{
			tempD[i][j] = 0.0;
			for (k=0;k<3;k++)
				tempD[i][j] += V[i][k]*TOne[k][j];
		}   
     for (i=0;i<3;i++)
		for (j=0;j<3;j++)
		{
			DSFor[i][j] = tempD[i][j];			
		}

	//get [DSFor]=[TSteelForMOne][V][TOne]
     for (i=0;i<3;i++)
		for (j=0;j<3;j++)
		{
			tempD[i][j] = 0.0;
			for (k=0;k<3;k++)
				tempD[i][j] += TSteelFor_One[i][k]*DSFor[k][j];
		}   
     for (i=0;i<3;i++)
		for (j=0;j<3;j++)
		{
			DSFor[i][j] = tempD[i][j];			
		}
    
     double strainS4_b; //uniaxial strain of steel in 4 direction
	 strainS4_b = DSFor[0][0]*Tstrain[0] + DSFor[0][1]*Tstrain[1] + DSFor[0][2]*Tstrain[2];
     

     //get stiffness
	 double tangentS4;
	 theMaterial[3]->setTrialStrain( strainS4_b );
	 tangentS4 = theMaterial[3]->getTangent();
	 stressS4 = theMaterial[3]->getStress();

	 for (j=0;j<3;j++)
	 {
		 DSFor[0][j] = rou4*tangentS4*DSFor[0][j];
		 DSFor[1][j] = 0.0;
		 DSFor[2][j] = 0.0;
	 }

     //get [DSFor]=[TMSFor][DsFor][TSteelFor_One][V][TOne]
     for (i=0;i<3;i++)
		for (j=0;j<3;j++)
		{
			tempD[i][j] = 0.0;
			for (k=0;k<3;k++)
				tempD[i][j] += TMSteelFor[i][k]*DSFor[k][j];
		}   
     for (i=0;i<3;i++)
		for (j=0;j<3;j++)
		{
			DSFor[i][j] = tempD[i][j];			
		}








     
    //****************** get tangent_matrix  ****************    
	// Get tangent_matrix
	for (i=0;i<3;i++)
	{
		for (j=0;j<3;j++)
		{
			D[i][j] = 0.0;
			D[i][j] = DC[i][j] + DSOne[i][j] + DSTwo[i][j] + DSThr[i][j] + DSFor[i][j];
			tangent_matrix(i,j) = D[i][j];
		}
	}
    tangent_matrix(0,2) = 0.5*D[0][2];
	tangent_matrix(1,2) = 0.5*D[1][2];
	tangent_matrix(2,2) = 0.5*D[2][2];
      

    //**************** get Tstress and stress_vec ****************
    Tstress[0] = pow(cos(citaR),2)*cigmaOneC + pow(sin(citaR),2)*cigmaTwoC
		         - 2*sin(citaR)*cos(citaR)*halfGammaOneTwo*GOneTwoC
			     + pow(cos(citaS1),2)*rou1*stressS1 + pow(cos(citaS2),2)*rou2*stressS2
				 + pow(cos(citaS3),2)*rou3*stressS3 + pow(cos(citaS4),2)*rou4*stressS4;

	Tstress[1] = pow(sin(citaR),2)*cigmaOneC + pow(cos(citaR),2)*cigmaTwoC
		         + 2*sin(citaR)*cos(citaR)*halfGammaOneTwo*GOneTwoC
			     + pow(sin(citaS1),2)*rou1*stressS1 + pow(sin(citaS2),2)*rou2*stressS2
				 + pow(sin(citaS3),2)*rou3*stressS3 + pow(sin(citaS4),2)*rou4*stressS4;

	Tstress[2] = cos(citaR)*sin(citaR)*cigmaOneC - cos(citaR)*sin(citaR)*cigmaTwoC
		         + (pow(cos(citaR),2)-pow(sin(citaR),2))*halfGammaOneTwo*GOneTwoC
			     + cos(citaS1)*sin(citaS1)*rou1*stressS1 + cos(citaS2)*sin(citaS2)*rou2*stressS2
				 + cos(citaS3)*sin(citaS3)*rou3*stressS3 + cos(citaS4)*sin(citaS4)*rou4*stressS4;

    stress_vec(0) = Tstress[0];
	stress_vec(1) = Tstress[1];
	stress_vec(2) = Tstress[2];

   


	// get principal stress direction

	double citaP;
	double temp_citaP;

    if ( fabs(Tstress[0]-Tstress[1]) < 1e-7 )
	{
	   citaP = 0.25*pi;	
	}
	else // Tstrain[0] != Tstrain[1]
	{
	    temp_citaP = 0.5 * atan(fabs(2.0*1000000.0*Tstress[2]/(1000000.0*Tstress[0]-1000000.0*Tstress[1]))); 
		if ( fabs(Tstress[2]) < 1e-7 ) 
		{
			citaP = 0;
		} else if ( (Tstress[0] > Tstress[1]) && ( Tstress[2] > 0) )
		{			
			citaP = temp_citaP;
		}
		else if ( (Tstress[0] > Tstress[1]) && ( Tstress[2] < 0) )
		{
			citaP = pi - temp_citaP;
		}
		else if ( (Tstress[0] < Tstress[1]) && ( Tstress[2] > 0) )
		{
			citaP = 0.5*pi - temp_citaP;
		}
		else if ( (Tstress[0] < Tstress[1]) && ( Tstress[2] < 0) )
		{
			citaP = 0.5*pi + temp_citaP;
		}
		else
		{
			opserr << "RAFourSteelRCPlaneStress::determineTrialStress: Failure to calculate principal stress direction\n";
			opserr << " Tstress[0] = " << Tstress[0] << endln;
			opserr << " Tstress[1] = " << Tstress[1] << endln;
			opserr << " Tstress[2] = " << Tstress[2] << endln;
		}
	}


    while (  (citaP - 0.5*pi) > 1e-8 ) {
		citaP = citaP-0.5*pi;		
	}

	citaStress = citaP; // assign for screen output




	
	return 0;
}
