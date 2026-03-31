// Written: Quan Gu and Zhijian Qiu
// Created: 2013.7
//
// Reference:JP Conte, MK. Jagannath, Seismic relibility analysis of concrete 
// gravity dams, A Report on Research, Rice University, 1995. 
//           EA de Souza Neto, D Peri´c, DRJ Owen, Computational methods for 
// plasticity, Theory and applications (see pages 357 to 366), 2008.
// 
// 3D J2 plasticity model with linear isotropic and kinematic hardening
//  
// -------------------

#include <math.h>   
#include <stdlib.h>
#include <PlaneStressSimplifiedJ2.h>
#include <Information.h>
#include <ID.h>
#include <MaterialResponse.h>
#include <Parameter.h>
#include <FEM_ObjectBroker.h>



Matrix PlaneStressSimplifiedJ2::tmpMatrix(3,3);
Vector PlaneStressSimplifiedJ2::tmpVector(3);

// --- element: eps(1,1),eps(2,2),eps(3,3),2*eps(1,2),2*eps(2,3),2*eps(1,3) ----
// --- material strain: eps(1,1),eps(2,2),eps(3,3),eps(1,2),eps(2,3),eps(1,3) , same sign ----

// be careful! Here we use  eps(1,1),eps(2,2),2*eps(1,2). i.e., the same as that of element. 

#include <SimplifiedJ2.h>
#include <elementAPI.h>

void *
OPS_PlaneStressSimplifiedJ2(void) {

  int tag;
  double K, G, sig0, H_kin, H_iso;

  int numArgs = OPS_GetNumRemainingInputArgs();
  if (numArgs < 6) {
    opserr << "ndMaterial PlaneStressSimplifiedJ2 incorrect num args: want tag G K sig0 H_kin H_iso\n";
    return 0;
  }

  int iData[1];
  double dData[5];

  int numData = 1;
  if (OPS_GetInt(&numData, iData) != 0) {
    opserr << "WARNING invalid integer values: nDMaterial PlaneStressSimplifiedJ2 \n";
    return 0;
  }  
  tag = iData[0]; 

  numData = 5;  
  if (OPS_GetDouble(&numData, dData) != 0) {
      opserr << "WARNING invalid double values: nDMaterial PlaneStressSimplifiedJ2 " << tag << endln;
    return 0;
  }  
  G = dData[0];
  K = dData[1];
  sig0 = dData[2];
  H_kin = dData[3];
  H_iso = dData[4];

  SimplifiedJ2 the3DMaterial(tag, 
			     3,
			     G,   
			     K,
			     sig0,
			     H_kin,
			     H_iso);
  
  
  NDMaterial *theMaterial = new PlaneStressSimplifiedJ2 (tag, 
							 2,
							 the3DMaterial);
  return theMaterial;
}



PlaneStressSimplifiedJ2::PlaneStressSimplifiedJ2 (int pTag, 
						   int nd, 
						   NDMaterial &passed3DMaterial)
  : NDMaterial(pTag,ND_TAG_PlaneStressSimplifiedJ2), the3DMaterial(0), stress(3),
    strain(3), Cstress(3), Cstrain(3),theTangent(3,3),
    savedStrain33(0.0), CsavedStrain33(0.0)
{
  this->ndm = 2;    
  the3DMaterial = passed3DMaterial.getCopy();  
  
  stress.Zero();
  strain.Zero();
  
  Cstress.Zero();
  Cstrain.Zero();
}

PlaneStressSimplifiedJ2::PlaneStressSimplifiedJ2()
  :NDMaterial(0,ND_TAG_PlaneStressSimplifiedJ2), the3DMaterial(0), stress(3),
   strain(3), Cstress(3), Cstrain(3),theTangent(3,3),
   savedStrain33(0.0), CsavedStrain33(0.0)
{

}

PlaneStressSimplifiedJ2::~PlaneStressSimplifiedJ2()
{
  if (the3DMaterial != 0)
    delete the3DMaterial;
};

     

	
int PlaneStressSimplifiedJ2::plastIntegrator(){

	int maxIter = 25;
	double tol = 1e-12;
	double e33 = CsavedStrain33;

	int debugFlag =0;
	static int counter =0;
	counter++;
//	opserr<<"counter:"<<counter<<endln;

	if (fabs(e33)>tol ) {
	 // opserr<<"testing planestress j2 part "<<endln;
	//  debugFlag =1;
	}

	static Vector strain3D(6);
	static Vector stress3D(6);
	static Matrix tangent3D(6,6);
	
	strain3D(0) = strain(0);
	strain3D(1) = strain(1);
	strain3D(2) = e33;
	strain3D(3) = strain(2);
	strain3D(4) = 0.0;
	strain3D(5) = 0.0;

	the3DMaterial->setTrialStrain(strain3D);
	stress3D = the3DMaterial->getStress();
	tangent3D = the3DMaterial->getTangent();


	int i =0;

// ------ debug ---------

	if (debugFlag ==1){
		opserr<<"iteration number:" <<i<<endln;	
		opserr<<"strain is:" <<strain3D<<endln;
		opserr<<"stress is:"<<stress3D<<endln;
		opserr<<"tangent is:"<< tangent3D<<endln;
	
	}

	double e33_old=e33+1.0;

	while (( fabs(e33-e33_old)>tol) &&( fabs(stress3D(2))>tol) &&(i<maxIter)) {

	    e33_old = e33;		
		e33 -= stress3D(2)/tangent3D(2,2);
		strain3D(2) = e33;
	    the3DMaterial->setTrialStrain(strain3D);
	    stress3D = the3DMaterial->getStress();
		tangent3D = the3DMaterial->getTangent();

	if (debugFlag ==1){
		opserr<<"iteration number:" <<i<<endln;	
		opserr<<"strain is:" <<strain3D<<endln;
		opserr<<"stress is:"<<stress3D<<endln;
		opserr<<"tangent is:"<< tangent3D<<endln;
	
	}
	//   opserr.precision(16);
	//	opserr<<"iteration number is" <<i;	
	//	opserr<<": strain_zz is:" <<strain3D(2)<< ", stress is:"<<stress3D(2)<<endln;




		i++;

	} 

	if (( fabs(e33-e33_old)>tol) &&(fabs(stress3D(2))>tol)) {
		opserr<<"Fatal: PlaneStressSimplifiedJ2::plastIntegrator() can not find e33!"<<endln;
		exit(-1);
	}

	// --------- update the stress and tangent -----
	savedStrain33 = e33;

//	opserr<<"Total iteration number:" <<i<<endln;	


   stress(0) = stress3D(0);
   stress(1) = stress3D(1);
   stress(2) = stress3D(3);
   

   double D22 = tangent3D(2,2);
   static Vector D12(3);
   static Vector D21(3);
   static Matrix D11(3,3);

 D11(0,0)=tangent3D(0,0);
 D11(0,1)=tangent3D(0,1);
 D11(0,2)=tangent3D(0,3);
 D11(1,0)=tangent3D(1,0);
 D11(1,1)=tangent3D(1,1);
 D11(1,2)=tangent3D(1,3);
 D11(2,0)=tangent3D(3,0);
 D11(2,1)=tangent3D(3,1);
 D11(2,2)=tangent3D(3,3);

D12(0) = tangent3D(0,2);
D12(1) = tangent3D(1,2);
D12(2) = tangent3D(3,2);

D21(0) = tangent3D(2,0);
D21(1) = tangent3D(2,1);
D21(2) = tangent3D(2,3);

for( int i=0; i<3; i++)
  for (int j=0; j<3; j++)
	  theTangent(i,j) = D11(i,j)-1.0/D22*D12(i)*D21(j);

	if (debugFlag ==1){
		opserr<<"Final 2D tangent is:"<< theTangent<<endln;
	
	}
 
	return 0;

};
 

int PlaneStressSimplifiedJ2::setTrialStrain (const Vector &pStrain){

 
    strain = pStrain;
 
  // ----- change to real strain instead of eng. strain

  // strain[2] /=2.0;     be careful!           
  
  this->plastIntegrator();



	return 0;

};   

int PlaneStressSimplifiedJ2::setTrialStrain(const Vector &v, const Vector &r){

	return this->setTrialStrain ( v);

};

int PlaneStressSimplifiedJ2::setTrialStrainIncr(const Vector &v){
	
	// ----- change to real strain instead of eng. strain
   // ---- since all strain in material is the true strain, not eng.strain. 

		strain[0] = Cstrain[0]+v[0];
		strain[1] = Cstrain[1]+v[1];
		strain[2] = Cstrain[2]+v[2];     //  no need to divide by 2.0;
	  
	 this->plastIntegrator();

	 return 0;

};

int PlaneStressSimplifiedJ2::setTrialStrainIncr(const Vector &v, const Vector &r){

 

	return this->setTrialStrainIncr(v);


};

     // Calculates current tangent stiffness.

const Matrix & PlaneStressSimplifiedJ2::getTangent (void){
		return theTangent;

};
const Matrix & PlaneStressSimplifiedJ2::getInitialTangent (void){

 

	return this->getTangent();

};
        
     
const Vector & PlaneStressSimplifiedJ2::getStress (void){

  return stress;

};

const Vector & PlaneStressSimplifiedJ2::getStrain (void){

	return strain; 
};

const Vector & PlaneStressSimplifiedJ2::getCommittedStress (void){ 

    return Cstress;
};

const Vector & PlaneStressSimplifiedJ2::getCommittedStrain (void){

    return Cstrain; 

};


int PlaneStressSimplifiedJ2::commitState (void){


	CsavedStrain33 = savedStrain33; 
	Cstress = stress;
	Cstrain = strain;
	//CcumPlastStrainDev = cumPlastStrainDev;

	return the3DMaterial->commitState();
};

int PlaneStressSimplifiedJ2::revertToLastCommit (void)
{
  savedStrain33 = CsavedStrain33;
  stress = Cstress;
  strain = Cstrain;

  return the3DMaterial->revertToLastCommit();
};

int PlaneStressSimplifiedJ2::revertToStart(void)
{
  CsavedStrain33 = 0.0;
  Cstress.Zero();
  Cstrain.Zero();

  return the3DMaterial->revertToStart();
}



NDMaterial * PlaneStressSimplifiedJ2::getCopy (void){
    PlaneStressSimplifiedJ2 * theJ2 = new PlaneStressSimplifiedJ2(this->getTag(),this->ndm, *the3DMaterial);
    return theJ2;
};

NDMaterial * PlaneStressSimplifiedJ2::getCopy (const char *type){
  if (strcmp(type,"PlaneStress") == 0) {
    PlaneStressSimplifiedJ2 * theJ2 = new PlaneStressSimplifiedJ2(this->getTag(),this->ndm, *the3DMaterial);
    return theJ2;
  } else {
    return 0;
  }
};
 


int PlaneStressSimplifiedJ2::sendSelf(int commitTag, Channel &theChannel)
{
  int res = 0;
  int dbTag = this->getDbTag();

  static ID idData(4);
  idData(0) = this->getTag();
  idData(1) = the3DMaterial->getClassTag();
  int matDbTag = the3DMaterial->getDbTag();
  if (matDbTag == 0) {
    matDbTag = theChannel.getDbTag();
    the3DMaterial->setDbTag(matDbTag);
  }
  idData(2) = matDbTag;
  idData(3) = ndm;

  res = theChannel.sendID(dbTag, commitTag, idData);
  if (res < 0) {
    opserr << "PlaneStressSimplifiedJ2::sendSelf -- could not send ID" << endln;
    return res;
  }
  
  static Vector data(7);

  data(0) = CsavedStrain33;
  for (int i = 0; i < 3; i++)
    data(1+i) = Cstress(i);
  for (int i = 0; i < 3; i++)
    data(1+3+i) = Cstrain(i);  

  res = theChannel.sendVector(dbTag, commitTag, data);
  if (res < 0) {
    opserr << "PlaneStressSimplifiedJ2::sendSelf -- could not send Vector" << endln;
    return res;
  }

  res = the3DMaterial->sendSelf(commitTag, theChannel);
  if (res < 0) { 
    opserr << "PlaneStressMaterial::sendSelf() - failed to send vector material" << endln;
    return res;
  }
  
  return res;
};  

int PlaneStressSimplifiedJ2::recvSelf(int commitTag, Channel &theChannel,
				      FEM_ObjectBroker &theBroker)
{
  int res = 0;
  int dbTag = this->getDbTag();

  static ID idData(4);
  res = theChannel.recvID(dbTag, commitTag, idData);
  if (res < 0) {
    opserr << "PlasticDamageConcrete3d::recvSelf -- could not receive ID" << endln;
    return res;
  }  

  this->setTag(idData(0));
  int matClassTag = idData(1);
  ndm = idData(3);

  // if the associated material has not yet been created or is of the wrong type
  // create a new material for recvSelf later
  if (the3DMaterial == 0 || the3DMaterial->getClassTag() != matClassTag) {
    if (the3DMaterial != 0)
      delete the3DMaterial;
    the3DMaterial = theBroker.getNewNDMaterial(matClassTag);
    if (the3DMaterial == 0) {
      opserr << "PlaneStressSimplifiedJ2::recvSelf() - failed to get a material of type: " << matClassTag << endln;
      return -1;
    }
  }
  the3DMaterial->setDbTag(idData(2));
  
  static Vector data(7);
  res = theChannel.recvVector(dbTag, commitTag, data);
  if (res < 0) {
    opserr << "PlaneStressSimplifiedJ2::recvSelf -- could not receive Vector" << endln;
    return res;
  }

  CsavedStrain33 = data(0);
  for (int i = 0; i < 3; i++)
    Cstress(i) = data(1+i);
  for (int i = 0; i < 3; i++)
    Cstrain(i) = data(1+3+i);

  res = the3DMaterial->recvSelf(commitTag, theChannel, theBroker);
  if (res < 0) { 
    opserr << "PlaneStressSimplifiedJ2::recvSelf() - failed to recv vector material" << endln;
    return res;
  }

  this->revertToLastCommit();
    
  return res;
};    
  
     
Response * PlaneStressSimplifiedJ2::setResponse (const char **argv, int argc, OPS_Stream &s){


  if (strcmp(argv[0],"stress") == 0 || strcmp(argv[0],"stresses") == 0)
		return new MaterialResponse(this, 1, stress);

  else if (strcmp(argv[0],"strain") == 0 || strcmp(argv[0],"strains") == 0)
		return new MaterialResponse(this, 2, strain);

  else if (strcmp(argv[0],"tangent") == 0 || strcmp(argv[0],"Tangent") == 0)
		return new MaterialResponse(this, 3, theTangent);

   else if (strcmp(argv[0],"strain33") == 0 || strcmp(argv[0],"Strain33") == 0)
		return new MaterialResponse(this, 4, savedStrain33 );

  else
		return 0;
	
}



int PlaneStressSimplifiedJ2::getResponse (int responseID, Information &matInfo){
		


	switch (responseID) {
		case -1:
			return -1;
		case 1:
			if (matInfo.theVector != 0)
				*(matInfo.theVector) =stress;
			return 0;

		case 2:
			if (matInfo.theVector != 0)
				*(matInfo.theVector) = strain;
			return 0;

		case 3:
			if (matInfo.theMatrix != 0)
				*(matInfo.theMatrix) = theTangent;
			return 0;

	 	case 4:
		  //if (matInfo.theDouble != 0)
			    matInfo.setDouble (savedStrain33);
			return 0;



		}
		
 

	return 0;
};

void PlaneStressSimplifiedJ2::Print(OPS_Stream &s, int flag){
	// -- to be implemented.
	return;
};


int PlaneStressSimplifiedJ2::setParameter(const char **argv, int argc, Parameter &param){
  // -- to be implemented.
  return 0;
};

int PlaneStressSimplifiedJ2::updateParameter(int responseID, Information &eleInformation){
  return 0;
};

