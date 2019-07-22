  
 
// Written: Quan Gu and Zhijian Qiu  
// Created: 2015/01/25 

//------------------------------------------

// stress and strain all defined as: elongation is negative; compression is positive (+)


//#include <math.h>  
#include <LinearCap.h>  
#include <Information.h>
#include <Parameter.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>  

#define ND_TAG_LinearCap  12345654342
#include <T2Vector.h>
#include <MaterialResponse.h>
#include <elementAPI.h>

#include <fstream>            // Quan Gu   2013 March   HK
using std::ofstream;          // Quan Gu   2013 March   HK
using std::ios;               // Quan Gu   2013 March   HK


Vector LinearCap::tempVector(6);
Matrix LinearCap::tempMatrix(6,6);

static int numLinearCap = 0;

void *
OPS_LinearCap(void) {

  if (numLinearCap == 0) {
    numLinearCap++;
    opserr << "LinearCap nDmaterial - Written: Quan Gu and Zhijian Qiu \n";
  }

  // Pointer to a uniaxial material that will be returned
  NDMaterial *theMaterial = 0;

  int numArgs = OPS_GetNumRemainingInputArgs();
  if (numArgs < 5) {
    opserr << "Want: nDMaterial LinearCap tag? ndm? rho? G? K? <theta? alpha? T? tol? >\n";
    return 0;	
  }

  int iData[2];
  double dData[7];

  // set the defaults
  double theta = 0.11;
  double alpha = 2.6614e7;
  double T = -2.0684e6; 
  double tol = 1.0e-10;

  dData[3] = theta;
  dData[4] = alpha;
  dData[5] = T;
  dData[6] = tol;

  int numData = 2;
  if (OPS_GetInt(&numData, iData) != 0) {
    opserr << "WARNING invalid nDMaterial LinearCap - problems reading first 2 integers" << endln;
    return 0;
  }
  numData = numArgs - 2;
  if (OPS_GetDouble(&numData, dData) != 0) {
    opserr << "WARNING invalid nDMaterial LinearCap - problems reading doubles" << endln;
    return 0;
  }

  theMaterial = new LinearCap(iData[0],
			      dData[1],
			      dData[2],
			      dData[0],
			      dData[3],
			      dData[4],
			      dData[5], 
			      iData[1],
			      dData[6]);
				      
  return theMaterial;
}

LinearCap::LinearCap( int    pTag,
		      double pG,
		      double pK,
		      double pRho,
		      double pTheta,
		      double pAlpha,
		      double pT,
		      int pNdm, 
		      double pTol_k)
		      
   :  NDMaterial(pTag,ND_TAG_LinearCap),CStrain(6),CPlastStrain(6),CStress(6),
			 strain(6),plastStrain(6), stress(6), stressDev(6), theTangent(6,6)
			 
 {
	shearModulus = pG;
	bulkModulus = pK;
	rho = pRho;
	theta = pTheta;
	alpha = pAlpha;
	T = pT;
// -------------2012,1,15 --------
	if (T>0) T = -T;

	ndm = pNdm;
	tol_k = pTol_k; 
	stressI1 =0.0;

	flag =1;
// --
	revertToStart();

	debug =0;
    SHVs =0;
    parameterID =0;

// -- theMode 
	theMode =-10;


};


 
LinearCap::LinearCap( const LinearCap & a)
 : NDMaterial(a.getTag(),ND_TAG_LinearCap),CStrain(6),CPlastStrain(6),CStress(6),
				 strain(6),plastStrain(6), stress(6), stressDev(6), theTangent(6,6)
				 
{

	shearModulus = a.shearModulus;
	bulkModulus = a.bulkModulus;
	rho = a.rho;
	theta = a.theta;
	alpha = a.alpha;
	T = a.T;
// -------------2012,1,15 --------
	if (T>0) T = -T;

	ndm = a.ndm;
	tol_k =a.tol_k; 
	stressI1 =0.0;

	flag =1;
// --
	revertToStart();

    SHVs =0;
    parameterID =0;

};

LinearCap::~LinearCap( ) {
	return;

};

double LinearCap::getRho(void) {
	return rho;
};

int LinearCap::setTrialStrain(const Vector &pStrain) {  // strain from element is eng. strain!
 
//  static Vector temp(6);
  if (ndm==3 && pStrain.Size()==6) 
    strain = -1.0*pStrain;
  else if (ndm==2 && pStrain.Size()==3) {
    strain[0] = -1.0*pStrain[0];
    strain[1] = -1.0*pStrain[1];
    strain[2] = 0.0;
    strain[3] = -1.0*pStrain[2];
    strain[4] = 0.0;
    strain[5] = 0.0;
  }
  else {
    opserr << "Fatal:LinearCap:: Material dimension is: " << ndm << endln;
    opserr << "But strain vector size is: " << pStrain.Size() << endln;
    exit(-1);
  }
  // ----- change to real strain instead of eng. strain

  for ( int i = 3; i<6; i++) 
	strain[i] /=2.0; 

  return 0;
	
};



int LinearCap::setTrialStrain(const Vector &pStrain, const Vector &r) {

	return setTrialStrain(pStrain);
	
};


int LinearCap::setTrialStrainIncr(const Vector &pStrainRate) {
 
 // ----- change to real strain instead of eng. strain
// ---- since all strain in material is the true strain, not eng.strain. 

	for (int i=0; i<3;i++) {
		tempVector(i) =  pStrainRate(i);
		tempVector(i+3) =  pStrainRate(i+3)/2.0;
	}
	
	if (ndm==3 && pStrainRate.Size()==6) 
       strain = CStrain-pStrainRate;

    else if (ndm==2 && pStrainRate.Size()==3) {
		strain[0] = CStrain[0]-pStrainRate[0];
		strain[1] = CStrain[1]-pStrainRate[1];
		strain[2] = 0.0;
		strain[3] = CStrain[2]-pStrainRate[2];
		strain[4] = 0.0;
		strain[5] = 0.0;
	}
	 else {
		opserr << "Fatal:LinearCap:: Material dimension is: " << ndm << endln;
		opserr << "But strain vector size is: " << pStrainRate.Size() << endln;
		exit(-1);
	  }



  return 0;


};


int LinearCap::setTrialStrainIncr(const Vector &pStrainRate, const Vector &r) {
	return setTrialStrainIncr(pStrainRate);
};


const Vector & LinearCap::getStrain(void) {       

	if (ndm==3){
		tempVector = -1.0*strain;
		return tempVector;
	
	}

  else {
    static Vector workV(3);   //, temp6(6);
    workV[0] = -1.0*strain[0];
    workV[1] = -1.0*strain[1];
    workV[2] = -1.0*strain[3];
    return workV;
  }



};

int LinearCap::commitState(void)  {
	//if (theMode!=4)
	//opserr<<" stress is:"<<stress<<endln;  //debug only!

	CStrain = strain;
	CStress = stress;
	CPlastStrain = plastStrain;

	return 0;

};

int LinearCap::revertToLastCommit(void)  {
	return 0;
};

int LinearCap::revertToStart(void)  {
	
	CStrain.Zero();
	CPlastStrain.Zero();
	CStress.Zero();
	strain.Zero();
	plastStrain.Zero();
	stress.Zero();

 	return 0;
};

NDMaterial * LinearCap::getCopy(void)  {
	
    LinearCap * copy = new LinearCap(*this);
	return copy;

};

NDMaterial * LinearCap::getCopy(const char *code)  {
    LinearCap * copy = new LinearCap(*this);
	return copy;
};

const char * LinearCap::getType(void) const  {
     return (ndm == 2) ? "PlaneStrain" : "ThreeDimensional";
};

int LinearCap::getOrder(void) const  {
	 return (ndm == 2) ? 3 : 6;
};

int LinearCap::sendSelf(int commitTag, Channel &theChannel)  {return 0;};

int LinearCap::recvSelf(int commitTag, Channel &theChannel,   
		   FEM_ObjectBroker &theBroker )  {return 0;};


Response * LinearCap::setResponse (const char **argv, int argc, OPS_Stream &matInformation)  {

  if (strcmp(argv[0],"stress") == 0 || strcmp(argv[0],"stresses") == 0)
		return new MaterialResponse(this, 1, stress);

  else if (strcmp(argv[0],"strain") == 0 || strcmp(argv[0],"strains") == 0)
		return new MaterialResponse(this, 2, strain);

  else if (strcmp(argv[0],"tangent") == 0 || strcmp(argv[0],"Tangent") == 0)
		return new MaterialResponse(this, 3, theTangent);

    
  else if (strcmp(argv[0],"plasticStrain") == 0 || strcmp(argv[0],"plasticStrains") == 0)
		return new MaterialResponse(this, 4, plastStrain);

 	
		return 0;
	

};

int LinearCap::getResponse (int responseID, Information &matInfo)  {
	
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
			if (matInfo.theVector != 0)
				*(matInfo.theVector) = plastStrain;
			return 0;

		}

	return 0;

};

void LinearCap::Print(OPS_Stream &s, int flag)  {return;};


// ------------------------------------------
  

// --------------
double LinearCap::failureEnvelop(double I){

	return alpha+theta*I; 
};


// --------------
double LinearCap::failureEnvelopDeriv(double I){
	return theta;
};



int LinearCap::findMode(double normS, double I1){


	int mode =-1;

	if ((I1 <= T) && (normS <= failureEnvelop(T)))
		mode =1;
	else if ((I1 <= T) && (normS >= failureEnvelop(T)) && (normS <= failureEnvelop(T)+2.0*shearModulus/(9.0*bulkModulus)*(T-I1)/failureEnvelopDeriv(T)))
	    mode =2;

	else if ((normS <= failureEnvelop(I1)) && (I1 >= T) )
        mode = 4;
	else if (normS >= failureEnvelop(T)+2.0*shearModulus/(9.0*bulkModulus)*(T-I1)/failureEnvelopDeriv(T)){

		mode = 3;
	}
	return mode;
};


const Matrix & LinearCap::getInitialTangent(void) { 
	return getTangent();
};



const Vector & LinearCap::getStress(void) {	


    double CPlastStrainI1 = CPlastStrain(0)+CPlastStrain(1)+CPlastStrain(2);
	Vector CPlastStrainDev = CPlastStrain;
	Vector unitVector2(6);
	for (int   i=0; i<3; i++) {
		unitVector2(i)=1.0;
		unitVector2(i+3)=0.0;
	}
    CPlastStrainDev.addVector(1.0, unitVector2, -CPlastStrainI1/3.0);
	


	double strainI1 = strain(0)+strain(1)+strain(2); 
	Vector strainDev = strain; 
	strainDev.addVector(1.0, unitVector2, -1.0*strainI1/3.0);

	Vector TStressDev = 2.0 * shearModulus * (strainDev - CPlastStrainDev);

	double TStressI1 = 3.0*bulkModulus * (strainI1 - CPlastStrainI1);

	double normTS = pow( TStressDev && TStressDev, 0.5); 
	double normS =0.0;

	theMode  = findMode(normTS, TStressI1);    // used everywhere in this code!!
   
	//if (theMode!=4) 
	//opserr<<" this mode = "<<theMode<<endln;

//	theMode = mode;  //  为了保存mode到成员函数，为了getTangent()调试用。

	double deltGammar1 =0.0;
	double deltGammar2 =0.0;

	Vector deltPlastStrainDev(6);
//-------------------------------------------

    Matrix dTstressDevdStrain(6,6);
	Vector dTstressI1dStrain(6);

	Matrix dstressDevdStrain(6,6);
	Vector dstressI1dStrain(6);

    Matrix I_dev(6,6);
	I_dev.Zero();

	for (int i=0; i<6; i++)  	I_dev(i,i) =1.0;

	for (int  i=0; i<3; i++)
		for (int j=0; j<3; j++)
			I_dev(i,j) -= 1.0/3.0;

    Vector I2(6);   // unit vector order 2
    I2.Zero();

	 for (int  i=0; i<3; i++)  
		 I2(i) = 1.0;

    dTstressDevdStrain.addMatrix(0.0, I_dev, 2.0 * shearModulus);
	dTstressI1dStrain.addVector(0.0, I2,3*bulkModulus);


	//-------------------------------
    if (theMode ==1) {
		deltGammar2 = (T-TStressI1)/(9.0*bulkModulus);
		this->stressI1 = T;
		stressDev = TStressDev;
		deltPlastStrainDev.Zero();
		deltPlastStrainI1 = -3.0*deltGammar2;

//-------------------------tangent-------------
		dstressDevdStrain = dTstressDevdStrain;
		dstressI1dStrain.Zero();


	} // theMode ==1
	
	else if (theMode ==2){

		deltGammar1 = (normTS-failureEnvelop(T))/(2.0*shearModulus);
		deltGammar2 = (T-TStressI1)/(9.0*bulkModulus) - deltGammar1*failureEnvelopDeriv(T); 
		this->stressI1 = T;
		normS =failureEnvelop(T);
		stressDev.addVector(0.0, TStressDev, normS/ normTS);
		deltPlastStrainDev.addVector( 0.0, TStressDev, deltGammar1/normTS);
		deltPlastStrainI1 = -3.0 * (deltGammar1 * failureEnvelopDeriv(T)+deltGammar2);

//-------------------------tangent-------------
		dstressI1dStrain.Zero();  

        Vector tmpVector(6);
		Matrix dndStrain(6,6); 
		Matrix tmpMatrix(6,6); 

        tmpVector.addMatrixVector(0.0,dTstressDevdStrain,TStressDev,1);

		for(int i=0; i<6; i++){
			 
			for ( int j=0; j<3; j++)
				tmpMatrix(i,j) = TStressDev(i)*tmpVector(j);
			for ( int  j=3; j<6; j++)
				tmpMatrix(i,j) = TStressDev(i)*tmpVector(j)*2.0;     // To be consistent with the transformation between 4th order tensor and matrix
	
		}

		dndStrain.addMatrix(0.0, tmpMatrix, -1/normTS/normTS/normTS);	 // -TStressDev @ TStressDev : dTstressDevdStrain /normTS/normTS/normTS   


		dndStrain.addMatrix(1.0, dTstressDevdStrain, 1/normTS);	 // TStressDev @ TStressDev : dTstressDevdStrain    
		

		dstressDevdStrain.addMatrix(1.0, dndStrain, normS);
}
 
	else if (theMode ==4) {
		normS = normTS;
		stressDev = TStressDev;
		stressI1 = TStressI1;
		deltPlastStrainDev.Zero();
		deltPlastStrainI1 = 0.0;

//-------------------------tangent-------------

		Matrix elasticTangent(6,6);
        elasticTangent.Zero();

        for(int i=0; i<3; i++)
	        for (int j=0; j<3; j++)
		       elasticTangent(i,j) = bulkModulus-2.0/3.0*shearModulus;     

        for(int  i=0; i<6; i++)
		       elasticTangent(i,i) += 2.0*shearModulus;       //定义出弹性矩阵
	    theTangent = elasticTangent;

	}
	
	else if (theMode ==3){
		
		deltGammar1 = (normTS-failureEnvelop(TStressI1))/(2.0*shearModulus+9*bulkModulus*theta*theta);
		normS = normTS - 2.0*shearModulus * deltGammar1;

		stressI1 = 9*theta*bulkModulus*deltGammar1 + TStressI1;
		stressDev.addVector(0.0, TStressDev, normS/ normTS);
		deltPlastStrainDev.addVector(0.0, stressDev, deltGammar1/normS);
		deltPlastStrainI1 = -3.0*deltGammar1*failureEnvelopDeriv(stressI1);

//-------------------------tangent-------------
		Vector dDeltGammar1dStrain(6);
		Vector dTstressNormdStrain(6);
		Vector dNormSdStrain(6);
		Matrix dNdStrain(6,6);


        //elasticTangent.Zero();
	    dTstressNormdStrain.addMatrixVector(0.0,dTstressDevdStrain, TStressDev, 1.0/normTS);

        dDeltGammar1dStrain.addVector(0.0, dTstressNormdStrain, 1.0);		
        dDeltGammar1dStrain.addVector(1.0, dTstressI1dStrain, -1.0*theta);// 8.8 Qiu
        dDeltGammar1dStrain.addVector(0.0, dDeltGammar1dStrain, 1.0/(2*shearModulus + 9*bulkModulus*theta*theta));
        
		dNormSdStrain = dTstressNormdStrain;
        dNormSdStrain.addVector(1.0, dDeltGammar1dStrain, -2*shearModulus);

        dNdStrain = dTstressDevdStrain; 
	    dNdStrain.addMatrix(0.0,dTstressDevdStrain, 1.0/normTS);
	    tempVector.addMatrixVector(0.0, dTstressDevdStrain, TStressDev, 1.0/normTS/normTS/normTS);	
	    
		for ( int m=0; m<6; m++){	 
			for ( int n=0; n<3; n++)				 
				 dNdStrain(m,n) -=TStressDev(m)*tempVector(n);
			for ( int n=3; n<6; n++)				 
				 dNdStrain(m,n) -=TStressDev(m)*tempVector(n)*2.0;
		}

        dstressDevdStrain.addMatrix(0.0,dNdStrain, normS);

		for (int  m=0; m<6; m++){	 
			for ( int n=0; n<3; n++)				 
				 dstressDevdStrain(m,n) +=TStressDev(m)*dNormSdStrain(n)/normTS;
			for (int  n=3; n<6; n++)				 
				 dstressDevdStrain(m,n) +=TStressDev(m)*dNormSdStrain(n)*2.0/normTS;
		}

        dstressI1dStrain = 9 * bulkModulus * theta *dDeltGammar1dStrain;
        dstressI1dStrain.addVector(1.0,dTstressI1dStrain, 1);

	} //theMode ==3


// --- compute strain and stress --
                                            
	double plastStrainI1 = CPlastStrainI1 + deltPlastStrainI1;
	Vector plastStrainDev = CPlastStrainDev + deltPlastStrainDev;

	plastStrain.addVector(0.0,  plastStrainDev, 1.0);
	plastStrain.addVector(1.0, unitVector2, plastStrainI1/3.0);

    stress.addVector(0.0, stressDev,1.0);
	stress.addVector(1.0, unitVector2, stressI1/3.0);

// --- theTangent------------------
	if (theMode==4){
		
	for (int i=0; i<6; i++)
	    for (int j=3; j<6; j++)
		  theTangent (i,j) /=2.0;		
	}
	else {
    theTangent = dstressDevdStrain;
	 for (int m=0; m<6; m++){	 
		 for ( int n=0; n<3; n++)				 
			 theTangent(m,n) +=1./3*dstressI1dStrain(n)*I2(m);
		 for (int  n=3; n<6; n++)				 
			 theTangent(m,n) +=1./3*dstressI1dStrain(n)*I2(m)*2.0;
	 }

	for (int i=0; i<6; i++)
	    for (int j=3; j<6; j++)
		  theTangent (i,j) /=2.0;
	}


  //opserr<<"stress in getstress is "<<stress<<endln;
	if (ndm==3){
		tempVector.addVector(0.0, stress,-1.0);
		return tempVector;
	}

  else {
    static Vector workV(3);//, temp6(6);
    workV[0] = -1.0*stress[0];
    workV[1] = -1.0*stress[1];
    workV[2] = -1.0*stress[3];
    return workV;
  }
	
	

}
  
const Matrix & LinearCap::getTangent(void) { 
/*	
	
	static Vector tempStress(6);
	static Matrix compTangent(6,6);
    compTangent.Zero();

	this->getStress(); 
	tempStress = stress;
	
	// --- store strain

	static Vector strain_save(6);
	strain_save = strain;

	
	double perturb = 0.001;

	for (int  i=0; i<6; i++) {
	
		if (fabs(strain(i))<1.0e-14) 
			perturb =0.0000001;
		else 
			perturb =0.0001*strain(i);    

		strain(i) += perturb;
		this->getStress();


		for (int j=0; j<6; j++)
			compTangent(j,i) = (stress(j)-tempStress(j))/perturb;
	
		// --- recover strain ----
		strain(i)  = strain_save(i); 

	}

   strain  = strain_save;
	this->getStress();
 

  for ( i=0; i<6; i++)
	  for (int j=3; j<6; j++)
		  compTangent (i,j) /=2.0;

  if (theMode==2){
 
   opserr<<"compTangent in getTangent "<< compTangent<<endln;    
   opserr<<"theTangent in getTangent "<< theTangent<<endln;
   opserr<<"----------------------------------------------------"<<endln;
 
	  }


  this->getStress();


*/

    //theTangent=compTangent;
	  

//----------end debug ----------

  if (ndm==3) 
    return theTangent;
  else {
    static Matrix workM(3,3);
    workM(0,0) = theTangent(0,0);
    workM(0,1) = theTangent(0,1);
    workM(0,2) = theTangent(0,3);
    workM(1,0) = theTangent(1,0);
    workM(1,1) = theTangent(1,1);
    workM(1,2) = theTangent(1,3);
    workM(2,0) = theTangent(3,0);
    workM(2,1) = theTangent(3,1);
    workM(2,2) = theTangent(3,3);
    return workM;
  }
}



///////////////////////// add sensitivity ///////////////////////////////////
int 
LinearCap::setParameter (const char **argv, int argc, Parameter &param)

{	
	if (argc < 1)
		return -1;

	if (strcmp(argv[0],"G") == 0) {
		return param.addObject(1, this);
	}

	if (strcmp(argv[0],"K") == 0) {
		return param.addObject(2, this);
	}

	if (strcmp(argv[0],"rho") == 0) {
		return param.addObject(3, this);
	}

	if (strcmp(argv[0],"theta") == 0) {
		return param.addObject(4, this);
	}


	if (strcmp(argv[0],"alpha") == 0) {
		return param.addObject(5, this);
	}

	if (strcmp(argv[0],"T") == 0) {
		return param.addObject(6, this);
	}
	
	else
		opserr << "WARNING: Could not set parameter in CapPlasticity. " << endln;
                
	return -1;
}



int 
LinearCap::updateParameter(int passedParameterID,Information &info)

{
	switch (passedParameterID) {
	case -1:
		return -1;

	case 1:
		this->shearModulus= info.theDouble; //  G
		break;

	case 2:
		this->bulkModulus  = info.theDouble;  //K
		break;

	case 3:
		this->rho= info.theDouble;  // rho
		break;

	case 4:
		this->theta= info.theDouble; // theta
		break;

	case 5:
		this->alpha= info.theDouble; // alpha
		break;

	case 6:
		this->T  = info.theDouble;

		
		// -------------2012,1,15 --------
	    if (T>0) 
			this->T = -T;

		break;
 
  
	default:
		return -1;
	}
    
 
	return 0;
	}


int LinearCap::activateParameter(int passedParameterID){
	parameterID = passedParameterID;
	return 0;						
	}

const Vector &  LinearCap::getStressSensitivity (int gradNumber, bool conditional){
    int static ii=0;
	ii++;
//if (ii==2165)
//opserr<<"error"<<endln;
	double GSensitivity = 0.0;
	double KSensitivity = 0.0;
	double rhoSensitivity = 0.0;  
	double thetaSensitivity = 0.0;
	double alphaSensitivity = 0.0;
	double TSensitivity = 0.0;

	if (parameterID == 1) {  
		GSensitivity = 1.0;
	}
	else if (parameterID == 2) {   
		KSensitivity = 1.0;
	}
	else if (parameterID == 3) {   
		rhoSensitivity = 1.0;
	}
	else if (parameterID == 4) {    
		thetaSensitivity = 1.0;
	}
	else if (parameterID == 5) {   
		alphaSensitivity = 1.0;
	}
	else if (parameterID == 6) {   
		TSensitivity = 1.0;
	}
	else {

		// Nothing random here, but may have to return something in any case

	}

  
  // --- define history variables -----------
	static Vector CStrainSensitivity(6);	CStrainSensitivity.Zero();
	static Vector CStressSensitivity(6);   	CStressSensitivity.Zero();
	Vector CPlastStrainSensitivity(6);      CPlastStrainSensitivity.Zero();
 	Vector deltPlastStrainDev(6);

	static Vector stressSensitivity(6);
    stressSensitivity.Zero(); 

if (SHVs !=0) {
		for (int i=0;i<6;i++)	{
           CStrainSensitivity(i) = (*SHVs)(i,(gradNumber));
           CStressSensitivity(i)= (*SHVs)(i+6,(gradNumber));
           CPlastStrainSensitivity(i) = (*SHVs)(i+12,(gradNumber));
		
		}


	}


    double CPlastStrainI1 = CPlastStrain(0)+CPlastStrain(1)+CPlastStrain(2);

	Vector CPlastStrainDev = CPlastStrain;
	Vector unitVector2(6);
	for (int   i=0; i<3; i++) {
		unitVector2(i)=1.0;
		unitVector2(i+3)=0.0;
	}

    CPlastStrainDev.addVector(1.0, unitVector2, -CPlastStrainI1/3.0);
	
	double strainI1 = strain(0)+strain(1)+strain(2); 
	Vector strainDev = strain;
	strainDev.addVector(1.0, unitVector2, -1.0*strainI1/3.0);

	Vector TStressDev = 2.0 * shearModulus * (strainDev - CPlastStrainDev);
	double TStressI1 = 3.0*bulkModulus * (strainI1 - CPlastStrainI1);

// -------- sensitivity ----------

	double CPlastStrainI1Sensitivity = CPlastStrainSensitivity(0)+CPlastStrainSensitivity(1)+CPlastStrainSensitivity(2);

	Vector CPlastStrainDevSensitivity = CPlastStrainSensitivity;

    CPlastStrainDevSensitivity.addVector(1.0, unitVector2, -CPlastStrainI1Sensitivity/3.0);

	double strainI1Sensitivity = 0.0;                                // conditional sensitivity

	Vector strainDevSensitivity(6);   strainDevSensitivity.Zero();   // conditional sensitivity

	Vector TStressDevSensitivity(6);
	
	TStressDevSensitivity.addVector(0.0, strainDev, 2.0*GSensitivity);
	TStressDevSensitivity.addVector(1.0, CPlastStrainDev, -2.0*GSensitivity);
	TStressDevSensitivity.addVector(1.0, strainDevSensitivity, 2.0*shearModulus);
	TStressDevSensitivity.addVector(1.0, CPlastStrainDevSensitivity, -2.0*shearModulus);

    double TStressI1Sensitivity = 3.0*KSensitivity * (strainI1 - CPlastStrainI1)+3.0*bulkModulus * (strainI1Sensitivity - CPlastStrainI1Sensitivity);
 
	// Following are variables needed in each mode ..

	Vector deltPlastStrainDevSensitivity(6);   deltPlastStrainDevSensitivity.Zero();
	double deltPlastStrainI1Sensitivity = 0.0; 
	Vector stressDevSensitivity(6);  stressDevSensitivity.Zero(); 

	// ----- to be solved in each mode ...
    double stressI1Sensitivity =0.0; 
	double deltGammar1Sensitivity = 0.0; 	
	double deltGammar2Sensitivity = 0.0;	


// ========  theMode is judged already ========

	double normTS = pow( TStressDev && TStressDev, 0.5);    // cannot use norm() 
	double normTSSensitivity = (TStressDev && TStressDevSensitivity) / normTS;   

	double normS =0.0;
	double normSSensitivity = 0.0; 

	double deltGammar1 =0.0;
	double deltGammar2 =0.0;

if (theMode ==1) {

		deltGammar2 = (T-TStressI1)/(9.0*bulkModulus);
		this->stressI1 = T;
		stressDev = TStressDev;
		deltPlastStrainDev.Zero();
		deltPlastStrainI1 = -3.0*deltGammar2;

   //==================== adding sensitivity here ============================

        deltGammar2Sensitivity = ((TSensitivity-TStressI1Sensitivity)*bulkModulus-KSensitivity*(T-TStressI1))/9./bulkModulus/bulkModulus;
		stressI1Sensitivity = TSensitivity;
		stressDevSensitivity = TStressDevSensitivity;
		deltPlastStrainI1Sensitivity = -3.*deltGammar2Sensitivity;
		deltPlastStrainDevSensitivity.Zero();

}


else if (theMode ==2){
		deltGammar1 = (normTS-failureEnvelop(T))/(2.0*shearModulus);
		deltGammar2 = (T-TStressI1)/(9.0*bulkModulus) - deltGammar1*failureEnvelopDeriv(T); 
		this->stressI1 = T;
		normS =failureEnvelop(T);
		stressDev.addVector(0.0, TStressDev, normS/ normTS);
		deltPlastStrainDev.addVector( 0.0, TStressDev, deltGammar1/normTS);
		deltPlastStrainI1 = -3.0 * (deltGammar1 * failureEnvelopDeriv(T)+deltGammar2);

    //==================== adding sensitivity here ===================================
		double failureEnvelopSensitivity = alphaSensitivity+thetaSensitivity*T + theta*TSensitivity;

		double failureEnvelopDerivSensitivity = 0.0;
			   failureEnvelopDerivSensitivity = thetaSensitivity;

	    deltGammar1Sensitivity = ((normTSSensitivity-failureEnvelopSensitivity)*shearModulus-GSensitivity*(normTS-failureEnvelop(T)))/2./shearModulus/shearModulus;
       
		
		
		deltGammar2Sensitivity = ((TSensitivity-TStressI1Sensitivity)*bulkModulus-KSensitivity*(T-TStressI1))/9./bulkModulus/bulkModulus
			                    -(deltGammar1Sensitivity*failureEnvelopDeriv(T)+deltGammar1*failureEnvelopDerivSensitivity);

		stressI1Sensitivity = TSensitivity;


		Vector nSensitivity(6);
			nSensitivity.addVector(0.0, TStressDevSensitivity,1.0/normTS);
			nSensitivity.addVector(1.0, TStressDev, -normTSSensitivity/normTS/normTS);

		stressDevSensitivity.addVector(0.0, TStressDev, failureEnvelopSensitivity/normTS);
		stressDevSensitivity.addVector(1.0, nSensitivity, failureEnvelop(T));

/*
		stressDevSensitivity.addVector(0.0, nSensitivity, normS);
		stressDevSensitivity.addVector(1.0, TStressDev, normSSensitivity/normTS);
*/		
		
            
		deltPlastStrainDevSensitivity.addVector(0.0, TStressDev, deltGammar1Sensitivity/normTS);
		deltPlastStrainDevSensitivity.addVector(1.0, nSensitivity, deltGammar1);
		deltPlastStrainI1Sensitivity=-3.*(deltGammar1Sensitivity*failureEnvelopDeriv(T)+deltGammar1*failureEnvelopDerivSensitivity+deltGammar2Sensitivity);


	} //theMode ==2

else if (theMode ==3){
		
		deltGammar1 = (normTS-failureEnvelop(TStressI1))/(2.0*shearModulus+9*bulkModulus*theta*theta);
		normS = normTS - 2.0*shearModulus * deltGammar1;

		stressI1 = 9*theta*bulkModulus*deltGammar1 + TStressI1;
		stressDev.addVector(0.0, TStressDev, normS/ normTS);
		deltPlastStrainDev.addVector(0.0, stressDev, deltGammar1/normS);
		deltPlastStrainI1 = -3.0*deltGammar1*failureEnvelopDeriv(stressI1);

    //==================== adding sensitivity here ===================================

        deltGammar1Sensitivity = (normTSSensitivity - alphaSensitivity - TStressI1Sensitivity * theta -TStressI1 * thetaSensitivity)/(2*shearModulus+9*bulkModulus*theta*theta)
			                   - (normTS-failureEnvelop(TStressI1))*(2*GSensitivity+9*KSensitivity*theta*theta+18*bulkModulus*theta*thetaSensitivity)/(2*shearModulus+9*bulkModulus*theta*theta)/(2*shearModulus+9*bulkModulus*theta*theta);
 
		normSSensitivity = normTSSensitivity-2.*GSensitivity*deltGammar1-2.*shearModulus*deltGammar1Sensitivity;  
		Vector nSensitivity(6);
		nSensitivity.addVector(0.0, TStressDevSensitivity,1.0/normTS);
		nSensitivity.addVector(1.0, TStressDev, -normTSSensitivity/normTS/normTS);
		stressI1Sensitivity = 9*bulkModulus*theta*deltGammar1Sensitivity + 9*KSensitivity*theta*deltGammar1+9*bulkModulus*thetaSensitivity*deltGammar1+TStressI1Sensitivity;

		stressDevSensitivity.addVector(0.0, TStressDev, normSSensitivity/normTS);
		stressDevSensitivity.addVector(1.0, nSensitivity, normS);
		deltPlastStrainDevSensitivity.addVector(0.0, TStressDev, deltGammar1Sensitivity/normTS);
		deltPlastStrainDevSensitivity.addVector(1.0, nSensitivity, deltGammar1);

        double failureEnvelopDerivSensitivity = thetaSensitivity;

        deltPlastStrainI1Sensitivity=-3.*(deltGammar1Sensitivity*failureEnvelopDeriv(stressI1)+failureEnvelopDerivSensitivity*deltGammar1);


  }
  
else if (theMode ==4){

		normS = normTS;
		stressDev = TStressDev;
		stressI1 = TStressI1;
		deltPlastStrainDev.Zero();
		deltPlastStrainI1 = 0.0;

      //==================== adding sensitivity here ===================================

		stressDevSensitivity = TStressDevSensitivity;  
		stressI1Sensitivity = TStressI1Sensitivity;
		deltPlastStrainDevSensitivity.Zero();
		deltPlastStrainI1Sensitivity = 0.0;
  

}

    double plastStrainI1 = CPlastStrainI1 + deltPlastStrainI1;
	Vector plastStrainDev = CPlastStrainDev + deltPlastStrainDev;

	plastStrain.addVector(0.0,  plastStrainDev, 1.0);
	plastStrain.addVector(1.0, unitVector2, plastStrainI1/3.0);

    stress.addVector(0.0, stressDev,1.0);
	stress.addVector(1.0, unitVector2, stressI1/3.0);
  
	stressSensitivity.addVector(0.0, stressDevSensitivity,1.0);
    stressSensitivity.addVector(1.0, unitVector2,stressI1Sensitivity/3.0);  ///  qiuzhijian--getstressentivity
   // if (theMode!=4)
     //opserr<<" stressSensitivity is "<<stressSensitivity<<endln;

	if (ndm==3){
		tempVector.addVector(0.0, stressSensitivity,-1.0);
		return tempVector;
	}

  else {
    static Vector workV(3);//, temp6(6);
    workV[0] = -1.0*stressSensitivity[0];
    workV[1] = -1.0*stressSensitivity[1];
    workV[2] = -1.0*stressSensitivity[3];
    return workV;
	
  }

	return 0;  
}
int LinearCap::commitSensitivity(const Vector & strainGradient, int gradNumber, int numGrads){
 
	double GSensitivity = 0.0;
	double KSensitivity = 0.0;
	double rhoSensitivity = 0.0;
	double thetaSensitivity = 0.0;	
	double alphaSensitivity = 0.0;
	double TSensitivity = 0.0;

	static Vector stressSensitivity(6);
    stressSensitivity.Zero();

    static Vector strainSensitivity(6);

    if (ndm==3 && strainGradient.Size()==6) 
         strainSensitivity = strainGradient;
    else if (ndm==2 && strainGradient.Size()==3) {
		strainSensitivity[0] = strainGradient[0];
		strainSensitivity[1] = strainGradient[1];
		strainSensitivity[2] = 0.0;
		strainSensitivity[3] = strainGradient[2];
		strainSensitivity[4] = 0.0;
		strainSensitivity[5] = 0.0;

	}
    else {
		opserr << "Fatal:CapPlasticity:: Material dimension is: " << ndm << endln;
		opserr << "But strain vector size is: " << strainGradient.Size() << endln;
		exit(-1);
	}

	strainSensitivity *=-1.0;
	for (int j=3; j<6;j++)
		strainSensitivity[j] /=2.0;
	Vector plastStrainSensitivity(6); 	plastStrainSensitivity.Zero();
	//sour tstsopur teif (theMode!=4){
	//opserr<<"strainSensitivity"<<strainSensitivity<<endln;
	//}

	if (parameterID == 1) {  
		GSensitivity = 1.0;
	}
	else if (parameterID == 2) {   
		KSensitivity = 1.0;
	}
	else if (parameterID == 3) {   
		rhoSensitivity = 1.0;
	}
	else if (parameterID == 4) {    
		thetaSensitivity = 1.0;
	}
	else if (parameterID == 5) {   
		alphaSensitivity = 1.0;
	}
	else if (parameterID == 6) {   
		TSensitivity = 1.0;
	}
	else {

		// Nothing random here, but may have to return something in any case

	}

// --- define history variables -----------
	static Vector CStrainSensitivity(6);	    CStrainSensitivity.Zero();
	static Vector CStressSensitivity(6);   	    CStressSensitivity.Zero();
	static Vector CPlastStrainSensitivity(6);   CPlastStrainSensitivity.Zero();

  	if (SHVs ==0) {

       SHVs = new Matrix(18,numGrads);
	   SHVs->Zero();
	}

	else{

			for (int i=0;i<6;i++)	{
           CStrainSensitivity(i) = (*SHVs)(i,(gradNumber));
           CStressSensitivity(i)= (*SHVs)(i+6,(gradNumber));
           CPlastStrainSensitivity(i) = (*SHVs)(i+12,(gradNumber));

			}

	}





// ----- follow the stress computational process ----------

	double CPlastStrainI1 = CPlastStrain(0)+CPlastStrain(1)+CPlastStrain(2);

	Vector CPlastStrainDev = CPlastStrain;

	Vector unitVector2(6);
	for ( int i=0; i<3; i++) {
		unitVector2(i)=1.0;
		unitVector2(i+3)=0.0;  
	}
    CPlastStrainDev.addVector(1.0, unitVector2, -CPlastStrainI1/3.0);
	
	double strainI1 = strain(0)+strain(1)+strain(2); 

	Vector strainDev = strain;
	Vector deltPlastStrainDev(6); 
	
	strainDev.addVector(1.0, unitVector2, -1.0*strainI1/3.0);

	Vector TStressDev = 2.0 * shearModulus * (strainDev - CPlastStrainDev);
 
	double TStressI1 = 3.0*bulkModulus * (strainI1 - CPlastStrainI1);

// -------- sensitivity ----------

	double CPlastStrainI1Sensitivity = CPlastStrainSensitivity(0)+CPlastStrainSensitivity(1)+CPlastStrainSensitivity(2);

	Vector CPlastStrainDevSensitivity = CPlastStrainSensitivity;
    CPlastStrainDevSensitivity.addVector(1.0, unitVector2, -CPlastStrainI1Sensitivity/3.0);

	double strainI1Sensitivity = strainSensitivity(0)+strainSensitivity(1)+strainSensitivity(2);     // unconditional sensitivity

	Vector strainDevSensitivity(6);   
	       strainDevSensitivity.addVector(0.0, strainSensitivity, 1.0);   // unconditional sensitivity
	       strainDevSensitivity.addVector(1.0, unitVector2, -1.0/3*strainI1Sensitivity);   // unconditional sensitivity

	Vector TStressDevSensitivity(6);
	
	TStressDevSensitivity.addVector(0.0, strainDev, 2.0*GSensitivity);
	TStressDevSensitivity.addVector(1.0, CPlastStrainDev, -2.0*GSensitivity);

	TStressDevSensitivity.addVector(1.0, strainDevSensitivity, 2.0*shearModulus);
	TStressDevSensitivity.addVector(1.0, CPlastStrainDevSensitivity, -2.0*shearModulus);

    double TStressI1Sensitivity = 3.0*KSensitivity * (strainI1 - CPlastStrainI1)+3.0*bulkModulus * (strainI1Sensitivity - CPlastStrainI1Sensitivity);
 
	// Following are variables needed in each mode ..
	Vector deltPlastStrainDevSensitivity(6);   deltPlastStrainDevSensitivity.Zero();
	double deltPlastStrainI1Sensitivity = 0.0; 
	Vector stressDevSensitivity(6);  stressDevSensitivity.Zero(); 

	// ----- to be solved in each mode ...
    double stressI1Sensitivity =0.0; 
	double hardening_kSensitivity =0.0;
	double deltGammar1Sensitivity = 0.0; 	
	double deltGammar2Sensitivity = 0.0;	


// ========  theMode is judged already ========

	double normTS = pow( TStressDev && TStressDev, 0.5);   
	double normTSSensitivity = (TStressDev && TStressDevSensitivity) / normTS;   

	double normS =0.0;
	double normSSensitivity = 0.0; 

	double deltGammar1 =0.0;
	double deltGammar2 =0.0;

	if (theMode ==1) {

		deltGammar2 = (T-TStressI1)/(9.0*bulkModulus);
		this->stressI1 = T;
		stressDev = TStressDev;
		deltPlastStrainDev.Zero();
		deltPlastStrainI1 = -3.0*deltGammar2;

   //==================== adding sensitivity here ============================

        deltGammar2Sensitivity = ((TSensitivity-TStressI1Sensitivity)*bulkModulus-KSensitivity*(T-TStressI1))/9./bulkModulus/bulkModulus;
		stressI1Sensitivity = TSensitivity;
		stressDevSensitivity = TStressDevSensitivity;
		deltPlastStrainI1Sensitivity = -3.*deltGammar2Sensitivity;
		deltPlastStrainDevSensitivity.Zero();
	}


else if (theMode ==2){

		deltGammar1 = (normTS-failureEnvelop(T))/(2.0*shearModulus);
		deltGammar2 = (T-TStressI1)/(9.0*bulkModulus) - deltGammar1*failureEnvelopDeriv(T); 
		this->stressI1 = T;
		normS =failureEnvelop(T);
		stressDev.addVector(0.0, TStressDev, normS/ normTS);
		deltPlastStrainDev.addVector( 0.0, TStressDev, deltGammar1/normTS);
		deltPlastStrainI1 = -3.0 * (deltGammar1 * failureEnvelopDeriv(T)+deltGammar2);

    //==================== adding sensitivity here ===================================
		double failureEnvelopSensitivity = alphaSensitivity+thetaSensitivity*T + theta*TSensitivity;

		double failureEnvelopDerivSensitivity = 0.0;
			   failureEnvelopDerivSensitivity = thetaSensitivity;

	    deltGammar1Sensitivity = ((normTSSensitivity-failureEnvelopSensitivity)*shearModulus-GSensitivity*(normTS-failureEnvelop(T)))/2./shearModulus/shearModulus;
       
		
		
		deltGammar2Sensitivity = ((TSensitivity-TStressI1Sensitivity)*bulkModulus-KSensitivity*(T-TStressI1))/9./bulkModulus/bulkModulus
			                    -(deltGammar1Sensitivity*failureEnvelopDeriv(T)+deltGammar1*failureEnvelopDerivSensitivity);

		stressI1Sensitivity = TSensitivity;


		Vector nSensitivity(6);
			nSensitivity.addVector(0.0, TStressDevSensitivity,1.0/normTS);
			nSensitivity.addVector(1.0, TStressDev, -normTSSensitivity/normTS/normTS);
/*
		stressDevSensitivity.addVector(0.0, nSensitivity, normS);
		stressDevSensitivity.addVector(1.0, TStressDev, normSSensitivity/normTS);
*/
		stressDevSensitivity.addVector(0.0, TStressDev, failureEnvelopSensitivity/normTS);
		stressDevSensitivity.addVector(1.0, nSensitivity, failureEnvelop(T));

		
            
		deltPlastStrainDevSensitivity.addVector(0.0, TStressDev, deltGammar1Sensitivity/normTS);
		deltPlastStrainDevSensitivity.addVector(1.0, nSensitivity, deltGammar1);

		deltPlastStrainI1Sensitivity=-3.*(deltGammar1Sensitivity*failureEnvelopDeriv(T)+deltGammar1*failureEnvelopDerivSensitivity+deltGammar2Sensitivity);


	} //theMode ==2

else if (theMode ==3){
		
		deltGammar1 = (normTS-failureEnvelop(TStressI1))/(2.0*shearModulus+9*bulkModulus*theta*theta);
		normS = normTS - 2.0*shearModulus * deltGammar1;

		stressI1 = 9*theta*bulkModulus*deltGammar1 + TStressI1;
		stressDev.addVector(0.0, TStressDev, normS/ normTS);
		deltPlastStrainDev.addVector(0.0, stressDev, deltGammar1/normS);
		deltPlastStrainI1 = -3.0*deltGammar1*failureEnvelopDeriv(stressI1);

    //==================== adding sensitivity here ===================================

        deltGammar1Sensitivity = (normTSSensitivity - alphaSensitivity - TStressI1Sensitivity * theta -TStressI1 * thetaSensitivity)/(2*shearModulus+9*bulkModulus*theta*theta)
			                   - (normTS-failureEnvelop(TStressI1))*(2*GSensitivity+9*KSensitivity*theta*theta+18*bulkModulus*theta*thetaSensitivity)/(2*shearModulus+9*bulkModulus*theta*theta)/(2*shearModulus+9*bulkModulus*theta*theta);
 
		normSSensitivity = normTSSensitivity-2.*GSensitivity*deltGammar1-2.*shearModulus*deltGammar1Sensitivity;  
		Vector nSensitivity(6);
		nSensitivity.addVector(0.0, TStressDevSensitivity,1.0/normTS);
		nSensitivity.addVector(1.0, TStressDev, -normTSSensitivity/normTS/normTS);
		stressI1Sensitivity = 9*bulkModulus*theta*deltGammar1Sensitivity + 9*KSensitivity*theta*deltGammar1+9*bulkModulus*thetaSensitivity*deltGammar1+TStressI1Sensitivity;

		stressDevSensitivity.addVector(0.0, TStressDev, normSSensitivity/normTS);
		stressDevSensitivity.addVector(1.0, nSensitivity, normS);
		deltPlastStrainDevSensitivity.addVector(0.0, TStressDev, deltGammar1Sensitivity/normTS);
		deltPlastStrainDevSensitivity.addVector(1.0, nSensitivity, deltGammar1);

        double failureEnvelopDerivSensitivity = thetaSensitivity;

        deltPlastStrainI1Sensitivity=-3.*(deltGammar1Sensitivity*failureEnvelopDeriv(stressI1)+failureEnvelopDerivSensitivity*deltGammar1);


  } 
else if (theMode ==4){

		normS = normTS;
		stressDev = TStressDev;
		stressI1 = TStressI1;
		deltPlastStrainDev.Zero();
		deltPlastStrainI1 = 0.0;

      //==================== adding sensitivity here ===================================

		stressDevSensitivity = TStressDevSensitivity;  
		stressI1Sensitivity = TStressI1Sensitivity;
		deltPlastStrainDevSensitivity.Zero();
		deltPlastStrainI1Sensitivity = 0.0;

}

    double plastStrainI1 = CPlastStrainI1 + deltPlastStrainI1;
	Vector plastStrainDev = CPlastStrainDev + deltPlastStrainDev;

	plastStrain.addVector(0.0,  plastStrainDev, 1.0);
	plastStrain.addVector(1.0, unitVector2, plastStrainI1/3.0);  
  
	double plastStrainI1Sensitivity = CPlastStrainI1Sensitivity + deltPlastStrainI1Sensitivity;
	Vector plastStrainDevSensitivity = CPlastStrainDevSensitivity + deltPlastStrainDevSensitivity;

	plastStrainSensitivity.addVector(0.0,  plastStrainDevSensitivity, 1.0);
	plastStrainSensitivity.addVector(1.0, unitVector2, plastStrainI1Sensitivity/3.0);

  

    stress.addVector(0.0, stressDev,1.0);
	stress.addVector(1.0, unitVector2, stressI1/3.0);

	stressSensitivity.addVector(0.0, stressDevSensitivity,1.0);
    stressSensitivity.addVector(1.0, unitVector2,stressI1Sensitivity/3.0); 


	for ( int i =0; i<6; i++) {
		(*SHVs)(i, (gradNumber)) = strainSensitivity(i);
        (*SHVs)(i+6,(gradNumber)) = stressSensitivity(i);
        (*SHVs)(i+12,(gradNumber)) = plastStrainSensitivity(i);
	}


	return 0;
}; 
