// $Revision: 1.1 $
// $Date: 2009-07-23 23:44:46 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/nD/soil/MultiYieldSurfaceClay.cpp,v $
                                                                        
// Written: ZHY
// Created: August 2000
// Consistent tangent and Sensitivity: Quan Gu, Jul. 2009
// Compared with PressureIndependMultiYield, this model use:
//                             1. consistent tangent instead of continuum tangent
//                             2. only plastic state ( no elast2plast transfer)
//                             3. original strain only ( setSubStrainRate() does not work)
//
// refer to "Finite element response sensitivity analysis of multi-yield-surface J2 
// plasticity model by direct differentiation method", Quan Gu, Joel P. Conte, Ahmed Elgamal
// and Zhaohui Yang, Computer Methods in Applied Mechanics and Engineering
// Volume 198, Issues 30-32, 1 June 2009, Pages 2272-2285 
//
//
// MultiYieldSurfaceClay.cpp
// -------------------
//
#include <math.h>
#include <stdlib.h>
#include <MultiYieldSurfaceClay.h>
#include <MultiYieldSurface.h>

#include <Information.h>
#include <Parameter.h>
#include <ID.h>
#include <MaterialResponse.h>
#include <string.h>
#include <elementAPI.h>

Matrix MultiYieldSurfaceClay::theTangent(6,6);
Matrix MultiYieldSurfaceClay::dTrialStressdStrain(6,6);    //classwide matrix
Matrix MultiYieldSurfaceClay::dContactStressdStrain(6,6);  //classwide matrix
Matrix MultiYieldSurfaceClay::dSurfaceNormaldStrain(6,6);  //classwide matrix
Vector MultiYieldSurfaceClay::dXdStrain(6);                // classwide Vector


Vector     MultiYieldSurfaceClay::temp6(6);    // classwide Vector
Vector     MultiYieldSurfaceClay::temp(6);     // classwide Vector
Vector     MultiYieldSurfaceClay::devia(6);    // classwide Vector


double delta(int i,int j);
 
T2Vector MultiYieldSurfaceClay::dCurrentStress;
T2Vector MultiYieldSurfaceClay::dTrialStress;
T2Vector MultiYieldSurfaceClay::dCurrentStrain;
T2Vector MultiYieldSurfaceClay::dSubStrainRate;
T2Vector MultiYieldSurfaceClay::dStrainRate;
T2Vector MultiYieldSurfaceClay::dContactStress;


T2Vector MultiYieldSurfaceClay::subStrainRate;
int MultiYieldSurfaceClay::matCount=0;
int* MultiYieldSurfaceClay::loadStagex=0;  //=0 if elastic; =1 if plastic
int* MultiYieldSurfaceClay::ndmx=0;  //num of dimensions (2 or 3)
double* MultiYieldSurfaceClay::rhox=0;
double* MultiYieldSurfaceClay::frictionAnglex=0;
double* MultiYieldSurfaceClay::peakShearStrainx=0;
double* MultiYieldSurfaceClay::refPressurex=0;
double* MultiYieldSurfaceClay::cohesionx=0;
double* MultiYieldSurfaceClay::pressDependCoeffx=0;
int*    MultiYieldSurfaceClay::numOfSurfacesx=0;
double* MultiYieldSurfaceClay::residualPressx=0;

void* OPS_MultiYieldSurfaceClay()
{
    const int numParam = 6;
    const int totParam = 10;

    int argc = OPS_GetNumRemainingInputArgs() + 2;

    char * arg[] = {"nd", "rho", "refShearModul", "refBulkModul",
		    "cohesi", "peakShearStra",
		    "frictionAng (=0)", "refPress (=100)", "pressDependCoe (=0.0)",
		    "numberOfYieldSurf (=20)"};
    if (argc < (3+numParam)) {
	opserr << "WARNING insufficient arguments\n";
	opserr << "Want: nDMaterial MultiYieldSurfaceClay tag? " << arg[0];
	opserr << "? "<< "\n";
	opserr << arg[1] << "? "<< arg[2] << "? "<< arg[3] << "? "<< "\n";
	opserr << arg[4] << "? "<< arg[5] << "? "<< arg[6] << "? "<< "\n";
	opserr << arg[7] << "? "<< arg[8] << "? "<< arg[9] << "? \n";
	return 0;
    }
    
    int tag;
    int numdata = 1;
    if (OPS_GetIntInput(&numdata, &tag) < 0) {
	opserr << "WARNING invalid MultiYieldSurfaceClay tag\n";
	return 0;
    }

    double param[10];
    param[6] = 0.0;
    param[7] = 100.;
    param[8] = 0.0;
    param[9] = 20;
    numdata = 10;
    if (OPS_GetDoubleInput(&numdata, param) < 0) {
	opserr << "WARNING invalid MultiYieldSurfaceClay double inputs\n";
	return 0;
    }

    static double * gredu = 0;
    // user defined yield surfaces
    if (param[9] < 0 && param[9] > -40) {
	param[9] = -int(param[9]);
	numdata = int(2*param[9]);
	gredu = new double[numdata];
	if (OPS_GetDoubleInput(&numdata, gredu) < 0) {
	    opserr << "WARNING invalid MultiYieldSurfaceClay double inputs\n";
	    return 0;
	}
    }

    MultiYieldSurfaceClay * temp =
	new MultiYieldSurfaceClay (tag, param[0], param[1], param[2],
				   param[3], param[4], param[5], param[6],
				   param[7], param[8], param[9], gredu);
    if (gredu != 0) {
	delete [] gredu;
	gredu = 0;
    }

    return temp;
}

MultiYieldSurfaceClay::MultiYieldSurfaceClay (int tag, int nd, 
							double r, double refShearModul,
							double refBulkModul, 
							double cohesi, double peakShearStra, 
							double frictionAng, double refPress, double pressDependCoe,
							int numberOfYieldSurf, double * gredu)
/* : NDMaterial(tag,ND_TAG_MultiYieldSurfaceClay), currentStress(),
   trialStress(), currentStrain(), strainRate(),consistentTangent(6,6)*/
   : NDMaterial(tag,ND_TAG_MultiYieldSurfaceClay),consistentTangent(6,6)
{

  if (nd !=2 && nd !=3) {
    opserr << "FATAL:MultiYieldSurfaceClay:: dimension error" << endln;
    opserr << "Dimension has to be 2 or 3, you give nd= " << nd << endln;
    exit(-1);
  }
  if (refShearModul <= 0) {
    opserr << "FATAL:MultiYieldSurfaceClay::MultiYieldSurfaceClay: refShearModulus <= 0" << endln;
    exit(-1);
  }
  if (refBulkModul <= 0) {
    opserr << "FATAL:MultiYieldSurfaceClay::MultiYieldSurfaceClay: refBulkModulus <= 0" << endln;
    exit(-1);
  }
  if (frictionAng < 0.) {
    opserr << "WARNING:MultiYieldSurfaceClay::MultiYieldSurfaceClay: frictionAngle < 0" << endln;
    opserr << "Will reset frictionAngle to zero." << endln;
    frictionAng = 0.;
  }
  if (frictionAng == 0. && cohesi <= 0. ) {
    opserr << "FATAL:MultiYieldSurfaceClay::MultiYieldSurfaceClay: frictionAngle && cohesion <= 0." << endln;
    exit(-1);
  }
  if (cohesi <= 0) {
    opserr << "WARNING:MultiYieldSurfaceClay::MultiYieldSurfaceClay: cohesion <= 0" << endln;
    opserr << "Will reset cohesion to zero." << endln;
    cohesi = 0.;
  }
  if (peakShearStra <= 0) {
    opserr << "FATAL:MultiYieldSurfaceClay::MultiYieldSurfaceClay: peakShearStra <= 0" << endln;
    exit(-1);
  }
  if (refPress <= 0) {
    opserr << "FATAL:MultiYieldSurfaceClay::MultiYieldSurfaceClay: refPress <= 0" << endln;
    exit(-1);
  }
  if (pressDependCoe < 0) {
    opserr << "WARNING:MultiYieldSurfaceClay::MultiYieldSurfaceClay: pressDependCoe < 0" << endln;
    opserr << "Will reset pressDependCoe to zero." << endln;
    pressDependCoe = 0.;
  }
  if (numberOfYieldSurf <= 0) {
    opserr << "WARNING:MultiYieldSurfaceClay::MultiYieldSurfaceClay: numberOfSurfaces <= 0" << endln;
    opserr << "Will use 10 yield surfaces." << endln;
    numberOfYieldSurf = 10;
  }
  if (numberOfYieldSurf > 100) {
    opserr << "WARNING:MultiYieldSurfaceClay::MultiYieldSurfaceClay: numberOfSurfaces > 100" << endln;
   // opserr << "Will use 100 yield surfaces." << endln;   // Quan 
   // numberOfYieldSurf = 100;
  }
  if (r < 0) {
    opserr << "WARNING:MultiYieldSurfaceClay::MultiYieldSurfaceClay: mass density < 0" << endln;
    opserr << "Will use rho = 0." << endln;
    r = 0.;
  }
// for sensitivity
	

	parameterID = 0;        
	SHVs = 0;
	myNumGrads=1;
	dCommittedMultiSurfaceSize=0;
//	dCommittedMultiSurfaceElastPlastModul=0;
	dCommittedMultiSurfacePlastModul=0;
	dMultiSurfaceCenter=0;
	dCommittedMultiSurfaceCenter=0;
//	dVolume=0.0;
	surfacesSensitivityMark=0;
//............more ..............


  if (matCount%20 == 0) {
     int * temp1 = loadStagex;
	 int * temp2 = ndmx;
	 double * temp3 = rhox;
     double * temp_6 = frictionAnglex;
     double * temp7 = peakShearStrainx;
     double * temp8 = refPressurex;
     double * temp9 = cohesionx;
     double * temp10 = pressDependCoeffx;
	 int * temp11 = numOfSurfacesx;
     double * temp12 = residualPressx;
     loadStagex = new int[matCount+20];
     ndmx = new int[matCount+20];
     rhox = new double[matCount+20];
     frictionAnglex = new double[matCount+20];
     peakShearStrainx = new double[matCount+20];
     refPressurex = new double[matCount+20];
	 cohesionx = new double[matCount+20];
     pressDependCoeffx = new double[matCount+20];
     numOfSurfacesx = new int[matCount+20];
     residualPressx = new double[matCount+20];

	 for (int i=0; i<matCount; i++) {
         loadStagex[i] = temp1[i];
		 ndmx[i] = temp2[i];
         rhox[i] = temp3[i];
         frictionAnglex[i] = temp_6[i];
         peakShearStrainx[i] = temp7[i];
         refPressurex[i] = temp8[i];
         cohesionx[i] = temp9[i];
         pressDependCoeffx[i] = temp10[i];
         numOfSurfacesx[i] = temp11[i];
         residualPressx[i] = temp12[i];
     }

	 if (matCount > 0) {
	     delete [] temp1; delete [] temp2; delete [] temp3;
	     delete [] temp_6; delete [] temp7; delete [] temp8; 
	     delete [] temp9; delete [] temp10; delete [] temp11; 
		 delete [] temp12;
     }
  }
// changed by guquan temporily !!!!!!! changed back  ----------
  ndmx[matCount] = nd; // we ignore 2d material       
//  ndmx[matCount] = 3;
// end changed by guquan ---------------------------------------

  loadStagex[matCount] = 0;   //default
  refShearModulus = refShearModul;
  refBulkModulus = refBulkModul;
  frictionAnglex[matCount] = frictionAng;
  peakShearStrainx[matCount] = peakShearStra;
  refPressurex[matCount] = -refPress;  //compression is negative
  cohesionx[matCount] = cohesi;
  pressDependCoeffx[matCount] = pressDependCoe;
  numOfSurfacesx[matCount] = numberOfYieldSurf;
  rhox[matCount] = r;

  e2p = 0;
  matN = matCount;
  matCount ++;

	theSurfaces = new MultiYieldSurface[numberOfYieldSurf+1]; //first surface not used, pointer array??
    committedSurfaces = new MultiYieldSurface[numberOfYieldSurf+1]; 
	activeSurfaceNum = committedActiveSurf = 0; 

  setUpSurfaces(gredu);  // residualPress is calculated inside.
  debugMarks=0;   // quan debug nov. 2005


  // === update to plastic now ==== 2009 July
  loadStagex[matN] = 1; 

}
   

MultiYieldSurfaceClay::MultiYieldSurfaceClay () 
 : NDMaterial(0,ND_TAG_MultiYieldSurfaceClay), 
   currentStress(), trialStress(), currentStrain(), 
  strainRate(), theSurfaces(0), committedSurfaces(0)
{
  //does nothing
  // === update to plastic now ==== 2009 July
  loadStagex[matN] = 1; 
}


MultiYieldSurfaceClay::MultiYieldSurfaceClay (const MultiYieldSurfaceClay & a)
 : NDMaterial(a.getTag(),ND_TAG_MultiYieldSurfaceClay), 
   currentStress(a.currentStress), trialStress(a.trialStress), 
  currentStrain(a.currentStrain), strainRate(a.strainRate),consistentTangent(6,6)
{
  matN = a.matN;
  e2p = a.e2p;
  refShearModulus = a.refShearModulus;
  refBulkModulus = a.refBulkModulus;
  

  int numOfSurfaces = numOfSurfacesx[matN];

  committedActiveSurf = a.committedActiveSurf;
  activeSurfaceNum = a.activeSurfaceNum; 

// for sensitivity
	

	parameterID = 0;        
	SHVs = 0;
	myNumGrads=1;
	dCommittedMultiSurfaceSize=0;
//	dCommittedMultiSurfaceElastPlastModul=0;
	dCommittedMultiSurfacePlastModul=0;
	dMultiSurfaceCenter=0;
	dCommittedMultiSurfaceCenter=0;
//	dVolume=0.0;
	surfacesSensitivityMark=0;
	debugMarks = a.debugMarks;
//............more ..............


  theSurfaces = new MultiYieldSurface[numOfSurfaces+1];  //first surface not used
  committedSurfaces = new MultiYieldSurface[numOfSurfaces+1];  
  for(int i=1; i<=numOfSurfaces; i++) {
    committedSurfaces[i] = a.committedSurfaces[i];  
    theSurfaces[i] = a.theSurfaces[i];  
  }
  
  // === update to plastic now ==== 2009 July
  loadStagex[matN] = 1; 
}


MultiYieldSurfaceClay::~MultiYieldSurfaceClay ()
{
  if (theSurfaces != 0) delete [] theSurfaces;
  if (committedSurfaces != 0) delete [] committedSurfaces;

// --------------- for sensitivity   ----------------------
  if (SHVs != 0)  	delete SHVs;    
  if (	surfacesSensitivityMark != 0)  delete [] surfacesSensitivityMark;
}


void MultiYieldSurfaceClay::elast2Plast(void)
{
  int loadStage = loadStagex[matN];
  double frictionAngle = frictionAnglex[matN];
  int numOfSurfaces = numOfSurfacesx[matN];

  if (loadStage != 1 || e2p == 1) return;
  e2p = 1;

  if (currentStress.volume() > 0. && frictionAngle > 0.) {
    //opserr << "WARNING:MultiYieldSurfaceClay::elast2Plast(): material in tension." << endln;
    currentStress.setData(currentStress.deviator(),0);
  }

  paramScaling();  // scale surface parameters corresponding to initial confinement

  // Active surface is 0, return
  if (currentStress.deviatorLength() == 0.) return;

  // Find active surface
  while (yieldFunc(currentStress, committedSurfaces, ++committedActiveSurf) > 0) {
    if (committedActiveSurf == numOfSurfaces) {
      //opserr <<"WARNING:MultiYieldSurfaceClay::elast2Plast(): stress out of failure surface"<<endln;
      deviatorScaling(currentStress, committedSurfaces, numOfSurfaces);
      initSurfaceUpdate();
      return;
    }
  } 
  committedActiveSurf--;
  initSurfaceUpdate();
}


int MultiYieldSurfaceClay::setTrialStrain (const Vector &strain)
{
  int ndm = ndmx[matN];

//  static Vector temp(6);
  if (ndm==3 && strain.Size()==6) 
    temp = strain;
  else if (ndm==2 && strain.Size()==3) {
    temp[0] = strain[0];
    temp[1] = strain[1];
    temp[2] = 0.0;
    temp[3] = strain[2];
    temp[4] = 0.0;
    temp[5] = 0.0;
  }
  else {
    opserr << "Fatal:D2PressDepMYS:: Material dimension is: " << ndm << endln;
    opserr << "But strain vector size is: " << strain.Size() << endln;
    exit(-1);
  }
	
  //strainRate.setData(temp-currentStrain.t2Vector(1),1);
  temp -= currentStrain.t2Vector(1);
  strainRate.setData(temp, 1);
	
  return 0;
}


int MultiYieldSurfaceClay::setTrialStrain (const Vector &strain, const Vector &rate)
{
  return setTrialStrain (strain);
}


int MultiYieldSurfaceClay::setTrialStrainIncr (const Vector &strain)
{
  int ndm = ndmx[matN];

//  static Vector temp(6);
  if (ndm==3 && strain.Size()==6) 
    temp = strain;
  else if (ndm==2 && strain.Size()==3) {
    temp[0] = strain[0];
    temp[1] = strain[1];
    temp[3] = strain[2];
  }
  else {
    opserr << "Fatal:D2PressDepMYS:: Material dimension is: " << ndm << endln;
    opserr << "But strain vector size is: " << strain.Size() << endln;
    exit(-1);
  }
  
  strainRate.setData(temp,1);
  return 0;
}


int MultiYieldSurfaceClay::setTrialStrainIncr (const Vector &strain, const Vector &rate)
{
  return setTrialStrainIncr(strain);
}

/*
const Matrix & MultiYieldSurfaceClay::getTangent (void)
{
  int loadStage = loadStagex[matN];
  int ndm = ndmx[matN];

  if (loadStage == 1 && e2p == 0) elast2Plast();

  if (loadStage!=1) {  //linear elastic
    for (int i=0;i<6;i++) 
      for (int j=0;j<6;j++) {
	theTangent(i,j) = 0.;
	if (i==j) theTangent(i,j) += refShearModulus;
	if (i<3 && j<3 && i==j) theTangent(i,j) += refShearModulus;
	if (i<3 && j<3) theTangent(i,j) += (refBulkModulus - 2.*refShearModulus/3.);
      }
  }
  else {
    double coeff;
    static Vector devia(6);
  

	if (activeSurfaceNum > 0) {
      devia = trialStress.deviator();
      devia -= theSurfaces[activeSurfaceNum].center();
	
      double size = theSurfaces[activeSurfaceNum].size();
      double plastModul = theSurfaces[activeSurfaceNum].modulus();
      coeff = 6.*refShearModulus*refShearModulus/(2.*refShearModulus+plastModul)/size/size; }

//		theTangent.addMatrix(0.0,consistentTangent,1.0); }


	 else 
			coeff = 0.;
    
		for (int i=0;i<6;i++) 
			for (int j=0;j<6;j++) {
			theTangent(i,j) = - coeff*devia[i]*devia[j];
			if (i==j) theTangent(i,j) += refShearModulus;
			if (i<3 && j<3 && i==j) theTangent(i,j) += refShearModulus;
			if (i<3 && j<3) theTangent(i,j) += (refBulkModulus - 2.*refShearModulus/3.);
			} //for
	 	 

  }

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
*/



const Matrix & MultiYieldSurfaceClay::getTangent (void)
{
  int loadStage = loadStagex[matN];
  int ndm = ndmx[matN];

  if (loadStage == 1 && e2p == 0) {
	  opserr << "FATAL:MultiYieldSurfaceClay::Can not deal with e2p" 
	       << endln; 
	  //exit(-1);
	double coeff;
	if (activeSurfaceNum > 0) {
      devia = trialStress.deviator();
      devia -= theSurfaces[activeSurfaceNum].center();
	
      double size = theSurfaces[activeSurfaceNum].size();
      double plastModul = theSurfaces[activeSurfaceNum].modulus();
      coeff = 6.*refShearModulus*refShearModulus/(2.*refShearModulus+plastModul)/size/size; 
	}
	else 
			coeff = 0.;
    
	for (int i=0;i<6;i++) 
		for (int j=0;j<6;j++) {
			theTangent(i,j) = - coeff*devia[i]*devia[j];
			if (i==j) theTangent(i,j) += refShearModulus;
			if (i<3 && j<3 && i==j) theTangent(i,j) += refShearModulus;
			if (i<3 && j<3) theTangent(i,j) += (refBulkModulus - 2.*refShearModulus/3.);
		
		} //for  
  
  
  } //  if (loadStage == 1 && e2p == 0) {


  if (loadStage!=1) {  //linear elastic
    {
	  opserr << "FATAL:MultiYieldSurfaceClay::can not deal with linear elastic" 
	       << endln; 
	  exit(-1);
	}
  }
  else {

//    double coeff;
 //   static Vector devia(6);
  
/*   guquan
	if (activeSurfaceNum > 0) {
      devia = trialStress.deviator();
      devia -= theSurfaces[activeSurfaceNum].center();
	
      double size = theSurfaces[activeSurfaceNum].size();
      double plastModul = theSurfaces[activeSurfaceNum].modulus();
      coeff = 6.*refShearModulus*refShearModulus/(2.*refShearModulus+plastModul)/size/size; }




	 else 
			coeff = 0.;
    
		for (int i=0;i<6;i++) 
			for (int j=0;j<6;j++) {
			theTangent(i,j) = - coeff*devia[i]*devia[j];
			if (i==j) theTangent(i,j) += refShearModulus;
			if (i<3 && j<3 && i==j) theTangent(i,j) += refShearModulus;
			if (i<3 && j<3) theTangent(i,j) += (refBulkModulus - 2.*refShearModulus/3.);
			} //for

//        opserr << "getTangent by zhaohui"<< theTangent << endln; 
 end guquan */
   
		theTangent.addMatrix(0.0,consistentTangent,1.0); 

//		opserr << "getTangent by Guquan"<< theTangent << endln; 

  }  //else

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




const Matrix & MultiYieldSurfaceClay::getInitialTangent (void)
{
  int ndm = ndmx[matN];

  for (int i=0;i<6;i++) 
    for (int j=0;j<6;j++) {
      theTangent(i,j) = 0.;
      if (i==j) theTangent(i,j) += refShearModulus;
      if (i<3 && j<3 && i==j) theTangent(i,j) += refShearModulus;
      if (i<3 && j<3) theTangent(i,j) += (refBulkModulus - 2.*refShearModulus/3.);
    }

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


const Vector & MultiYieldSurfaceClay::getStress (void)
{
  int loadStage = loadStagex[matN];
  int numOfSurfaces = numOfSurfacesx[matN];
  int ndm = ndmx[matN];

  int i;
  if (loadStage == 1 && e2p == 0) elast2Plast();

  if (loadStage!=1) {  //linear elastic
    //trialStrain.setData(currentStrain.t2Vector() + strainRate.t2Vector());
    getTangent();
    static Vector a(6);
    a = currentStress.t2Vector();
	a.addMatrixVector(1.0, theTangent, strainRate.t2Vector(1), 1.0);
    trialStress.setData(a);
  }

  else {
    for (i=1; i<=numOfSurfaces; i++) theSurfaces[i] = committedSurfaces[i];
    activeSurfaceNum = committedActiveSurf;
    subStrainRate = strainRate;
	// output strainRate for debug

//   opserr << "step0. strainRate is," << strainRate.t2Vector()<< endln;
//	opserr << "step0. deviator of strainRate is," << strainRate.deviator() << endln;

    setTrialStress(currentStress);
    if (isLoadReversal()) {
      updateInnerSurface();
      activeSurfaceNum = 0;
    }
    int numSubIncre = setSubStrainRate();
    
    for (i=0; i<numSubIncre; i++) {
      if (i==0)  
	setTrialStress(currentStress);
      else 
	setTrialStress(trialStress);

      if (activeSurfaceNum==0 && !isCrossingNextSurface()) continue;
      if (activeSurfaceNum==0) activeSurfaceNum++;
//	  opserr << "enter into plasticity "<< endln;
      stressCorrection(0);
      updateActiveSurface();
    }

	//--------- add consistent tangent part code -------------------
	// unitTensor = I*I 
	static Matrix unitTensor(6,6);
	static Matrix tempTangent(6,6);
	unitTensor.Zero();
	for(int i=0;i<3;i++){
		for(int j=0;j<3;j++){
			unitTensor(i,j)=1.0;
		}
	}
	// change derivation to deviatoric strain into deviation to totall strain
	//tempTangent=dTrialStressdStrain*unitTensor;    cause problem !

	doubledotMatrixProduct(tempTangent,dTrialStressdStrain,unitTensor);
	tempTangent /= -3.0;
	dTrialStressdStrain.addMatrix(1.0, tempTangent,1.0);

//	opserr << "step6. getStress, change into total strain trialStress is," << trialStress.deviator()<< endln;
//	opserr << "step6. getStress, change into total strain dTrialStressdStrain is," << dTrialStressdStrain<< endln;


    //--------------------------------------------------------------------


    //volume stress change
    double volum = refBulkModulus*(strainRate.volume()*3.);
    volum += currentStress.volume();
    // if (volum > 0) volum = 0.;


//	Vector temp(6);
	temp.addVector(0.0,trialStress.deviator(),1.0);
    trialStress.setData(temp,volum);


	// ----------- add consistent tangent part code ----------------------

	dTrialStressdStrain.addMatrix(1.0,unitTensor,refBulkModulus);

//	opserr << "step7. getStress, deal with volume, trialStress is," << trialStress.t2Vector()<< endln;
//	opserr << "step7. getStress, deal with volume, dTrialStressdStrain is," << dTrialStressdStrain<< endln;


	// save consistent tangent into consistangent
	consistentTangent.addMatrix(0.0,dTrialStressdStrain,1.0);
	
	/* recover to strain
	for ( i=0;i<6;i++){
		for(int j=3;j<6;j++){
			consistentTangent(i,j) *=2.;
		}
	}  */


  }

//	static Vector temp6(6);
	temp6.addVector(0.0,trialStress.t2Vector(),1.0);
  if (ndm==3)
	return temp6;
//	  return trialStress.t2Vector();
  else {
    static Vector workV(3);
    workV[0] = trialStress.t2Vector()[0];
    workV[1] = trialStress.t2Vector()[1];
    workV[2] = trialStress.t2Vector()[3];
    return workV;
  }
}



const Vector & MultiYieldSurfaceClay::getStrain (void)
{
  return getCommittedStrain();
}


int MultiYieldSurfaceClay::commitState (void)
{
  int loadStage = loadStagex[matN];
  int numOfSurfaces = numOfSurfacesx[matN];

  currentStress = trialStress;
  
  //currentStrain = T2Vector(currentStrain.t2Vector() + strainRate.t2Vector());
// static Vector temp(6);
  temp = currentStrain.t2Vector();
  temp += strainRate.t2Vector();
  currentStrain.setData(temp);
  temp.Zero();
  strainRate.setData(temp);
  
  if (loadStage==1) {
    committedActiveSurf = activeSurfaceNum;
//	opserr<<"committedActiveSurface is:"<<activeSurfaceNum<<endln;


    for (int i=1; i<=numOfSurfaces; i++) committedSurfaces[i] = theSurfaces[i];
  }

  return 0;
}
 

int MultiYieldSurfaceClay::revertToLastCommit (void)
{
  return 0;
}


int MultiYieldSurfaceClay::revertToStart (void)
{
    activeSurfaceNum = committedActiveSurf = 0; 
	currentStrain.Zero();
	currentStress.Zero();
	

    trialStress.Zero();
    strainRate.Zero();
    subStrainRate.Zero();
	devia.Zero();
	for (int i = 0; i<=numOfSurfacesx[matN]; i++){
		theSurfaces[i].setCenter(devia);
		committedSurfaces[i].setCenter(devia);
	}

// AddingSensitivity:BEGIN /////////////////////////////////
	if (SHVs != 0) 
		SHVs->Zero();
	// .... more .....
	parameterID = 0;        
/*
	dCommittedMultiSurfaceSize=0;
	dCommittedMultiSurfacePlastModul=0;
	dMultiSurfaceCenter=0;
	dCommittedMultiSurfaceCenter=0;
	dVolume=0.0;
*/	
	int numOfSurfaces = numOfSurfacesx[matN];
	for (int i=0; i<numOfSurfaces+1; i++){
		for(int j=0;j<myNumGrads;j++){
			if (dMultiSurfaceCenter !=0)
				for(int k=0;k<6;k++) dMultiSurfaceCenter[k+i*6+j*6*(numOfSurfaces+1)]=0.0;
			if (dCommittedMultiSurfaceSize !=0)
			    dCommittedMultiSurfaceSize[i+j*(numOfSurfaces+1)]=0.0;
			if (dCommittedMultiSurfacePlastModul !=0)
				dCommittedMultiSurfacePlastModul[i+j*(numOfSurfaces+1)]=0.0;
			for(int k=0;k<6;k++) 
				if (dCommittedMultiSurfaceCenter !=0)
					dCommittedMultiSurfaceCenter[k+i*6+j*6*(numOfSurfaces+1)]=0.0;
		}
	} //for

	if (surfacesSensitivityMark !=0)
		for (int i=0; i<myNumGrads; i++)
			surfacesSensitivityMark[i]=0;


// AddingSensitivity:END //////////////////////////////////

	return 0;
}


NDMaterial * MultiYieldSurfaceClay::getCopy (void)
{
  MultiYieldSurfaceClay * copy = new MultiYieldSurfaceClay(*this);

  return copy;
}


NDMaterial * MultiYieldSurfaceClay::getCopy (const char *code)
{
  if (strcmp(code,"MultiYieldSurfaceClay") == 0 || strcmp(code,"PlaneStrain") == 0
      || strcmp(code,"ThreeDimensional") == 0) {
    MultiYieldSurfaceClay * copy = new MultiYieldSurfaceClay(*this);
    return copy;
  }

  return 0;
}


const char * MultiYieldSurfaceClay::getType (void) const
{
  int ndm = ndmx[matN];

  return (ndm == 2) ? "PlaneStrain" : "ThreeDimensional";
}


int MultiYieldSurfaceClay::getOrder (void) const
{
  int ndm = ndmx[matN];

  return (ndm == 2) ? 3 : 6;
}




int MultiYieldSurfaceClay::sendSelf(int commitTag, Channel &theChannel)
{
  int loadStage = loadStagex[matN];
  int ndm = ndmx[matN];
  int numOfSurfaces = numOfSurfacesx[matN];
  double rho = rhox[matN];
  double frictionAngle = frictionAnglex[matN];
  double peakShearStrain = peakShearStrainx[matN];
  double refPressure = refPressurex[matN];
  double cohesion = cohesionx[matN];
  double pressDependCoeff = pressDependCoeffx[matN];
  double residualPress = residualPressx[matN];

  int i, res = 0;

  static ID idData(5);
  idData(0) = this->getTag();
  idData(1) = numOfSurfaces;
  idData(2) = loadStage;
  idData(3) = ndm;
  idData(4) = matN;

  res += theChannel.sendID(this->getDbTag(), commitTag, idData);
  if (res < 0) {
    opserr << "PressureDependMultiYield::sendSelf -- could not send ID\n";
    return res;
  }

  Vector data(23+numOfSurfaces*8);
//  static Vector temp(6);
  data(0) = rho;
  data(1) = refShearModulus;
  data(2) = refBulkModulus;
  data(3) = frictionAngle;
  data(4) = peakShearStrain;
  data(5) = refPressure;
  data(6) = cohesion;
  data(7) = pressDependCoeff;
  data(8) = residualPress;
  data(9) = e2p;
  data(10) = committedActiveSurf;
	
  temp = currentStress.t2Vector();
  for(i = 0; i < 6; i++) data(i+11) = temp[i];
  
  temp = currentStrain.t2Vector();
  for(i = 0; i < 6; i++) data(i+17) = temp[i];
  
  for(i = 0; i < numOfSurfaces; i++) {
    int k = 23 + i*8;
    data(k) = committedSurfaces[i+1].size();
    data(k+1) = committedSurfaces[i+1].modulus();
    temp = committedSurfaces[i+1].center();
    data(k+2) = temp(0);
    data(k+3) = temp(1);
    data(k+4) = temp(2);
    data(k+5) = temp(3);
    data(k+6) = temp(4);
    data(k+7) = temp(5);
  }

  res += theChannel.sendVector(this->getDbTag(), commitTag, data);
  if (res < 0) {
    opserr << "PressureDependMultiYield::sendSelf -- could not send Vector\n";
    return res;
  }
  
  return res;
}


int MultiYieldSurfaceClay::recvSelf(int commitTag, Channel &theChannel, 
					 FEM_ObjectBroker &theBroker)    
{
  int i, res = 0;

  static ID idData(5);

  res += theChannel.recvID(this->getDbTag(), commitTag, idData);
  if (res < 0) {
    opserr << "PressureDependMultiYield::recvSelf -- could not recv ID\n";
    return res;
  }

  this->setTag((int)idData(0));
  int numOfSurfaces = idData(1);
  int loadStage = idData(2);
  int ndm = idData(3);
  matN = idData(4);

  Vector data(23+idData(1)*8);
//  static Vector temp(6);
  
  res += theChannel.recvVector(this->getDbTag(), commitTag, data);
  if (res < 0) {
    opserr << "PressureDependMultiYield::recvSelf -- could not recv Vector\n";
    return res;
  }
    
  double rho = data(0);
  double refShearModulus = data(1);
  double refBulkModulus = data(2);
  double frictionAngle = data(3);
  double peakShearStrain = data(4);
  double refPressure = data(5);
  double cohesion = data(6);
  double pressDependCoeff = data(7);
  double residualPress = data(8);
  e2p = data(9);
  committedActiveSurf = data(10);
  
  for(i = 0; i < 6; i++) temp[i] = data(i+11);
  currentStress.setData(temp);
  
  for(i = 0; i < 6; i++) temp[i] = data(i+17);
  currentStrain.setData(temp);

  if (committedSurfaces != 0) {
    delete [] committedSurfaces;
    delete [] theSurfaces;
  }

  theSurfaces = new MultiYieldSurface[numOfSurfaces+1]; //first surface not used
  committedSurfaces = new MultiYieldSurface[numOfSurfaces+1]; 
  for (i=1; i<=numOfSurfaces; i++) {
    committedSurfaces[i] = MultiYieldSurface();
  }
  
  for(i = 0; i < numOfSurfaces; i++) {
    int k = 23 + i*8;
    temp(0) = data(k+2);
    temp(1) = data(k+3);
    temp(2) = data(k+4);
    temp(3) = data(k+5);
    temp(4) = data(k+6);
    temp(5) = data(k+7);
    committedSurfaces[i+1].setData(temp, data(k), data(k+1));
  }
  
  loadStagex[matN] = loadStage;
  ndmx[matN] = ndm;
  numOfSurfacesx[matN] = numOfSurfaces;
  rhox[matN] = rho;
  frictionAnglex[matN] = frictionAngle;
  peakShearStrainx[matN] = peakShearStrain;
  refPressurex[matN] = refPressure;
  cohesionx[matN] = cohesion;
  pressDependCoeffx[matN] = pressDependCoeff;
  residualPressx[matN] = residualPress;

  return res;
}


Response*
MultiYieldSurfaceClay::setResponse (const char **argv, int argc, OPS_Stream &theOutput)
{
  if (strcmp(argv[0],"stress") == 0 || strcmp(argv[0],"stresses") == 0)
		return new MaterialResponse(this, 1, this->getCommittedStress());

  else if (strcmp(argv[0],"strain") == 0 || strcmp(argv[0],"strains") == 0)
		return new MaterialResponse(this, 2, this->getCommittedStrain());
    
	else if (strcmp(argv[0],"tangent") == 0)
		return new MaterialResponse(this, 3, this->getTangent());
    
	else if (strcmp(argv[0],"backbone") == 0) {
	    int numOfSurfaces = numOfSurfacesx[matN];
        static Matrix curv(numOfSurfaces+1,(argc-1)*2);
		  for (int i=1; i<argc; i++)
		   	curv(0,(i-1)*2) = atoi(argv[i]);
		return new MaterialResponse(this, 4, curv);
	}
	else if (strcmp(argv[0],"stressSensitivity") == 0 || strcmp(argv[0],"stresssensitivity") == 0) {
		int gradientNum = atoi(argv[1]);
			return new MaterialResponse(this, gradientNum+100, this->getCommittedStressSensitivity(1));  // use only size of matrix

	}
	else if (strcmp(argv[0],"strainSensitivity") == 0 || strcmp(argv[0],"strainsensitivity") == 0) {
		int gradientNum = atoi(argv[1]);
		return new MaterialResponse(this, gradientNum+500, this->getCommittedStrainSensitivity(1));
	}
  // use only size of matrix

	






	else
		return 0;
}


int MultiYieldSurfaceClay::getResponse (int responseID, Information &matInfo)
{
	switch (responseID) {
		case -1:
			return -1;
		case 1:
			if (matInfo.theVector != 0)
				*(matInfo.theVector) = getCommittedStress();
			return 0;
		case 2:
			if (matInfo.theVector != 0)
				*(matInfo.theVector) = getCommittedStrain();
			return 0;
		case 3:
			if (matInfo.theMatrix != 0)
				*(matInfo.theMatrix) = getTangent();
			return 0;
		case 4:
			if (matInfo.theMatrix != 0) 
				getBackbone(*(matInfo.theMatrix));
			return 0;


// for sensitivity use Quan Gu Fe.2.05
		default:    // sensitivity
			if (responseID>100 && responseID<500) {
//				case responseID:   //G,k tao..

				if (matInfo.theVector != 0) 
					
					*(matInfo.theVector) = this->getCommittedStressSensitivity(responseID-100);
					return 0;

			} // if

			else if (responseID>500) {

//				case responseID:   //

					if (matInfo.theVector != 0) 
						*(matInfo.theVector) = this->getCommittedStrainSensitivity(responseID-500);
					return 0;

			} // if
			return -1;




	}
}


void MultiYieldSurfaceClay::getBackbone (Matrix & bb)
{
  double residualPress = residualPressx[matN];
  double refPressure = refPressurex[matN];
  double pressDependCoeff =pressDependCoeffx[matN];
  int numOfSurfaces = numOfSurfacesx[matN];

  double vol, conHeig, scale, factor, shearModulus, stress1, 
		     stress2, strain1, strain2, plastModulus, elast_plast, gre;

	for (int k=0; k<bb.noCols()/2; k++) {
		vol = bb(0,k*2);
		if (vol<=0.) {
			opserr <<k<< "\nNDMaterial " <<this->getTag()
			  <<": invalid confinement for backbone recorder, " << vol << endln;
			continue;
		}
		conHeig = vol + residualPress;
		scale = -conHeig / (refPressure-residualPress);
		factor = pow(scale, pressDependCoeff); 
		shearModulus = factor*refShearModulus;
		for (int i=1; i<=numOfSurfaces; i++) {
			if (i==1) {
				stress2 = theSurfaces[i].size()*factor/1.732050807568877;
				strain2 = stress2/shearModulus;
				bb(1,k*2) = strain2; bb(1,k*2+1) = shearModulus;
			} else {
				stress1 = stress2; strain1 = strain2;
				plastModulus = factor*theSurfaces[i-1].modulus();
				elast_plast = 2*shearModulus*plastModulus/(2*shearModulus+plastModulus);
				stress2 = factor*theSurfaces[i].size()/1.732050807568877;
			  strain2 = 2*(stress2-stress1)/elast_plast + strain1;
				gre = stress2/strain2;
        bb(i,k*2) = strain2; bb(i,k*2+1) = gre;
			}
		}
	}

}

void MultiYieldSurfaceClay::Print(OPS_Stream &s, int flag )
{
  s << "MultiYieldSurfaceClay" << endln;
}


const Vector & MultiYieldSurfaceClay::getCommittedStress (void)
{
	int ndm = ndmx[matN];
	int numOfSurfaces = numOfSurfacesx[matN];

	double scale = sqrt(3./2.)*currentStress.deviatorLength()/committedSurfaces[numOfSurfaces].size();
	if (loadStagex[matN] != 1) scale = 0.;
	if (ndm==3) {
		static Vector temp7(7);
//		static Vector temp6(6);
		temp6 = currentStress.t2Vector();
    temp7[0] = temp6[0];
    temp7[1] = temp6[1];
    temp7[2] = temp6[2];
    temp7[3] = temp6[3];
    temp7[4] = temp6[4];
    temp7[5] = temp6[5];
    temp7[6] = scale;
		return temp7;
	}
  else {
    static Vector temp3(3);
	//, temp6(6);
		temp6 = currentStress.t2Vector();
    temp3[0] = temp6[0];
    temp3[1] = temp6[1];
//    temp5[2] = temp6[2];
    temp3[2] = temp6[3];
//    temp5[4] = scale;
    return temp3;
  }
}


const Vector & MultiYieldSurfaceClay::getCommittedStrain (void)
{	
	int ndm = ndmx[matN];

  if (ndm==3)
    return currentStrain.t2Vector(1);
  else {
    static Vector workV(3);//, temp6(6);
		temp6 = currentStrain.t2Vector(1);
    workV[0] = temp6[0];
    workV[1] = temp6[1];
    workV[2] = temp6[3];
    return workV;
  }
}


// NOTE: surfaces[0] is not used 
void MultiYieldSurfaceClay::setUpSurfaces (double * gredu)
{ 
    double residualPress = residualPressx[matN];
    double refPressure = refPressurex[matN];
    double pressDependCoeff =pressDependCoeffx[matN];
    int numOfSurfaces = numOfSurfacesx[matN];
    double frictionAngle = frictionAnglex[matN];
	double cohesion = cohesionx[matN];
    double peakShearStrain = peakShearStrainx[matN];

	double  stress1, stress2, strain1, strain2, size, elasto_plast_modul, plast_modul;
	double pi = 3.14159265358979;
	double refStrain=0;
	double peakShear=0;
	double coneHeight=0;

	if (gredu == 0) {  //automatic generation of surfaces
	  if (frictionAngle > 0) {
	    double sinPhi = sin(frictionAngle * pi/180.);
	    double Mnys = 6.*sinPhi/(3.-sinPhi);
	    residualPress = 3.* cohesion / (sqrt(2.) * Mnys);
		  coneHeight = - (refPressure - residualPress);
		  peakShear = sqrt(2.) * coneHeight * Mnys / 3.; 
		  refStrain = (peakShearStrain * peakShear) 
			            / (refShearModulus * peakShearStrain - peakShear);
		}

	  else if (frictionAngle == 0.) { // cohesion = peakShearStrength
		  peakShear = cohesion;
		  refStrain = (peakShearStrain * peakShear) 
			            / (refShearModulus * peakShearStrain - peakShear);
		  residualPress = 0.;
		}
	
	  double stressInc = peakShear / numOfSurfaces;

#if !_DLL
	  ofstream out("init_surface.out");
	  out << "strain             stress" << endln;
#endif
	  for (int ii=1; ii<=numOfSurfaces; ii++){
        stress1 = ii * stressInc; 
				stress2 = stress1 + stressInc;
        strain1 = stress1 * refStrain / (refShearModulus * refStrain - stress1);
#if !_DLL
		out << strain1 << "      " << stress1 << endln;
#endif

		strain2 = stress2 * refStrain / (refShearModulus * refStrain - stress2);
        if (frictionAngle > 0.) size = 3. * stress1 / sqrt(2.) / coneHeight;
        else if (frictionAngle == 0.) size = 3. * stress1 / sqrt(2.);
 
        elasto_plast_modul = 2.*(stress2 - stress1)/(strain2 - strain1);

        if ( (2.*refShearModulus - elasto_plast_modul) <= 0) 
					plast_modul = UP_LIMIT;
        else 
					plast_modul = (2.*refShearModulus * elasto_plast_modul)/
                        (2.*refShearModulus - elasto_plast_modul);
        if (plast_modul < 0) plast_modul = 0;
        if (plast_modul > UP_LIMIT) plast_modul = UP_LIMIT;
        if (ii==numOfSurfaces) plast_modul = 0;

//		    static Vector temp(6);
			temp.Zero();
        committedSurfaces[ii] = MultiYieldSurface(temp,size,plast_modul);
		}  // ii
#if !_DLL
	  out.close();
#endif
	} 
		else {  //user defined surfaces
			if (frictionAngle > 0) {   // ignore user defined frictionAngle 
			  int ii = 2*(numOfSurfaces-1);
				double tmax = refShearModulus*gredu[ii]*gredu[ii+1];
				double Mnys = -(sqrt(3.) * tmax - 2. * cohesion) / refPressure;
				if (Mnys <= 0) {   // also ignore user defined cohesion
					 cohesion = sqrt(3.)/2 * tmax;
					 frictionAngle = 0.;  
						 coneHeight = 1.;
						 residualPress = 0.;
				  } else {
					 double sinPhi = 3*Mnys /(6+Mnys);
					   if (sinPhi<0. || sinPhi>1.) {
								 opserr <<"\nNDMaterial " <<this->getTag()<<": Invalid friction angle, please modify ref. pressure or G/Gmax curve."<<endln;
					   exit(-1);
					} 
    		   residualPress = 2. * cohesion / Mnys;
			   if (residualPress < 0.01) residualPress = 0.01; 
					coneHeight = - (refPressure - residualPress);
					frictionAngle = asin(sinPhi)*180/pi;
						}
				}  else if (frictionAngle == 0.) {   // ignore user defined cohesion
							int ii = 2*(numOfSurfaces-1);
							double tmax = refShearModulus*gredu[ii]*gredu[ii+1];
					cohesion = sqrt(3.)/2 * tmax;
							coneHeight = 1.;
							residualPress = 0.;
				}

			opserr << "\nNDMaterial " <<this->getTag()<<": Friction angle = "<<frictionAngle
													   <<", Cohesion = "<<cohesion<<"\n"<<endln;

			if (frictionAngle == 0.) pressDependCoeff = 0.; // ignore user defined pressDependCoeff

				for (int i=1; i<numOfSurfaces; i++) {
					int ii = 2*(i-1);
					strain1 = gredu[ii]; 
			  stress1 = refShearModulus*gredu[ii+1]*strain1; 
					strain2 = gredu[ii+2]; 
			  stress2 = refShearModulus*gredu[ii+3]*strain2; 

			  size = sqrt(3.) * stress1 / coneHeight;
			  elasto_plast_modul = 2.*(stress2 - stress1)/(strain2 - strain1);
					if ( (2.*refShearModulus - elasto_plast_modul) <= 0) 
							plast_modul = UP_LIMIT;
			  else 
							plast_modul = (2.*refShearModulus * elasto_plast_modul)/
								(2.*refShearModulus - elasto_plast_modul);
			  if (plast_modul <= 0) {
						opserr << "\nNDMaterial " <<this->getTag()<<": Surface " << i 
							   << " has plastic modulus < 0.\n Please modify G/Gmax curve.\n"<<endln;
				exit(-1);
			  }
			  if (plast_modul > UP_LIMIT) plast_modul = UP_LIMIT;

		//		  static Vector temp(6);
				  temp.Zero();
			  committedSurfaces[i] = MultiYieldSurface(temp,size,plast_modul);

					if (i==(numOfSurfaces-1)) {
						plast_modul = 0;
						size = sqrt(3.) * stress2 / coneHeight;
				committedSurfaces[i+1] = MultiYieldSurface(temp,size,plast_modul);
				}
			}
	  }  

  residualPressx[matN] = residualPress;
  frictionAnglex[matN] = frictionAngle;
  cohesionx[matN] = cohesion;
}



double MultiYieldSurfaceClay::yieldFunc(const T2Vector & stress, 
											 const MultiYieldSurface * surfaces, int surfaceNum)
{
//	static Vector temp(6);
	//temp = stress.deviator() - surfaces[surfaceNum].center();
	temp = stress.deviator();
	temp -= surfaces[surfaceNum].center();

	double sz = surfaces[surfaceNum].size();
	return 3./2.*(temp && temp) - sz * sz;
}


void MultiYieldSurfaceClay::deviatorScaling(T2Vector & stress, const MultiYieldSurface * surfaces, 
																			int surfaceNum, int count)
{
//	changed by guquan Jan 31 2004
	return;
//end
	count++;
	int numOfSurfaces = numOfSurfacesx[matN];
    
	double diff = yieldFunc(stress, surfaces, surfaceNum);

	if ( surfaceNum < numOfSurfaces && diff < 0. ) {
		
		opserr <<"deviatorScaling called, less than 0" << endln;

		double sz = surfaces[surfaceNum].size();
		double deviaSz = sqrt(sz*sz + diff);
//		static Vector devia(6);
		devia = stress.deviator(); 
//		static Vector temp(6);
		temp = devia - surfaces[surfaceNum].center();
		double coeff = (sz-deviaSz) / deviaSz;
		if (coeff < 1.e-13) coeff = 1.e-13;
		devia.addVector(1.0, temp, coeff);
		stress.setData(devia, stress.volume());
		deviatorScaling(stress, surfaces, surfaceNum, count);  // recursive call
	}

	if (surfaceNum==numOfSurfaces && fabs(diff) > LOW_LIMIT) {
		opserr <<"deviatorScaling called,bigger than bound" << endln;
		double sz = surfaces[surfaceNum].size();
		static Vector newDevia(6);
		newDevia.addVector(0.0, stress.deviator(), sz/sqrt(diff+sz*sz));
		stress.setData(newDevia, stress.volume());
	}
}


void MultiYieldSurfaceClay::initSurfaceUpdate()
{
	if (activeSurfaceNum == 0) return; 

	int numOfSurfaces = numOfSurfacesx[matN];

//	static Vector devia(6);
	devia = currentStress.deviator();
	double Ms = sqrt(3./2.*(devia && devia));
	static Vector newCenter(6);

	if (activeSurfaceNum < numOfSurfaces) { // failure surface can't move
		//newCenter = devia * (1. - committedSurfaces[activeSurfaceNum].size() / Ms); 
		newCenter.addVector(0.0, devia, 1.0-committedSurfaces[activeSurfaceNum].size()/Ms);
		committedSurfaces[activeSurfaceNum].setCenter(newCenter);
	}

	for (int i=1; i<activeSurfaceNum; i++) {
	  newCenter = devia * (1. - committedSurfaces[i].size() / Ms);
	  committedSurfaces[i].setCenter(newCenter); 
	}
}


void MultiYieldSurfaceClay::paramScaling(void)
{
	int numOfSurfaces = numOfSurfacesx[matN];
	double frictionAngle = frictionAnglex[matN];
    double residualPress = residualPressx[matN];
    double refPressure = refPressurex[matN];
    double pressDependCoeff =pressDependCoeffx[matN];

	if (frictionAngle == 0.) return;

	double conHeig = - (currentStress.volume() - residualPress);
	double scale = -conHeig / (refPressure-residualPress);
           
	scale = pow(scale, pressDependCoeff); 
	refShearModulus *= scale;
   	refBulkModulus *= scale;

	double plastModul, size;
//	static Vector temp(6);
	temp.Zero();
	for (int i=1; i<=numOfSurfaces; i++) {
	  plastModul = committedSurfaces[i].modulus() * scale;
	  size = committedSurfaces[i].size() * conHeig;
	  committedSurfaces[i] =  MultiYieldSurface(temp,size,plastModul);
	}

}


void MultiYieldSurfaceClay::setTrialStress(T2Vector & stress)
{
//  static Vector devia(6);
  //devia = stress.deviator() + subStrainRate.deviator()*2.*refShearModulus;
  devia = stress.deviator();
  devia.addVector(1.0, subStrainRate.deviator(), 2.*refShearModulus);
  
  trialStress.setData(devia, 0.0);
  //-------- for consistent tangent---------------------------------------
  //-------- dTrialStressdStrain=2*G*I  (4th order tensor)----------------
  

 /* int i,j,k,l,m,n;
  for( i=0;i<3;i++){
	  for( j=i;j<3;j++){
		   m=(-i*i-5*i+7*j-j*j+2*i*j)/2;
		   for(k=0;k<3;k++){
			   for(l=k;l<3;l++){
				   n=(-k*k-5*k+7*l-l*l+2*k*l)/2;
                   dTrialStressdStrain(m,n)=2*refShearModulus*(delta(i,k)*delta(j,l)+delta(i,l)*delta(j,k))/2.0;
			   }
		   }
	  }
  }
*/
  // above is only to check. to avoid this stupid computation, we use:
  dTrialStressdStrain.Zero();
  for(int i=0;i<3;i++){ 
	  dTrialStressdStrain(i,i)=2*refShearModulus;
	  dTrialStressdStrain(i+3,i+3)=refShearModulus;
  }
//  opserr << "step1. setTrialStress, trialStress is" << trialStress.deviator() << endln;
//  opserr << "step1. setTrialStress, dTrialStressdStrain is" << dTrialStressdStrain << endln;

}


int MultiYieldSurfaceClay::setSubStrainRate(void)
{
    int numOfSurfaces = numOfSurfacesx[matN];

	if (activeSurfaceNum==numOfSurfaces) return 1;

	//if (strainRate==T2Vector()) return 0;


// march 4 friday night 2005 for bugs possible cause not convergence 
//	if (strainRate.isZero()) return 0;
// end



	double elast_plast_modulus;
	if (activeSurfaceNum==0) 
	  elast_plast_modulus = 2*refShearModulus;
	else {
	  double plast_modulus = theSurfaces[activeSurfaceNum].modulus();
	  elast_plast_modulus = 2*refShearModulus*plast_modulus 
	    / (2*refShearModulus+plast_modulus);
	}
	static Vector incre(6);
	//incre = strainRate.deviator()*elast_plast_modulus;
	incre.addVector(0.0, strainRate.deviator(),elast_plast_modulus);

	static T2Vector increStress;
	increStress.setData(incre, 0);
	double singleCross = theSurfaces[numOfSurfaces].size() / numOfSurfaces;
	double totalCross = 3.*increStress.octahedralShear() / sqrt(2.);
	int numOfSub = totalCross/singleCross + 1;
	if (numOfSub > numOfSurfaces) numOfSub = numOfSurfaces;
	//incre = strainRate.t2Vector() / numOfSub;
	//---- changed by guquan -------------------------------
	numOfSub=1;
	//---- end change --------------------------------------
	incre = strainRate.t2Vector();
	incre /= numOfSub;
	subStrainRate.setData(incre);

	return numOfSub;
}


void
MultiYieldSurfaceClay::getContactStress(T2Vector &contactStress)
{
	static Vector center(6);
	center = theSurfaces[activeSurfaceNum].center(); 
//	static Vector devia(6);
	static Vector tempStress(6);
	
	static Vector dKdStrain(6);
	static Matrix tempTangent(6,6);
	
	//devia = trialStress.deviator() - center;
	devia = trialStress.deviator();
	devia -= center;
    
	// save for use later
	tempStress=devia;
	
	double Ms = sqrt(3./2.*(devia && devia));
	//devia = devia * theSurfaces[activeSurfaceNum].size() / Ms + center;
	devia *= theSurfaces[activeSurfaceNum].size() / Ms;
	devia += center;

	contactStress.setData(devia,0.0); 

	
	
//  -------- compute for consistent Tangent -------------------

	doubledotProduct(dKdStrain,tempStress,dTrialStressdStrain);
	dKdStrain *= 3/(2.*Ms);
  // output for debug
//	opserr << "step2.1 getContactStress, K is" << Ms << endln;
//	opserr << "step2.1 getContactStress, dKdStrain is" << dKdStrain << endln;
	
	tempTangent.Zero();
	tensorProduct(tempTangent,tempStress,dKdStrain);
	dContactStressdStrain.addMatrix(0.0,dTrialStressdStrain,theSurfaces[activeSurfaceNum].size() / Ms);
	dContactStressdStrain.addMatrix(1.0,tempTangent,-1.0*theSurfaces[activeSurfaceNum].size() / Ms/Ms);
  
	// output for debug
//	opserr << "step2.2 getContactStress, contactStress is" << contactStress.deviator() << endln;
//	opserr << "step2.2 getContactStress, dContactStressdStrain is" << dContactStressdStrain << endln;

}


int MultiYieldSurfaceClay::isLoadReversal(void)
{
  if(activeSurfaceNum == 0) return 0;

  static Vector surfaceNormal(6);
  getSurfaceNormal(currentStress, surfaceNormal);
 
  //(((trialStress.deviator() - currentStress.deviator()) && surfaceNormal) < 0) 
  // return 1;
  static Vector a(6);
  a = trialStress.deviator();
  a-= currentStress.deviator();
  if((a && surfaceNormal) < 0) 
    return 1;

  return 0;   
}


void
MultiYieldSurfaceClay::getSurfaceNormal(const T2Vector & stress, Vector &surfaceNormal)
{
  //Q = stress.deviator() - theSurfaces[activeSurfaceNum].center();
  // return Q / sqrt(Q && Q);

  static Vector tempStress(6),tempProduct(6);
  static Matrix tempTangent(6,6);


  surfaceNormal = stress.deviator();
  surfaceNormal -= theSurfaces[activeSurfaceNum].center();
  // save for use later
  tempStress.addVector(0.0,surfaceNormal,1.0);
  double k=1/sqrt(surfaceNormal && surfaceNormal);

  surfaceNormal /= sqrt(surfaceNormal && surfaceNormal);



  //------------- Consistent Tangent-------------------------------------
  tempProduct.Zero();
  dSurfaceNormaldStrain.Zero();
  doubledotProduct(tempProduct,tempStress,dContactStressdStrain);
  tensorProduct(dSurfaceNormaldStrain,tempStress,tempProduct);
  dSurfaceNormaldStrain *= -1.0*k*k*k;
  dSurfaceNormaldStrain.addMatrix(1.0,dContactStressdStrain,k);


  	// output for debug
//	opserr << "step3 getSurfaceNormal, surfaceNormal is" << surfaceNormal << endln;
//	opserr << "step3 getSurfaceNormal, dSurfaceNormaldStrain is" << dSurfaceNormaldStrain<< endln;

  
}


double MultiYieldSurfaceClay::getLoadingFunc(const T2Vector & contactStress, 
									 const Vector & surfaceNormal, int crossedSurface)
{
  double loadingFunc;
  double temp1 = 2. * refShearModulus ;
  double temp2 = theSurfaces[activeSurfaceNum].modulus();

  Vector tempStress(6);
  Matrix tempTangent(6,6);

  //for crossing first surface
  double temp_1 = temp1 + temp2;
  //loadingFunc = (surfaceNormal && (trialStress.deviator()-contactStress.deviator()))/temp_1;
//  static Vector temp(6);
  temp =trialStress.deviator();
  temp -= contactStress.deviator();
  loadingFunc = (surfaceNormal && temp)/temp_1;
   //for crossing more than one surface
  if(crossedSurface) {
    double temp3 = theSurfaces[activeSurfaceNum-1].modulus();
    loadingFunc *= (temp3 - temp2)/temp3;
  }

  // ------------ consistent tangent -----------------------------
    tempStress.Zero();
	dXdStrain.Zero();
	tempTangent.addMatrix(0.0,dTrialStressdStrain,1.0);
	tempTangent.addMatrix(1.0,dContactStressdStrain,-1.0);
	doubledotProduct(tempStress,surfaceNormal,tempTangent);
	
	doubledotProduct(dXdStrain,temp,dSurfaceNormaldStrain);
	dXdStrain.addVector(1.0,tempStress,1.0);
	dXdStrain /=temp_1;

   if(crossedSurface) {
    double temp3 = theSurfaces[activeSurfaceNum-1].modulus();
    dXdStrain *= (temp3 - temp2)/temp3;
  }

     	// output for debug
//	opserr << "step4 getLoadFunc, X is" << loadingFunc << endln;
//	opserr << "step4 getLoadingFunc, dXdStrain is" << dXdStrain<< endln;

 //---------------------------------------------------------------
  return loadingFunc;
}


void MultiYieldSurfaceClay::stressCorrection(int crossedSurface)
{
	static T2Vector contactStress;
	this->getContactStress(contactStress);
	static Vector surfaceNormal(6);
	this->getSurfaceNormal(contactStress, surfaceNormal);
	double loadingFunc = getLoadingFunc(contactStress, surfaceNormal, crossedSurface);
//	static Vector devia(6);

	Matrix tempTangent(6,6);

	//devia = trialStress.deviator() - surfaceNormal * 2 * refShearModulus * loadingFunc;
	devia.addVector(0.0, surfaceNormal, -2*refShearModulus*loadingFunc);
	devia += trialStress.deviator();

	//----------- add for consistent tangen ---------------------------
	tensorProduct(tempTangent,surfaceNormal,dXdStrain);
	dTrialStressdStrain.addMatrix(1.0,tempTangent,-2.0*refShearModulus);
	dTrialStressdStrain.addMatrix(1.0,dSurfaceNormaldStrain,-2.*refShearModulus*loadingFunc);


    trialStress.setData(devia, 0.0);



	// output strainRate for debug
//	opserr << "step5.stressCorrection, trialStress is," << trialStress.deviator() << endln;
//   opserr << "step5.stressCorrection, dTrialStressdStrain is," << dTrialStressdStrain << endln;



	deviatorScaling(trialStress, theSurfaces, activeSurfaceNum);

	if (isCrossingNextSurface()) {
		activeSurfaceNum++;
		stressCorrection(1);  //recursive call
	}
}


void MultiYieldSurfaceClay::updateActiveSurface(void)
{
  int numOfSurfaces = numOfSurfacesx[matN];

  if (activeSurfaceNum == numOfSurfaces) return;

	double A, B, C, X;
	static T2Vector direction;
	static Vector t1(6);
	static Vector t2(6);
//	static Vector temp(6);
	static Vector center(6);
	center = theSurfaces[activeSurfaceNum].center();
	double size = theSurfaces[activeSurfaceNum].size();
	static Vector outcenter(6);
	outcenter= theSurfaces[activeSurfaceNum+1].center();
	double outsize = theSurfaces[activeSurfaceNum+1].size();


	//t1 = trialStress.deviator() - center;
	//t2 = center - outcenter;
	t1 = trialStress.deviator();
	t1 -= center;
	t2 = center;
	t2 -= outcenter;

	A = t1 && t1;
	B = 2. * (t1 && t2);
	C = (t2 && t2) - 2./3.* outsize * outsize;
	X = secondOrderEqn(A,B,C,0);
	if ( fabs(X-1.) < LOW_LIMIT ) X = 1.;
	if (X < 1.){
	  opserr << "FATAL:MultiYieldSurfaceClay::updateActiveSurface(): error in Direction of surface motion." 
	       << endln; 
	  exit(-1);
	}

	//temp = (t1 * X + center) * (1. - size / outsize) - (center - outcenter * size / outsize);
	temp = center;
	temp.addVector(1.0, t1, X);
	temp *= (1.0 - size/outsize);
	t2 = center;
	t2.addVector(1.0, outcenter, -size/outsize);
	temp -= t2;

	direction.setData(temp);

	if (direction.deviatorLength() < LOW_LIMIT) return;

	temp = direction.deviator();  
	A = temp && temp;
	B = - 2 * (t1 && temp);
	if (fabs(B) < LOW_LIMIT) B = 0.; 
	C = (t1 && t1) - 2./3.* size * size;
	if ( fabs(C) < LOW_LIMIT || fabs(C)/(t1 && t1) < LOW_LIMIT ) return;

// guquan 2005 Apr. 21
//	opserr<<"Low_limit is " << LOW_LIMIT << endln;
	if (fabs(C)<1.0e-12)   C=0;
// end guquan

	if (B > 0. || C < 0.) {
	  opserr << "FATAL:MultiYieldSurfaceClay::updateActiveSurface(): error in surface motion.\n" 
	       << "A= " <<A <<" B= " <<B <<" C= "<<C <<" (t1&&t1)= "<<(t1&&t1) <<endln; 
	  exit(-1);
	}
	X = secondOrderEqn(A,B,C,1);  

	//center += temp * X;
	center.addVector(1.0, temp, X);
	theSurfaces[activeSurfaceNum].setCenter(center);
}      


void MultiYieldSurfaceClay::updateInnerSurface(void)
{
	if (activeSurfaceNum <= 1) return;

//	static Vector devia(6);
	devia = currentStress.deviator();
	static Vector center(6);
	center = theSurfaces[activeSurfaceNum].center();
	double size = theSurfaces[activeSurfaceNum].size();
	static Vector newcenter(6);

	for (int i=1; i<activeSurfaceNum; i++) {
		//newcenter = devia - (devia - center) * theSurfaces[i].size() / size;
		newcenter = center; 
		newcenter -= devia;
		newcenter *= theSurfaces[i].size()/size;
		newcenter += devia;

		theSurfaces[i].setCenter(newcenter);
	}
}


int MultiYieldSurfaceClay:: isCrossingNextSurface(void)
{
  int numOfSurfaces = numOfSurfacesx[matN];
  if (activeSurfaceNum == numOfSurfaces) return 0;  

  if(yieldFunc(trialStress, theSurfaces, activeSurfaceNum+1) > 0) return 1;
  
  return 0;
}
 



///////////////////////// add sensitivity ///////////////////////////////////
int  
MultiYieldSurfaceClay:: setParameter(const char **argv, int argc, Parameter &param)
{	if (argc < 1)
		return -1;

	if (strcmp(argv[0],"G") == 0) {
		return param.addObject(1, this);
	}

	if (strcmp(argv[0],"cohesion") == 0) {
		return param.addObject(2, this);
	}

	if (strcmp(argv[0],"K") == 0) {
		return param.addObject(3, this);
	}
	else
		opserr << "WARNING: Could not set parameter in MultiYieldSurfaceClay. " << endln;
                
	return -1;
}
	
int 
MultiYieldSurfaceClay::updateParameter(int passedParameterID, 
												 Information &info){
//     opserr << "updatePatameter is called " << endln;
     // exit(-1); // make sure it is not called
	// I still can not understand when this is called and for what purpose.

// --------switch 6  used! --------------------------------
	switch (passedParameterID) {
	case -1:
		return -1;
	case 1:
		this->refShearModulus= info.theDouble; // 
		break;
	case 2:
		cohesionx[matN]=info.theDouble;
//		this->peakShearStrainx[matN]  = info.theDouble;
		break;
	case 3:
		this->refBulkModulus  = info.theDouble;
		break;

	default:
		return -1;
	}
    

	this->setUpSurfaces(0);
	return 0;
	}

int MultiYieldSurfaceClay::activateParameter(int passedParameterID){
	parameterID = passedParameterID;
	return 0;						
	}




void MultiYieldSurfaceClay::setUpSurfacesSensitivity(int passedGradNumber )
{ 

// --- to be consistent with MHS new framework ---- 2009
// ==== Not change here but when call it. 
//	passedGradNumber +=1;

//  opserr <<"for checking  setupSurfaceSensi version 1.27"<<endln;
// do some checking job	
	if( surfacesSensitivityMark==0) {
	opserr << "surfacesSensitivityMark not exist !!!" << endln;
	exit (-1);
	}
	
	for (int i=1; i<passedGradNumber; i++) {
		if (this-> isSurfacesSensitivitySetUp(i) != 1) {
				opserr << "previous Grad's surfacesSensitivityMark not set up !!!" << endln;
				exit (-1);
		}
	}

    int numOfSurfaces = numOfSurfacesx[matN];
    double frictionAngle = frictionAnglex[matN];
	
	double cohesion = cohesionx[matN];
    double peakShearStrain = peakShearStrainx[matN];

	double  stress1, stress2, strain1, strain2, size, elasto_plast_modul, plast_modul;
	double pi = 3.14159265358979;
	double refStrain, peakShear;


	double dPeakShear;
	double dRefShearModulus;
	double dRefStrain;
    double stressInc;
    double dStressInc;
	double dStress1;
	double dStress2;
	double dStrain1;
	double dStrain2;
	double dSize;
	double dElasto_plast_modul;
	double dPlast_modul;




	if (frictionAngle != 0) {
	    opserr << "Fatal: can not deal with frictionAngle != 0 case now" << endln;
		exit(-1);
	  }
      
// sensitivity variables: sita=G: in [ii,0]    peakShear = cohesion : in [ii,1]

//for(int grad_Number=0; grad_Number<myNumGrads;grad_Number++)	  

int grad_Number=passedGradNumber-1;

//{	  
 dPeakShear=0.0;
 dRefShearModulus=0.0;




 // -----------switch 7 here --------------------------
 // place 4 for add if grad_Number ! =0 
	if (parameterID !=0) {   // doing nothing
		if(parameterID==1) dRefShearModulus=1.0;
		if(parameterID==2) dPeakShear=1.0;
//		if(parameterID==3) ; //for Bulkmodulus, doing nothing
	}

     if (frictionAngle == 0.) { // cohesion = peakShearStrength
		  peakShear = cohesion;
		  refStrain = (peakShearStrain * peakShear) 
			            / (refShearModulus * peakShearStrain - peakShear);

		  dRefStrain=1.0/pow(refShearModulus*peakShearStrain-peakShear,2)*
			  ((dPeakShear*peakShearStrain)*(refShearModulus * peakShearStrain - peakShear)
			  -peakShearStrain * peakShear*(dRefShearModulus*peakShearStrain-dPeakShear));

		}
	
	  stressInc = peakShear / numOfSurfaces;
	  dStressInc = dPeakShear / numOfSurfaces;

	  for (int ii=1; ii<=numOfSurfaces; ii++){   // Note: surface 0 is never used!!
        stress1 = ii * stressInc;
		dStress1=ii*dStressInc;

		stress2 = stress1 + stressInc;
		dStress2=dStress1+dStressInc;
        
		strain1 = stress1 * refStrain / (refShearModulus * refStrain - stress1);
		dStrain1=1.0/pow((refShearModulus * refStrain - stress1),2)*(
			(dStress1 * refStrain+stress1*dRefStrain)*(refShearModulus * refStrain - stress1)
			-(stress1 * refStrain)*(dRefShearModulus * refStrain+refShearModulus * dRefStrain - dStress1));
        
		strain2 = stress2 * refStrain / (refShearModulus * refStrain - stress2);
		dStrain2=1.0/pow((refShearModulus * refStrain - stress2),2)*(
			(dStress2 * refStrain+stress2*dRefStrain)*(refShearModulus * refStrain - stress2)
			-(stress2 * refStrain)*(dRefShearModulus * refStrain+refShearModulus * dRefStrain - dStress2));
        

		if (frictionAngle > 0.) {
				opserr << "Fatal: can not deal with frictionAngle != 0 case now" << endln;
				exit(-1);}
        else if (frictionAngle == 0.) {
				size = 3. * stress1 / sqrt(2.);
				dSize=3.*dStress1/sqrt(2.);}
 
        elasto_plast_modul = 2.*(stress2 - stress1)/(strain2 - strain1);
		dElasto_plast_modul=2./pow((strain2 - strain1),2)*((dStress2 - dStress1)*(strain2 - strain1)
			-(stress2 - stress1)*(dStrain2 - dStrain1));

			if ( (2.*refShearModulus - elasto_plast_modul) <= 0)  {
			 opserr << "Fatal: can not deal with plast_model< 0 case now" << endln;
		     exit(-1);
			}
			else {
			plast_modul = (2.*refShearModulus * elasto_plast_modul)/
                        (2.*refShearModulus - elasto_plast_modul);
			dPlast_modul = 2./pow((2.*refShearModulus - elasto_plast_modul),2)*(
				(dRefShearModulus * elasto_plast_modul+refShearModulus * dElasto_plast_modul)*
				(2.*refShearModulus - elasto_plast_modul)-(refShearModulus * elasto_plast_modul)*
				(2.*dRefShearModulus - dElasto_plast_modul));
			}

        if ((plast_modul < 0)&&(ii !=numOfSurfaces)) {
			 opserr << "Fatal: can not deal with plast_model<0 case now" << endln;
		     exit(-1);
		}
        if (plast_modul > UP_LIMIT){ 
			opserr << "Fatal: can not deal with plast_modul > UP_LIMIT" << endln;
			exit(-1);
		}
        if (ii==numOfSurfaces) {
			plast_modul = 0;
			dPlast_modul = 0;
		}

		if (dCommittedMultiSurfaceSize == 0){ //create and initialize array. 
		dMultiSurfaceCenter = new double [6*(numOfSurfaces+1)*myNumGrads];
		dCommittedMultiSurfaceSize=new double [(numOfSurfaces+1)*myNumGrads];
		dCommittedMultiSurfacePlastModul=new double [(numOfSurfaces+1)*myNumGrads];
		dCommittedMultiSurfaceCenter = new double [6*(numOfSurfaces+1)*myNumGrads];
		
			for (int i=0; i<numOfSurfaces+1; i++){
				for(int j=0;j<myNumGrads;j++){
					for(int k=0;k<6;k++) dMultiSurfaceCenter[k+i*6+j*6*(numOfSurfaces+1)]=0.0;
					dCommittedMultiSurfaceSize[i+j*(numOfSurfaces+1)]=0.0;
					dCommittedMultiSurfacePlastModul[i+j*(numOfSurfaces+1)]=0.0;
					for(int k=0;k<6;k++) dCommittedMultiSurfaceCenter[k+i*6+j*6*(numOfSurfaces+1)]=0.0;
				}
			} //for
		} //if




	// Note: surface 0 is never used !!
		dCommittedMultiSurfaceSize[ii+grad_Number*(numOfSurfaces+1)]=dSize;
		dCommittedMultiSurfacePlastModul[ii+grad_Number*(numOfSurfaces+1)]=dPlast_modul;
			if ((ii==1) &&(grad_Number==0))
			{
				  opserr.setPrecision(16);

//				  opserr << "step1. setUpSurfaceSensitivity, size is:"<<size<< endln;
//				  opserr << "step1. setUpSurfaceSensitivity, dsize is:"<<dSize<< endln;
//				  opserr << "step1. setUpSurfaceSensitivity, plastModul  is:"<<plast_modul<< endln;
//				  opserr << "step1. setUpSurfaceSensitivity, dPlastModul is:"<<dPlast_modul<< endln;

			}

 		}  // ii
	
//}//grad_Number
//*/
	  this-> setSurfacesSensitivitySetUpMark(passedGradNumber);

//		opserr<<"*/*/*/*/" <<endln;
} 


void MultiYieldSurfaceClay::setTrialStressSensitivity(T2Vector & stress,T2Vector & dStress)
{
  static 
	Vector /* devia(6),*/  dTempStress(6);
  //devia = stress.deviator() + subStrainRate.deviator()*2.*refShearModulus;
  devia = stress.deviator();
  devia.addVector(1.0, subStrainRate.deviator(), 2.*refShearModulus);
  trialStress.setData(devia, 0.0);
//  opserr << "step3. setTrialStressSensitivity, trialStress is:"<< endln;
//  opserr << trialStress.deviator()<< endln;

// ---------sensitivity part -------------------------------------------
// double dPeakShear=0.0;
 double dRefShearModulus=0.0;

 //-----------------------switch 1 here ------------------------
 // place 1 need to add if parameterID !=0  bug here
 if (parameterID !=0) {  //doing nothing
  if(parameterID==1) dRefShearModulus=1.0;
//  if(parameterID==2) dPeakShear=1.0;  doing nothing
//  if(parameterID==3) ; //for Bulkmodulus, doing nothing
 
 }   // end change

//	dTrialStress=dStress+2*dRefShearModulus*subStrainRate.deviator()+2*G*dStrain.deviator();
   dTempStress=dStress.deviator();
   dTempStress.addVector(1.0, subStrainRate.deviator(), 2.*dRefShearModulus);
   dTempStress.addVector(1.0, dSubStrainRate.deviator(), 2.*refShearModulus);		
   dTrialStress.setData(dTempStress,0.0);	   
//  opserr << "step3. setTrialStressSensitivity, dTrialStress is:"<< endln;
//  opserr << dTrialStress.deviator()<< endln;

}



void MultiYieldSurfaceClay::updateInnerSurfaceSensitivity(void)
{
	if (activeSurfaceNum <= 1) return;

	int numOfSurfaces=numOfSurfacesx[matN];
	
//	static 		Vector devia(6);
	devia = currentStress.deviator();
	static 
	     Vector center(6);
	center = theSurfaces[activeSurfaceNum].center();
	double size = theSurfaces[activeSurfaceNum].size();
	static 
		Vector newcenter(6);


	static	
		T2Vector dStress;   
	static	
		Vector dTempStress(6),dCenter(6),tempStress(6);
	double dSize;
	
	for (int ii=1; ii<activeSurfaceNum; ii++) {
		//newcenter = devia - (devia - center) * theSurfaces[ii].size() / size;
		devia = currentStress.deviator();
		newcenter = center; 
		newcenter -= devia;
		newcenter *= theSurfaces[ii].size()/size;
		newcenter += devia;

		theSurfaces[ii].setCenter(newcenter);

//        opserr << "step2. updateInnerSurfaceSensitivity, theSurfaces "<<ii<<" is:"<< endln;
//        opserr << newcenter<< endln;	
		
		
		//----------------- sensitivity part---------------------------------------
	
		dSize=dCommittedMultiSurfaceSize[ii+(gradNumber-1)*(numOfSurfaces+1)];
		for (int N=0;N<6;N++) dCenter(N)=dMultiSurfaceCenter[N+activeSurfaceNum*6+(gradNumber-1)*6*(numOfSurfaces+1)];
		dTempStress=dCurrentStress.deviator();
		dTempStress.addVector(theSurfaces[ii].size(),dCenter,-1.0*theSurfaces[ii].size());
		
		tempStress=devia;
		tempStress.addVector(dSize,center,-1.0*dSize);
		
		tempStress.addVector(size,dTempStress,size);

		devia.addVector(1.0,center,-1.0); 
		
		devia *= theSurfaces[ii].size();
		devia *= dCommittedMultiSurfaceSize[activeSurfaceNum+(gradNumber-1)*(numOfSurfaces+1)];

		tempStress.addVector(1./pow(size,2),devia,-1./pow(size,2)); ///????

		dTempStress=dCurrentStress.deviator();
		dTempStress.addVector(1.0,tempStress,-1.0);

		for (int N=0;N<6;N++) dMultiSurfaceCenter[N+ii*6+(gradNumber-1)*6*(numOfSurfaces+1)]=dTempStress(N);

//       opserr << "step2. updateInnerSurfaceSensitivity, Sensitiviry of theSurfaces "<<ii<<" is:"<< endln;
//       opserr << dTempStress<< endln;

	
	} //for ii
	
}



int MultiYieldSurfaceClay::setSubStrainRateSensitivity(void)
{
    int numOfSurfaces = numOfSurfacesx[matN];

	if (activeSurfaceNum==numOfSurfaces) return 1;



	double elast_plast_modulus;
	if (activeSurfaceNum==0) 
	  elast_plast_modulus = 2*refShearModulus;
	else {
	  double plast_modulus = theSurfaces[activeSurfaceNum].modulus();
	  elast_plast_modulus = 2*refShearModulus*plast_modulus 
	    / (2*refShearModulus+plast_modulus);
	}
	static 
		Vector incre(6);
	//incre = strainRate.deviator()*elast_plast_modulus;
	incre.addVector(0.0, strainRate.deviator(),elast_plast_modulus);

	static
		T2Vector increStress;
	increStress.setData(incre, 0);
	double singleCross = theSurfaces[numOfSurfaces].size() / numOfSurfaces;
	double totalCross = 3.*increStress.octahedralShear() / sqrt(2.);
	int numOfSub = totalCross/singleCross + 1;
	if (numOfSub > numOfSurfaces) numOfSub = numOfSurfaces;
	//incre = strainRate.t2Vector() / numOfSub;
	//---- changed by guquan --------------------
	numOfSub=1;

	incre = strainRate.t2Vector();
	incre /= numOfSub;
	subStrainRate.setData(incre);

// -------------- sensitivity part ---------------------------------


	static	
		Vector dTempStrain(6);
    dTempStrain=dStrainRate.deviator();
	dTempStrain/=numOfSub;
	dSubStrainRate.setData(dTempStrain,0.0);


	return numOfSub;
}




void
MultiYieldSurfaceClay::getContactStressSensitivity(T2Vector &contactStress)
{
	static 
		Vector center(6);
	center = theSurfaces[activeSurfaceNum].center(); 
	//static 		Vector devia(6);
	//devia = trialStress.deviator() - center;
	devia = trialStress.deviator();
	devia -= center;

	double Ms = sqrt(3./2.*(devia && devia));
	//devia = devia * theSurfaces[activeSurfaceNum].size() / Ms + center;
	devia *= theSurfaces[activeSurfaceNum].size() / Ms;
	devia += center;

	contactStress.setData(devia,0.0); 
//    opserr << "step41. getContactStressSensitivity, contactStress is:"<< endln;
//    opserr << contactStress.deviator()<< endln;


// -------------- sensitivity part -------------------------------------


	
	static 
		Vector dTempStress(6),dCenter(6),tempStress(6);
	double dMs;
    int numOfSurfaces = numOfSurfacesx[matN];

	devia = trialStress.deviator();
	devia -= center;

	dTempStress=dTrialStress.deviator();
	for (int N=0;N<6;N++) dCenter(N)=dMultiSurfaceCenter[N+activeSurfaceNum*6+(gradNumber-1)*6*(numOfSurfaces+1)];
	dTempStress.addVector(1.0,dCenter,-1.0);
	
	dMs=3./(2.*Ms)*(dTempStress&&devia);
	
	tempStress=devia*Ms*dCommittedMultiSurfaceSize[activeSurfaceNum+(gradNumber-1)*(numOfSurfaces+1)];
	devia *=theSurfaces[activeSurfaceNum].size()*dMs;
	tempStress.addVector(1.0,devia,-1.0);
	
	dTempStress *= Ms*theSurfaces[activeSurfaceNum].size();
	
	tempStress.addVector(1.0,dTempStress,1.0);
	tempStress /=pow(Ms,2);
	tempStress +=dCenter;


	dContactStress.setData(tempStress,0.0);
//    opserr << "step41. getContactStressSensitivity, dContactStress is:"<< endln;
//    opserr << dContactStress.deviator()<< endln;

}




void
MultiYieldSurfaceClay::getSurfaceNormalSensitivity(const T2Vector & stress,
														const T2Vector & dStress1,
														Vector &surfaceNormal,
														Vector &dSurfaceNormal )
{
  //Q = stress.deviator() - theSurfaces[activeSurfaceNum].center();
  // return Q / sqrt(Q && Q);
	static 
	   Vector dCenter(6),tempStress(6);
	double tempNormal,temp;	
    int numOfSurfaces = numOfSurfacesx[matN];

	surfaceNormal = stress.deviator();
	surfaceNormal -= theSurfaces[activeSurfaceNum].center();

	tempStress = surfaceNormal;
	tempNormal= sqrt(surfaceNormal && surfaceNormal);

	surfaceNormal /= tempNormal;
// 	opserr << "step42. getSurfaceNormalSensitivity, surfaceNormal is:"<< endln;
//    opserr <<  surfaceNormal<< endln;


//------------------------------ sensitivity part ------------------------------
	for (int N=0;N<6;N++) dCenter(N)=dMultiSurfaceCenter[N+activeSurfaceNum*6+(gradNumber-1)*6*(numOfSurfaces+1)];
	
	dSurfaceNormal=dStress1.deviator();
	
	dSurfaceNormal.addVector(1.0,dCenter,-1.0);
	
	temp=(dSurfaceNormal && tempStress)/tempNormal;

	dSurfaceNormal.addVector(tempNormal,tempStress,-1.0*temp);
	dSurfaceNormal /= pow(tempNormal,2);

//	opserr << "step42. getSurfaceNormalSensitivity, dSurfaceNormal is:"<< endln;
//    opserr <<  dSurfaceNormal<< endln;

}








double MultiYieldSurfaceClay::getLoadingFuncSensitivity(const T2Vector & contactStress, 
									 const Vector & surfaceNormal,const Vector & dSurfaceNormal,
									 int crossedSurface)
{
  double loadingFunc;
  double temp1 = 2. * refShearModulus ;
  double temp2 = theSurfaces[activeSurfaceNum].modulus();

  //for crossing first surface
  double temp_1 = temp1 + temp2;
  //loadingFunc = (surfaceNormal && (trialStress.deviator()-contactStress.deviator()))/temp_1;
  //static 
//	  Vector temp(6);
  temp =trialStress.deviator();
  temp -= contactStress.deviator();
  loadingFunc = (surfaceNormal && temp)/temp_1;
   //for crossing more than one surface
  if(crossedSurface) {
    double temp3 = theSurfaces[activeSurfaceNum-1].modulus();
    loadingFunc *= (temp3 - temp2)/temp3;
  }

//    opserr << "step43. getLoadingFunc, loadingFunc is:"<< endln;
//    opserr <<  loadingFunc<< endln;

//------------ sensitivity part --------------------------------

 double dPlastModulus,dRefShearModulus;
 double temp4,temp5;
 static 
	 Vector dTempStress(6);
 int numOfSurfaces = numOfSurfacesx[matN];
 
 dPlastModulus=dCommittedMultiSurfacePlastModul[activeSurfaceNum+(gradNumber-1)*(numOfSurfaces+1)];
 dRefShearModulus=0.0;
 
 // -------switch 2 here --------------------
 //place 2 need add if (parameterID 1=0)
 if (parameterID !=0){
	if (parameterID==1) dRefShearModulus=1.0;  
//	if (parameterID==2); //  dPeakShear=1.0;  doing nothing, dPeakShear never used
 }//end change

 temp4=(dSurfaceNormal && temp);
 dTempStress=dTrialStress.deviator();
 dTempStress.addVector(1.0,dContactStress.deviator(),-1.0);
 
 temp4=temp4+(surfaceNormal && dTempStress);
 
 temp5=(surfaceNormal && temp);
 temp5 *=(dPlastModulus+2.*dRefShearModulus);
 dLoadingFunc = (temp4*temp_1-temp5)/pow(temp_1,2);

 if(crossedSurface) {
    double temp3 = theSurfaces[activeSurfaceNum-1].modulus();
	double dTemp3=dCommittedMultiSurfacePlastModul[(activeSurfaceNum-1)+(gradNumber-1)*(numOfSurfaces+1)];
	dLoadingFunc *=(temp3-temp2)/temp3;
	dLoadingFunc +=loadingFunc*(-1.0*dPlastModulus*temp3+temp2*dTemp3)/pow(temp3,2);

  }

//    opserr << "step43. getLoadingFunc, dLoadingFunc is:"<< endln;
//    opserr <<  dLoadingFunc<< endln;

  return loadingFunc;
}




void MultiYieldSurfaceClay::updateActiveSurfaceSensitivity(void)
{
  int numOfSurfaces = numOfSurfacesx[matN];

  if (activeSurfaceNum == numOfSurfaces) return;

	double A, B, C, X;
	static 
		T2Vector direction;
	static 
		Vector t1(6);
	static 
		Vector t2(6);
	//static
//		Vector temp(6);
	static
		Vector center(6);
	center = theSurfaces[activeSurfaceNum].center();
	double size = theSurfaces[activeSurfaceNum].size();
	static
		Vector outcenter(6);
	outcenter= theSurfaces[activeSurfaceNum+1].center();
	double outsize = theSurfaces[activeSurfaceNum+1].size();

	int mark=0.;
    //----------------- step1.-----------------------------
	
	//t1 = trialStress.deviator() - center;
	//t2 = center - outcenter;
	
	t1 = trialStress.deviator();
	t1 -= center;
	t2 = center;
	t2 -= outcenter;

	A = t1 && t1;
	B = 2. * (t1 && t2);
	C = (t2 && t2) - 2./3.* outsize * outsize;
	X = secondOrderEqn(A,B,C,0);
	if ( fabs(X-1.) < LOW_LIMIT ) {X = 1.; mark=1;}
	if (X < 1.){
	  opserr << "FATAL:MultiYieldSurfaceClay::updateActiveSurface(): error in Direction of surface motion." 
	       << endln; 
	  exit(-1);
	}

   // ------- sensitivity for step 1 .------------------------

	static 
		Vector dTempStress(6),dCenter(6),tempStress(6),dOutCenter(6);

	Vector dTempStress1(6),dMu(6);
	double	dA,dB,dC,dX,dSize,dOutSize;
	dSize   = dCommittedMultiSurfaceSize[activeSurfaceNum+(gradNumber-1)*(numOfSurfaces+1)];
	dOutSize= dCommittedMultiSurfaceSize[activeSurfaceNum+1+(gradNumber-1)*(numOfSurfaces+1)];
	static
		T2Vector dDirection;

	for (int N=0;N<6;N++) dCenter(N)=dMultiSurfaceCenter[N+activeSurfaceNum*6+(gradNumber-1)*6*(numOfSurfaces+1)];
	for (int N=0;N<6;N++)     dOutCenter(N)=dMultiSurfaceCenter[N+(activeSurfaceNum+1)*6+(gradNumber-1)*6*(numOfSurfaces+1)];
	
	dTempStress=dTrialStress.deviator();
	dTempStress.addVector(1.0,dCenter,-1.0);

	dOutCenter.addVector(-1.0,dCenter,1.0);
	
	dA=2.*(dTempStress && t1);
	dB=2.*(dOutCenter && t1)+2.*(t2 &&dTempStress);
	dC=2.*(dOutCenter && t2)-4./3.*outsize*dOutSize;
	
	dX=(-dA*X*X-dB*X-dC)/(2.*A*X+B);
//	if ( fabs(X-1.) < LOW_LIMIT ) dX = 0.; 
	if (mark==1) {dX=0.;mark=0;}


	//---------step2. ----------------------------------------
	//temp = (t1 * X + center) * (1. - size / outsize) - (center - outcenter * size / outsize);
	temp = center;
	temp.addVector(1.0, t1, X);
    
//	opserr << "step51. updateActiveSurface. TaoStar is:"<< endln;
//    opserr <<  temp << endln;

	tempStress=temp;  // save for sensitivity purpose 

	temp *= (1.0 - size/outsize);
	t2 = center;
	t2.addVector(1.0, outcenter, -size/outsize);
	temp -= t2;

	direction.setData(temp);
	// added by guquan and maybe cause problem.
	temp=direction.deviator();

	if (direction.deviatorLength() < LOW_LIMIT) return;


	//-------------- sensitivity part----------------------------


	//---- dTempStress1=dCenter+dX*t1+X*dTempStress
	dTempStress1=dCenter;
	dTempStress1.addVector(1.0,t1,dX);
	dTempStress1.addVector(1.0,dTempStress,X);

//	opserr << "step51. updateActiveSurface. dTaoStar is:"<< endln;
//    opserr <<  dTempStress1 << endln;


	tempStress.addVector(1.0,outcenter,-1.0);
	double tt=(dSize*outsize-size*dOutSize)/(pow(outsize,2));

	dMu=dTempStress1;
	dMu.addVector(1.0,dCenter,-1.0);
	dMu.addVector(1.0,tempStress,-1.0*tt);

	for (int N=0;N<6;N++) dOutCenter(N)=dMultiSurfaceCenter[N+(activeSurfaceNum+1)*6+(gradNumber-1)*6*(numOfSurfaces+1)];
	dTempStress1.addVector(1.0,dOutCenter,-1.0);

	dMu.addVector(1.0,dTempStress1,-1.0*size/outsize);

	dDirection.setData(dMu);
	dMu=dDirection.deviator();

/* 	opserr << "step52. updateActiveSurface. Mu is:"<< endln;
    opserr <<  direction.deviator() << endln;
    opserr << "step52. updateActiveSurface. dMu is:"<< endln;
    opserr <<  dDirection.deviator() << endln;
*/

	//----------- step 3. -----------------------------------

	temp = direction.deviator();  
	A = temp && temp;
	B = - 2 * (t1 && temp);
	if (fabs(B) < LOW_LIMIT) { B = 0.; mark=1;}; 
	C = (t1 && t1) - 2./3.* size * size; 
	if ( fabs(C) < LOW_LIMIT || fabs(C)/(t1 && t1) < LOW_LIMIT ) return;

	if (B > 0. || C < 0.) {
	  opserr << "FATAL:MultiYieldSurfaceClay::updateActiveSurface(): error in surface motion.\n" 
	       << "A= " <<A <<" B= " <<B <<" C= "<<C <<" (t1&&t1)= "<<(t1&&t1) <<endln; 
	  exit(-1);
	}
	X = secondOrderEqn(A,B,C,1);  


	//---- sensitivity part --------------------------------

	dA=2.*(dMu && temp);
	dB= -2.* (dMu && t1)-2.*(temp && dTempStress);
	if (mark==1) {dB=0.; mark=0;}
	dC=2.* (dTempStress && t1)-4./3.*size*dSize;
	dX=(-dA*X*X-dB*X-dC)/(2*A*X+B);
	
	
	//center += temp * X;

	// ------- step 4. -------------------------------------

	center.addVector(1.0, temp, X);
	theSurfaces[activeSurfaceNum].setCenter(center);

	// ----------sensitivity part --------------------------

	dCenter.addVector(1.0,temp,dX);
	dCenter.addVector(1.0,dMu,X);
	for (int N=0;N<6;N++) dMultiSurfaceCenter[N+activeSurfaceNum*6+(gradNumber-1)*6*(numOfSurfaces+1)]=dCenter(N);

/* 	opserr << "step53. updateActiveSurface. center(m) is:"<< endln;
    opserr <<  center<< endln;
    opserr << "step53. updateActiveSurface. dCenter is:"<< endln;
    opserr <<  dCenter << endln;
*/

}      




void MultiYieldSurfaceClay::stressCorrectionSensitivity(int crossedSurface)
{
	// most probably, this part wrong.
	static
		T2Vector contactStress;
	this->getContactStressSensitivity(contactStress);
	static
		Vector surfaceNormal(6),dSurfaceNormal(6);
	this->getSurfaceNormalSensitivity(contactStress,dContactStress, surfaceNormal,dSurfaceNormal);
	double loadingFunc = getLoadingFuncSensitivity(contactStress, surfaceNormal,dSurfaceNormal,crossedSurface);

	//static		Vector devia(6);

	//devia = trialStress.deviator() - surfaceNormal * 2 * refShearModulus * loadingFunc;
	devia.addVector(0.0, surfaceNormal, -2*refShearModulus*loadingFunc);
	devia += trialStress.deviator();




	//----- sensitivity part ---------------------------------------
//	double dPeakShear=0.0;
	double dRefShearModulus=0.0;

// ------------ switch 3 ----------------------------
	//place 3 need to add if parameterID !=0
	if (parameterID !=0) {  //do nothing
		if(parameterID==1)  dRefShearModulus=1.0;  
//		if(parameterID==2)  dPeakShear=1.0;  
//		if(parameterID==3) ; //for Bulkmodulus, doing nothing

	} //end change
	static
		Vector dTempStress(6);
	dTempStress=surfaceNormal;
	dTempStress *= (dLoadingFunc*2.*refShearModulus+loadingFunc*2.*dRefShearModulus);
	dTempStress.addVector(1.0,dSurfaceNormal,loadingFunc*2.*refShearModulus);
	dTempStress.addVector(-1.0,dTrialStress.deviator(),1.0);
	dTrialStress.setData(dTempStress,0.0);

	//----------------------------------------------------------
	trialStress.setData(devia, 0.0);
/*
 	opserr << "step44. trialStress is:"<< endln;
    opserr <<  trialStress.deviator()<< endln;
 	opserr << "step44. dTrialStress is:"<< endln;
    opserr <<  dTrialStress.deviator()<< endln;
*/	
	deviatorScaling(trialStress, theSurfaces, activeSurfaceNum);

	if (isCrossingNextSurface()) {
		activeSurfaceNum++;
		stressCorrectionSensitivity(1);  //recursive call
	}
}


const Vector & 
MultiYieldSurfaceClay::getStressSensitivity(int passedGradNumber, 
												 bool conditional){ 


// gradNumber=passedGradNumber;

// --- to be consistent with MHS new framework ---- 2009
  passedGradNumber = passedGradNumber+1; 
  gradNumber=passedGradNumber;



/*
opserr << "---------------------------------"<< endln;
opserr << "| Begin getStressSensitivity    |"<< endln;
opserr << "---------------------------------"<< endln;
*/
// --------- sensitivity part ------------------------------------
 static 
	 T2Vector dLastStrain;   
 static 
	 Vector dTempStress(6),dTempStrain(6);
 double dLastStrainVolume;
 

 //  static Vector dStrainRate(6);


//-----gu2004.10.-- plan to change 2. Using isSurfacesSensitivitySetUp(gradNumber) to check whether this
// GradNumber's surfaceSensitivity is set up  -----

 if (this->isSurfacesSensitivitySetUp(passedGradNumber)!=1) {  // 0 means not set yet.
	 if (this->isSurfacesSensitivitySetUp(passedGradNumber)==0)   
// something wrong with setUpSurfacesSensitivity????  no
		 this->setUpSurfacesSensitivity(passedGradNumber); // 
	else {
		opserr << "isSurfacesSensitivitySetUp(gradNumber) Not Exist!" << endln;
		exit (-1);
	} //else
  }
 
 if (SHVs==0) { 
	  dTempStrain.Zero();
	  dTempStress.Zero(); 
  }


  else if (SHVs !=0){
		for (int i=0;i<6;i++)	{
			dTempStrain(i)=(*SHVs)(i,gradNumber-1);
			dTempStress(i)=(*SHVs)(i+6,gradNumber-1);
			}
	    
	}
  
//============================
  
 /* if (debugMarks ==1){
	  opserr.setPrecision(20);
	  opserr<<"getStressSens step1, dTempStress is "<<dTempStrain<<" and dTempStrain is "<<dTempStrain<<endln;
	  opserr<<"getStressSens step1, dTempStress is "<<dTempStress<<" and dTempStrain is "<<dTempStress<<endln;
  }   ==> always 0. correct  */  
 // ===========================
  // get volume part stress and strian firstly
  double dLastStressVolume=(dTempStress(0)+dTempStress(1)+dTempStress(2))/3.;
         dLastStrainVolume=(dTempStrain(0)+dTempStrain(1)+dTempStrain(2))/3.;
  
  // get deviator part only

   dCurrentStress.setData(dTempStress,0.0);    
   dLastStrain.setData(dTempStrain,0.0);  
 
 /*  if (debugMarks ==1){
	  opserr.setPrecision(20);
	  opserr<<"getStressSens step1, dCurrentStress is "<<dCurrentStress.t2Vector()<<" and dLastStrain is "<<dLastStrain.t2Vector()<<endln;

  }  yes, 0 */

// --------------------------------------------------------	
//------------------ Program ------------------------------
  int loadStage = loadStagex[matN];
  int numOfSurfaces = numOfSurfacesx[matN];
  int ndm = ndmx[matN];

  int i;
  if (loadStage == 1 && e2p == 0) //elast2PlastSensitivity();
  {  opserr << "Fatal: can not deal with elast2plast right now" << endln;
 //    exit(-1);

  }
  if (loadStage!=1) {  //linear elastic
     opserr << "Fatal: can not deal with linear elastic material right now" << endln;
     exit(-1);
  }
 
  else {
	  for (i=1; i<=numOfSurfaces; i++)  theSurfaces[i] = committedSurfaces[i];

	  for (i=1; i<=numOfSurfaces; i++){
			for(int j=0;j<myNumGrads;j++){

				for(int k=0;k<6;k++) dMultiSurfaceCenter[k+i*6+j*6*(numOfSurfaces+1)]=dCommittedMultiSurfaceCenter[k+i*6+j*6*(numOfSurfaces+1)];
/*				if (debugMarks ==1){
					 opserr.setPrecision(20);
					 for(int k=0;k<6;k++) opserr<< dMultiSurfaceCenter[k+i*6+j*6*(numOfSurfaces+1)]<<"  ";
						 opserr<<endln;
				} //if debug     wrong, corrected!!!!*/
			}
		} //for

    activeSurfaceNum = committedActiveSurf;
    subStrainRate = strainRate;

//   dSubStrainRate.deviator=dCurrentStrain.deviator()-dStrain.deviator();

    dCurrentStrain.Zero();
    dTempStrain=dCurrentStrain.deviator();
	dTempStrain.addVector(1.0,dLastStrain.deviator(),-1.0);
    dStrainRate.setData(dTempStrain,0.0);
	dSubStrainRate.setData(dTempStrain,0.0);

// ----------------------------------------------------------------	
/*  === already wrong in 1e-14 scale.
	  if (debugMarks ==1){
		  opserr.setPrecision(20);
		  opserr<<"getStressSens step1, dTrialStress is "<<dTrialStress.t2Vector()<<endln;
	  }  */

	setTrialStress(currentStress);
    if (isLoadReversal()) {

      updateInnerSurfaceSensitivity(); 


      activeSurfaceNum = 0;
    }

    int numSubIncre = setSubStrainRateSensitivity();
    
    for (i=0; i<numSubIncre; i++) {
      if (i==0) 

	   setTrialStressSensitivity(currentStress,dCurrentStress);
      else 

	   setTrialStressSensitivity(trialStress,dTrialStress); 

      if (activeSurfaceNum==0 && !isCrossingNextSurface()) continue;
      if (activeSurfaceNum==0) activeSurfaceNum++;

	  stressCorrectionSensitivity(0); 

    }
    
	//----------volume stress change and its sensitivity-----
    double volum = refBulkModulus*(strainRate.volume()*3.);
    volum += currentStress.volume();

//    Vector temp(6);
	temp.addVector(0.0,trialStress.deviator(),1.0);
	trialStress.setData(temp,volum);


//    opserr<<" ========== activeSurfaceNum:"<<activeSurfaceNum<<endln;

	double dRefBulkModulus =0.0;

//---- switch 4 ----------------	
	if (parameterID !=0) {  


		if(parameterID==3) dRefBulkModulus=1.0 ;
  
	} //end change	


//   
	double dCurrentStressVolume=dLastStressVolume-3.0*refBulkModulus*dLastStrainVolume+3*dRefBulkModulus*subStrainRate.volume();

//	dTrialStress.setData(dTrialStress.deviator(),dCurrentStressVolume);	

	Vector dTemp6(dTrialStress.deviator());
	dTrialStress.setData(dTemp6,dCurrentStressVolume);	


  }
//  opserr << "---------------------------------"<< endln;
//  opserr << "| End getStressSensitivity    |"<< endln;
//  opserr << "---------------------------------"<< endln;



//	static Vector temp6(6);

	temp6.addVector(0.0,dTrialStress.t2Vector(),1.0);





  if (ndm==3)
    //return dTrialStress.t2Vector();   ~~~~~guquan ~~~~~ change by following line
	return temp6;
  else {
    static Vector workV(3);
    workV[0] = temp6[0];
    workV[1] = temp6[1];
    workV[2] = temp6[3];
    return workV;
  }    

	}







	
int MultiYieldSurfaceClay::commitSensitivity (const Vector & strainSens, int passedGradNumber, int numGrads) {

//	gradNumber=passedGradNumber;  

// ---- to be consistent with MHS framework ----
		
    passedGradNumber=passedGradNumber+1;
	gradNumber = passedGradNumber;


  int ndm = ndmx[matN];

  static Vector strainSensitivity(6);
  if (ndm==3 && strainSens.Size()==6) 
    strainSensitivity = strainSens;
  else if (ndm==2 && strainSens.Size()==3) {
    strainSensitivity[0] = strainSens[0];
    strainSensitivity[1] = strainSens[1];
    strainSensitivity[2] = 0.0;
    strainSensitivity[3] = strainSens[2];
    strainSensitivity[4] = 0.0;
    strainSensitivity[5] = 0.0;
  }
  else {
    opserr << "Fatal:D2PressDepMYS:: Material dimension is: " << ndm << endln;
    opserr << "But strain vector size is: " << strainSens.Size() << endln;
    exit(-1);
  }




	if (SHVs == 0) {
		SHVs = new Matrix(12,numGrads);


		(*SHVs).Zero();





//-----gu2004.10.-- plan to change 1. extend to real size. myNumGrads by default is set to 1.-----
//
//  what array we need extend ?		
//  1. SHVs (12, NumGrads)   ---- no need !!

//	2.	dMultiSurfaceCenter = new double [6*(numOfSurfaces+1)*myNumGrads];
//	3.	dCommittedMultiSurfaceCenter = new double [6*(numOfSurfaces+1)*myNumGrads];

//	4.	dCommittedMultiSurfaceSize=new double [(numOfSurfaces+1)*myNumGrads];
//	5.	dCommittedMultiSurfacePlastModul=new double [(numOfSurfaces+1)*myNumGrads];
//	6. 	surfacesSensitivityMark = new int (myNumGrads);
//
//		
	if (numGrads>myNumGrads) {


	
	int numOfSurfaces = numOfSurfacesx[matN];
    double * dTemp;
	int * dTemp1;
	dTemp=new double [6*(numOfSurfaces+1)*myNumGrads];
    int i;
// ----------- copy old data	----------------------------------------------
	for( i=0; i<6*(numOfSurfaces+1)*myNumGrads;i++) dTemp[i]=dMultiSurfaceCenter[i];
	delete [] dMultiSurfaceCenter;
	dMultiSurfaceCenter = new double [6*(numOfSurfaces+1)*numGrads];


    for( i=0; i<6*(numOfSurfaces+1)*myNumGrads;i++) dMultiSurfaceCenter[i]=dTemp[i];
    for( i=6*(numOfSurfaces+1)*myNumGrads; i<6*(numOfSurfaces+1)*numGrads;i++) dMultiSurfaceCenter[i]=0;
   
//-------------------------------------------------------------------------------
	for( i=0; i<6*(numOfSurfaces+1)*myNumGrads;i++) dTemp[i]=dCommittedMultiSurfaceCenter[i];
	delete [] dCommittedMultiSurfaceCenter;
	dCommittedMultiSurfaceCenter = new double [6*(numOfSurfaces+1)*numGrads];
    for( i=0; i<6*(numOfSurfaces+1)*myNumGrads;i++) dCommittedMultiSurfaceCenter[i]=dTemp[i];
    for( i=6*(numOfSurfaces+1)*myNumGrads; i<6*(numOfSurfaces+1)*numGrads;i++) dCommittedMultiSurfaceCenter[i]=0;
  
	

	delete [] dTemp;
	
//--------------------------------------------------------------------------------
	dTemp=new double [(numOfSurfaces+1)*myNumGrads];
	for( i=0;i<(numOfSurfaces+1)*myNumGrads;i++) dTemp[i]=dCommittedMultiSurfaceSize[i];
	delete [] dCommittedMultiSurfaceSize;
	dCommittedMultiSurfaceSize=new double [(numOfSurfaces+1)*numGrads];
    for( i=0;i<(numOfSurfaces+1)*myNumGrads;i++) dCommittedMultiSurfaceSize[i]=dTemp[i];
    for( i=(numOfSurfaces+1)*myNumGrads;i<(numOfSurfaces+1)*numGrads;i++) dCommittedMultiSurfaceSize[i]=0;

//--------------------------------------------------------------------------------
	
	for( i=0;i<(numOfSurfaces+1)*myNumGrads;i++) dTemp[i]=dCommittedMultiSurfacePlastModul[i];
	delete [] dCommittedMultiSurfacePlastModul;
	dCommittedMultiSurfacePlastModul=new double [(numOfSurfaces+1)*numGrads];
    for( i=0;i<(numOfSurfaces+1)*myNumGrads;i++) dCommittedMultiSurfacePlastModul[i]=dTemp[i];
    for( i=(numOfSurfaces+1)*myNumGrads;i<(numOfSurfaces+1)*numGrads;i++) dCommittedMultiSurfacePlastModul[i]=0;

	delete [] dTemp;

//----------------------------------------------------------------------------------
	dTemp1=new int [myNumGrads];
	for( i=0; i<myNumGrads;i++) dTemp1[i]=surfacesSensitivityMark[i];
	delete [] surfacesSensitivityMark;
	surfacesSensitivityMark = new int [numGrads];
	for( i=0; i<myNumGrads;i++) surfacesSensitivityMark[i]=dTemp1[i];
	for( i=myNumGrads;i<numGrads;i++) surfacesSensitivityMark[i]=0;

	delete [] dTemp1;


//--------------------------------------------------------------------------------
	myNumGrads=numGrads;


	}   // end 	if (numGrads>myNumGrads)



			
} // end if (SHVs ==0)


//debugMarks=1;

// --------- sensitivity part ------------------------------------
 static
	 T2Vector dLastStrain;   
 static
	 Vector dTempStress(6),dTempStrain(6);
 double dLastStrainVolume;
 double dCurrentStrainVolume;

		for (int i=0;i<6;i++)	{
			dTempStrain(i)=(*SHVs)(i,gradNumber-1);
			dTempStress(i)=(*SHVs)(i+6,gradNumber-1);
			}
        
//============================
/*		if (debugMarks ==1){
			opserr.setPrecision(20);
			opserr<<"commitStressSens step1a, dTempStress is "<<dTempStress<<" and dTempStrain is "<<dTempStrain<<endln;
		}  */
 // ===========================
		

 // get volume part stress and strian firstly
  double dLastStressVolume=(dTempStress(0)+dTempStress(1)+dTempStress(2))/3.;
   dLastStrainVolume=(dTempStrain(0)+dTempStrain(1)+dTempStrain(2))/3.;
   dCurrentStrainVolume=(strainSensitivity(0)+strainSensitivity(1)+strainSensitivity(2))/3.;


//    opserr << "-------------------------------------------"<< endln;
//    opserr << "| step AAA the currentStrainSensitivity(total)|"<< endln;
//    opserr << "-------------------------------------------"<< endln;
//	opserr << strainSensitivity<< endln;
 
	static Vector tempaa(6);
	tempaa = currentStrain.t2Vector();
	 tempaa += strainRate.t2Vector();
  
//    opserr << "| step AAA the currentStrain(total)  is   |"<< endln;
//    opserr << "-------------------------------------------"<< endln;
 
  // get deviator part only


   dCurrentStress.setData(dTempStress,0.0);    
   dLastStrain.setData(dTempStrain,0.0);  



// --------------------------------------------------------	
//------------------ Program ------------------------------
  int loadStage = loadStagex[matN];
  int numOfSurfaces = numOfSurfacesx[matN];
//  int ndm = ndmx[matN];

  if (loadStage == 1 && e2p == 0) //elast2PlastSensitivity();
  {  opserr << "Fatal: can not deal with elast2plast right now" << endln;
     exit(-1);
  }
  if (loadStage!=1) {  //linear elastic
     opserr << "Fatal: can not deal with linear elastic material right now" << endln;
     exit(-1);
  }
 
  else {
	  
int i;
	  for (i=1; i<=numOfSurfaces; i++)  theSurfaces[i] = committedSurfaces[i];


	  for (i=1; i<=numOfSurfaces; i++){
			for(int j=0;j<myNumGrads;j++){
				
				for(int k=0;k<6;k++) dMultiSurfaceCenter[k+i*6+j*6*(numOfSurfaces+1)]=dCommittedMultiSurfaceCenter[k+i*6+j*6*(numOfSurfaces+1)];

			}
		} //for

    activeSurfaceNum = committedActiveSurf;
    subStrainRate = strainRate;


//   dSubStrainRate.deviator=dCurrentStrain.deviator()-dStrain.deviator();
//    dCurrentStrain.setData(strainGradient);  anyproblem?

	dCurrentStrain.setData(strainSensitivity,1);
    dTempStrain=dCurrentStrain.deviator();
	dTempStrain.addVector(1.0,dLastStrain.deviator(),-1.0);
    dStrainRate.setData(dTempStrain,0.0);
	dSubStrainRate.setData(dTempStrain,0.0);
// ----------------------------------------------------------------	



	setTrialStress(currentStress);
    if (isLoadReversal()) {

      updateInnerSurfaceSensitivity();
      activeSurfaceNum = 0;
    }
    int numSubIncre = setSubStrainRateSensitivity();

    for (i=0; i<numSubIncre; i++) {
      if (i==0) 

	   setTrialStressSensitivity(currentStress,dCurrentStress);

      else 

	   setTrialStressSensitivity(trialStress,dTrialStress);
      if (activeSurfaceNum==0 && !isCrossingNextSurface()) continue;
      if (activeSurfaceNum==0) activeSurfaceNum++;

      stressCorrectionSensitivity(0);
      updateActiveSurfaceSensitivity();



    }
    
	//----------volume stress change and its sensitivity-----
    double volum = refBulkModulus*(strainRate.volume()*3.);


    volum += currentStress.volume();

//	Vector temp(6);
	temp.addVector(0.0,trialStress.deviator(),1.0);
    trialStress.setData(temp,volum);

	double dRefBulkModulus =0.0;
	if (parameterID !=0) {  

// ----------- switch 5 here ---------------------
	if(parameterID==3)  dRefBulkModulus =1.0;
		
  	}

//   
//	double dCurrentStressVolume=dLastStressVolume-3.0*refBulkModulus*dLastStrainVolume+3*dRefBulkModulus*subStrainRate.volume();
	

	
	//-------------------------------------------------------
	double dCurrentStressVolume=dLastStressVolume+3.0*refBulkModulus*(dCurrentStrainVolume-dLastStrainVolume)+3*dRefBulkModulus*subStrainRate.volume();

//	Vector temp6(6);
	temp6.addVector(0.0,dTrialStress.deviator(),1.0);
	dTrialStress.setData(temp6,dCurrentStressVolume);


  }
// -----------------------add committed part------------------
//    opserr << "---------------------------------"<< endln;
//    opserr << "| End commitStressSensitivity    |"<< endln;
//    opserr << "---------------------------------"<< endln;

  	  for (int i=1; i<=numOfSurfaces; i++){
			for(int j=0;j<myNumGrads;j++){
				for(int k=0;k<6;k++) dCommittedMultiSurfaceCenter[k+i*6+j*6*(numOfSurfaces+1)]=dMultiSurfaceCenter[k+i*6+j*6*(numOfSurfaces+1)];
			}
		} //for

	  dCurrentStress=dTrialStress;

	  dTempStrain=dCurrentStrain.t2Vector(0);
	  dTempStress=dCurrentStress.t2Vector(0);

 /* if (debugMarks ==1){
	  opserr.setPrecision(20);
	  opserr<<"commitStressSens step2, dTempStrain is "<<dTempStrain<<endln;
	  opserr<<"commitStressSens step2, dTempStress is "<<dTempStress<<endln;
  }
*/

	  for(int i=0;i<6;i++){
	  
		  	(*SHVs)(i,gradNumber-1)=dTempStrain(i);
			(*SHVs)(i+6,gradNumber-1)=dTempStress(i);

		 }
	  
return 0;
//*/    //guquan
	}




int MultiYieldSurfaceClay::isSurfacesSensitivitySetUp(int passedGradNumber) {

// 0: Not set yet
// 1: done

// ---- to be consistent with MHS new framework ---- 2009
// === Not added by one here. But when call this function.
//    passedGradNumber +=1;

	if (surfacesSensitivityMark==0) {  // create surfacesSensitivityMark array
		surfacesSensitivityMark = new int [myNumGrads];
		for (int i=0;i<myNumGrads;i++)		surfacesSensitivityMark[i]=0;
	}



	if (surfacesSensitivityMark[passedGradNumber-1] ==0) return 0;
	else if (surfacesSensitivityMark[passedGradNumber-1] ==1) return 1;
	else {
		opserr << "MultiYieldSurfaceClay::isSurfacesSensitivitySetup, surfacesSensitivityMark(passsedGradNumber) NOT EXIST!"<< endln;	
		exit (-1);
	}
	
}


void MultiYieldSurfaceClay::setSurfacesSensitivitySetUpMark(int passedGradNumber){
	if (surfacesSensitivityMark[passedGradNumber-1] !=0) {
	opserr <<"Error! this surface sensitivity already set" << endln;
	exit(-1);
	}
	surfacesSensitivityMark[passedGradNumber-1]=1;

}

//Feb .1 05  need to change
const Vector &MultiYieldSurfaceClay::getCommittedStressSensitivity(int GradientNumber){


	int ndm = ndmx[matN];
//	static Vector temp6(6);
	temp6.Zero();
	int i;

	if (SHVs !=0){
		for(i=0;i<6;i++)  temp6(i)=(*SHVs)(i+6,GradientNumber-1);	 
	}
	if (ndm==3)
		return temp6;
	
	else {
    static Vector workV(3);
	
    workV[0] = temp6[0];
    workV[1] = temp6[1];
    workV[2] = temp6[3];
    return workV;
  }


}





const Vector &MultiYieldSurfaceClay::getCommittedStrainSensitivity(int GradientNumber){

	int ndm = ndmx[matN];
//	static Vector temp6(6);
	temp6.Zero();
	int i;
	// obtain data from SHVs into temp need attation that T2Vector(0) or T2Vector(1)
	if (SHVs !=0){
		for(i=0;i<6;i++) 	temp6(i)=(*SHVs)(i,GradientNumber-1);	
	}
	if (ndm==3)
		return temp6;
	
	else {
    static Vector workV(3);
	
    workV[0] = temp6[0];
    workV[1] = temp6[1];
    workV[2] = temp6[3];
    return workV;
  }


}




///////////////////////// end sensitivity /////////////////////////////////////


/*  plan to change into classwide
static Vector     temp6(6);
static Vector     temp(6);
static Vector     devia(6);

  */

