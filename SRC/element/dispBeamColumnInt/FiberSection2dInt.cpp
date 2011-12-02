// $Source: /usr/local/cvs/OpenSees/SRC/element/dispBeamColumnInt/FiberSection2dInt.cpp,v $

// $Revision: 1.4 $

// $Date: 2010-03-25 23:59:38 $



// Created: 07/04

// Modified by: LMS 

// Description: This file contains the class implementation of FiberSection2dInt.Based on FiberSection2d.cpp.

#include <stdlib.h>
#include <math.h>

#define PI (3.14159265)

#include <Channel.h>
#include <Vector.h>
#include <Matrix.h>
#include <MatrixUtil.h>
#include <Fiber.h>
#include <classTags.h>
#include "FiberSection2dInt.h"
#include <ID.h>

#include <FEM_ObjectBroker.h>

#include <Information.h>

#include <MaterialResponse.h>

#include <UniaxialMaterial.h>



#include <float.h>



#define MyMIN(a, b)  (((a) < (b)) ? (a) : (b))

#define MyMAX(a, b)  (((a) > (b)) ? (a) : (b))



ID FiberSection2dInt::code(2);			



// constructors:

FiberSection2dInt::FiberSection2dInt(int tag, int num, Fiber **fibers, int Hnum, Fiber **Hfibers, int strip1, double t1, int strip2, double t2, int strip3, double t3): 

  SectionForceDeformation(tag, SEC_TAG_FiberSection2dInt),

  numFibers(num), theMaterials1(0), theMaterials2(0), matData(0),

  numHFibers(Hnum), theHMaterials(0), matHData(0), 

  NStrip(strip1+strip2+strip3), NStrip1(strip1), tavg1(t1), NStrip2(strip2), tavg2(t2), NStrip3(strip3), tavg3(t3), 

  StripCenterLoc(100), StripLoc(100,1000), FiberLoc(1000),

  yBar(0.0), ymax(0.0), ymin(0.0), e(3), eCommit(3), s(0), ks(0), sigmaY(0), tau(0), alpha(0), alphaCommit(0), 

  iterFile(0),exf(0),e1f(0),e2f(0),eyf(0),sxf(0),s1f(0),s2f(0),syf(0)

{

  

	if (numFibers != 0) {

    theMaterials1 = new UniaxialMaterial *[numFibers];

    theMaterials2 = new UniaxialMaterial *[numFibers];



    if (theMaterials1 == 0) {

      opserr << "FiberSection2dInt::FiberSection2dInt -- failed to allocate Material pointers";

      exit(-1);

    }



    matData = new double [numFibers*2];



    if (matData == 0) {

      opserr << "FiberSection2dInt::FiberSection2dInt -- failed to allocate double array for material data\n";

      exit(-1);

    }





    double Qz = 0.0;

    double A  = 0.0;

	 ymax=-10000.0;	//initial guess

	 ymin=10000.0;



    

    for (int i = 0; i < numFibers; i++) {

      Fiber *theFiber = fibers[i];

      double yLoc, zLoc, Area;

      theFiber->getFiberLocation(yLoc, zLoc);

      Area = theFiber->getArea();

      A  += Area;

      Qz += yLoc*Area;

      matData[i*2] = -yLoc;

      matData[i*2+1] = Area;

      UniaxialMaterial *theMat1 = theFiber->getMaterial();		

      theMaterials1[i] = theMat1->getCopy();

 

      UniaxialMaterial *theMat2 = theFiber->getMaterial();

      theMaterials2[i] = theMat2->getCopy();





      if (theMaterials1[i] == 0) {

	opserr << "FiberSection2dInt::FiberSection2dInt -- failed to get copy of a Material\n";

	exit(-1);

      }



	if (-yLoc>ymax) ymax=-yLoc;

	if (-yLoc<ymin) ymin=-yLoc;



    }    



    yBar = -Qz/A;  

  }





 



  if (numHFibers != 0) {



    theHMaterials = new UniaxialMaterial *[numHFibers*NStrip];





    if (theHMaterials == 0) {

      opserr << "FiberSection2dInt::FiberSection2dInt -- failed to allocate HMaterial pointers";

      exit(-1);

    }

    

    matHData = new double [numHFibers*2];

    

    if (matHData == 0) {

      opserr << "FiberSection2dInt::FiberSection2dInt -- failed to allocate double array for Hmaterial data\n";

      exit(-1);

    }

    

    

    double HQz = 0.0;

    double HA  = 0.0;

    

    for (int i = 0; i < numHFibers; i++) {

      Fiber *theHFiber = Hfibers[i];

      double yHLoc, zHLoc, HArea;

      theHFiber->getFiberLocation(yHLoc, zHLoc);

      HArea = theHFiber->getArea();

      HA  += HArea;

      HQz += yHLoc*HArea;

      matHData[i*2] = -yHLoc;

      matHData[i*2+1] = HArea;

      UniaxialMaterial *theHMat = theHFiber->getMaterial();

      

      for (int jj = 0; jj < NStrip; jj++) {

	theHMaterials[i*numHFibers + jj] = theHMat->getCopy();



	if (theHMaterials[i*numHFibers + jj] == 0) {

	  opserr << "FiberSection2dInt::FiberSection2dInt -- failed to get copy of a HMaterial\n";

	  exit(-1);

	}

      }

    }    



  }





  double YLoc[100];

  

  int ycount = 0;

  int loc = 0;

  for (int i = 0; i < numFibers; i++) {

    double y = matData[loc++];

    double A = matData[loc++];

    if (i==0) {

      YLoc[0] = y;

      ycount += 1; 

    }

    else {

      if (fabs(YLoc[ycount-1] - y) >= DBL_EPSILON) {

	YLoc[ycount] = y;

	ycount += 1;

      }

    }

    FiberLoc(i)=ycount-1;

  }

  

  if (ycount != NStrip) {

    opserr <<  "\n Failed - Not consistent number of fibers \n";

    exit(-1);

  }

	

  for (int j = 0; j < NStrip; j++) StripCenterLoc(j) = +(YLoc[j] - yBar);

  

  for (int k = 0; k < NStrip; k++) {

    int count=0;

    double Ac=0.0;

    for (int i = 0; i < numFibers; i++) {

      if (FiberLoc(i)==k){

	count++;	

	StripLoc(k,count+1)=i;

	Ac += matData[2*i+1];

      }

    }

    StripLoc(k,0)=count;	//num fibers in strip

    StripLoc(k,1)=Ac;		//total concrete area in strip

  }

  

  

  for (int kk = 0; kk < NStrip; kk++) {

    exCommit[kk]=0.0;

  }

  

  

  s = new Vector(sData, 3);

  ks = new Matrix(kData, 3, 3);



  sData[0] = 0.0;

  sData[1] = 0.0;

  sData[2] = 0.0;





  kData[0] = 0.0;

  kData[1] = 0.0;

  kData[2] = 0.0;

  kData[3] = 0.0;

  kData[4] = 0.0;

  kData[5] = 0.0;

  kData[6] = 0.0;

  kData[7] = 0.0;

  kData[8] = 0.0;



  code(0) = SECTION_RESPONSE_P;							

  code(1) = SECTION_RESPONSE_MZ;

  code(2) = SECTION_RESPONSE_VY;





// AddingSensitivity:BEGIN ////////////////////////////////////

	parameterID = 0;

// AddingSensitivity:END //////////////////////////////////////



}



// constructor for blank object that recvSelf needs to be invoked upon

FiberSection2dInt::FiberSection2dInt():

  SectionForceDeformation(0, SEC_TAG_FiberSection2dInt),

  numFibers(0), theMaterials1(0), theMaterials2(0), matData(0),

  numHFibers(0), theHMaterials(0), matHData(0), NStrip1(0), tavg1(0.0), NStrip2(0), tavg2(0.0), NStrip3(0), tavg3(0.0), StripCenterLoc(100), StripLoc(100,1000), FiberLoc(1000),

  yBar(0.0), ymax(0.0), ymin(0.0), e(3), eCommit(3), s(0), ks(0), sigmaY(0), tau(0), alpha(0), alphaCommit(0), iterFile(0), exf(0), e1f(0),e2f(0),eyf(0), sxf(0), s1f(0),s2f(0),syf(0)

{

  s = new Vector(sData, 3);

  ks = new Matrix(kData, 3, 3);



  sData[0] = 0.0;

  sData[1] = 0.0;

  sData[2] = 0.0;



  kData[0] = 0.0;

  kData[1] = 0.0;

  kData[2] = 0.0;

  kData[3] = 0.0;

  kData[4] = 0.0;

  kData[5] = 0.0;

  kData[6] = 0.0;

  kData[7] = 0.0;

  kData[8] = 0.0;



  code(0) = SECTION_RESPONSE_P;

  code(1) = SECTION_RESPONSE_MZ;

  code(2) = SECTION_RESPONSE_VY;





// AddingSensitivity:BEGIN ////////////////////////////////////

	parameterID = 0;

// AddingSensitivity:END //////////////////////////////////////

}



int

FiberSection2dInt::addFiber(Fiber &newFiber)		

{

  // need to create larger arrays

  int newSize = numFibers+1;

  UniaxialMaterial **newArray1 = new UniaxialMaterial *[newSize]; 

  UniaxialMaterial **newArray2 = new UniaxialMaterial *[newSize]; 



  double *newMatData = new double [2 * newSize];

  if (newArray1 == 0 || newMatData == 0) {

    opserr <<"FiberSection2d::addFiber -- failed to allocate Fiber pointers\n";

    return -1;

  }



  // copy the old pointers and data

  int i;

  for (i = 0; i < numFibers; i++) {

    newArray1[i] = theMaterials1[i];

    newArray2[i] = theMaterials2[i];



    newMatData[2*i] = matData[2*i];

    newMatData[2*i+1] = matData[2*i+1];

  }



  // set the new pointers and data

  double yLoc, zLoc, Area;

  newFiber.getFiberLocation(yLoc, zLoc);

  Area = newFiber.getArea();

  newMatData[numFibers*2] = -yLoc;

  newMatData[numFibers*2+1] = Area;

  UniaxialMaterial *theMat = newFiber.getMaterial();

  newArray1[numFibers] = theMat->getCopy();

  newArray2[numFibers] = theMat->getCopy();



  if (newArray1[numFibers] == 0) {

    opserr <<"FiberSection2dInt::addFiber -- failed to get copy of a Material\n";

    delete [] newMatData;

    return -1;

  }



  numFibers++;



  if (theMaterials1 != 0) {

    delete [] theMaterials1;

	delete [] theMaterials2;

    delete [] matData;

  }



  theMaterials1 = newArray1;

  theMaterials2 = newArray2;

  matData = newMatData;



  double Qz = 0.0;

  double A  = 0.0;



  ymax=-10000.0;	//initial guess

  ymin=10000.0;

  

  // Recompute centroid

  for (i = 0; i < numFibers; i++) {

    yLoc = -matData[2*i];

    Area = matData[2*i+1];

    A  += Area;

    Qz += yLoc*Area;

	if (-yLoc>ymax) ymax=-yLoc;

	if (-yLoc<ymin) ymin=-yLoc;



  }



  yBar = -Qz/A;







  double YLoc[100];

  

  int ycount = 0;

  int loc = 0;

  for (int ii = 0; ii < numFibers; ii++) {

    double y = matData[loc++];

    double A = matData[loc++];

    if (ii==0) {

      YLoc[0] = y;

      ycount += 1; 

    }

    else {

      if (fabs(YLoc[ycount-1] - y) >= DBL_EPSILON) {

	YLoc[ycount] = y;

	ycount += 1;

      }

    }

    FiberLoc(ii)=ycount-1;

  }

  

  if (ycount != NStrip) {

    opserr <<  "\n Failed - Not consistent number of fibers \n";

    exit(-1);

  }

  

  for (int j = 0; j < NStrip; j++) StripCenterLoc(j) = YLoc[j] - yBar;

  

  

  for (int k = 0; k < NStrip; k++) {

    int count=0;

    double Ac=0.0;

    for (int i = 0; i < numFibers; i++) {

      if (FiberLoc(i)==k){

	count++;	

	StripLoc(k,count+1)=i;

	Ac += matData[2*i+1];

      }

    }

    StripLoc(k,0)=count;	//num fibers in strip

    StripLoc(k,1)=Ac;		//total concrete area in strip

  }

  

  return 0;

}







// destructor:

FiberSection2dInt::~FiberSection2dInt()

{

  if (theMaterials1 != 0) {

    for (int i = 0; i < numFibers; i++)

      if (theMaterials1[i] != 0) {

	delete theMaterials1[i];

	delete theMaterials2[i];

      }

    delete [] theMaterials1;

    delete [] theMaterials2;

  }

  

  if (matData != 0)

    delete [] matData;

  

  

  if (theHMaterials != 0) {

    for (int i = 0; i < numHFibers; i++)

      if (theHMaterials[i*numHFibers +0] != 0)

	for (int jj = 0; jj < NStrip; jj++) 

	  if (theHMaterials[i * numHFibers + jj] != 0)

	    delete theHMaterials[i * numHFibers + jj];

    

    delete [] theHMaterials;

  }

  

  if (matHData != 0)

    delete [] matHData;

  

  

  if (s != 0)

    delete s;

  

  if (ks != 0)

    delete ks;

  

  if (sigmaY != 0)

    delete sigmaY;

  

  if (tau != 0)

    delete tau;

  

  if (alpha != 0)

    delete alpha;

  

  if (alphaCommit != 0)

    delete alphaCommit;

  

  if (iterFile != 0)

    delete iterFile;

  

  if (exf != 0)

    delete exf;

  

  if (e1f != 0)

    delete e1f;

  

  if (e2f != 0)

    delete e2f;

  

  if (eyf != 0)

    delete eyf;

  

  if (sxf != 0)

    delete sxf;

  

  if (s1f != 0)

    delete s1f;

  

  if (s2f != 0)

    delete s2f;

  

  if (syf != 0)

    delete syf;

}





int FiberSection2dInt::setTrialSectionDeformation (const Vector &deforms)	

{

  return 0;

}



int FiberSection2dInt::revertToLastCommit (void)	

{

  return 0;

}



int FiberSection2dInt::revertToStart (void)	

{

  return 0;

}





int

FiberSection2dInt::commitState(void)			

{

  return 0;

}





const Vector&

FiberSection2dInt::getSigmaY(void)	

{

  return *sigmaY;

}



const Vector&

FiberSection2dInt::getTau(void)	

{

  return *tau;

}



const Vector&

FiberSection2dInt::getAlpha(void)	

{

  return *alpha;

}



const Vector&

FiberSection2dInt::getIter(void)	

{

  return *iterFile;

}



const Vector&

FiberSection2dInt::getEX(void)	

{

  return *exf;

}



const Vector&

FiberSection2dInt::getEY(void)	

{

  return *eyf;

}



const Vector&

FiberSection2dInt::getE1(void)	

{

  return *e1f;

}



const Vector&

FiberSection2dInt::getE2(void)	

{

  return *e2f;

}







const Vector&

FiberSection2dInt::getSX(void)	

{

  return *sxf;

}



const Vector&

FiberSection2dInt::getSY(void)	

{

  return *syf;

}



const Vector&

FiberSection2dInt::getS1(void)

{

  return *s1f;

}



const Vector&

FiberSection2dInt::getS2(void)	

{

  return *s2f;

  

}







void

FiberSection2dInt::beta(double e0, double e2, double &sc1, double &tc1, double &tc12, double &beta)	

{

	

double Kc=0.27*(-e2/e0-0.37);	//compression softening parameter by Vecchio and Collins (1993)

beta= 1.0/(1.0+Kc);

double delbeta=beta*beta*0.27/e0;



if ((beta>1.0)||(e2<0.0)) {

    beta=1.0;

    delbeta=0.0;

}



if (sc1>0.0) {

    beta=1.0;

    delbeta=0.0;

}



tc12 = delbeta*sc1;

sc1 = beta*sc1;

tc1 = beta*tc1;



}





const Vector &
FiberSection2dInt::getSectionDeformation (void)
{
  return e;
}



int FiberSection2dInt::setTrialSectionDeformationB (const Vector &deforms, double L)

{

  int res = 0;


  
  e = deforms;	// axial strain, curvature and shear strain
  
  
  kData[0] = 0.0; kData[1] = 0.0; kData[2] = 0.0; kData[3] = 0.0; kData[4] = 0.0; kData[5] = 0.0;
  kData[6] = 0.0; kData[7] = 0.0; kData[8] = 0.0;
  sData[0] = 0.0; sData[1] = 0.0; sData[2] = 0.0;

  double d0 = deforms(0);						// axial strain
  double d1 = deforms(1);						// curvature
  double d2 = deforms(2);						// shear strain
  double iterMax = 200;
  double tol = 1e-8;
  

  for (int jj = 0; jj < NStrip; jj++) {

    

    double tavg;

    if (jj<NStrip1)

      tavg=tavg1;

    else {

      if (jj<NStrip2+NStrip1)

	tavg=tavg2;

      else

	tavg=tavg3;

    }

    

    double ey=d0 + StripCenterLoc(jj)*d1;

    double gamma=d2;

    

    if ((fabs(gamma) <= DBL_EPSILON) && (fabs(ey) <= DBL_EPSILON) ) {

      *ks = this->	getInitialTangent();

      sData[0] = 0.0;

      sData[1] = 0.0;

      sData[2] = 0.0;

      sigmaY = new Vector(sy, NStrip);

      tau = new Vector(txy, NStrip);

      alpha = new Vector(alfa, NStrip);

      alphaCommit = new Vector(alfaCommit, NStrip);

      iterFile = new Vector(iterOut, NStrip);

      exf = new Vector(exCommit, NStrip);

      eyf = new Vector(eyCommit, NStrip);

      e1f = new Vector(e1Commit, NStrip);

      e2f = new Vector(e2Commit, NStrip);

      

      sxf = new Vector(sxCommit, NStrip);

      syf = new Vector(syCommit, NStrip);

      s1f = new Vector(s1Commit, NStrip);

      s2f = new Vector(s2Commit, NStrip);

      

      break;

    }

    

    //	double ex=0;

    double ex=exCommit[jj];

    

    double XPrevPos = -1.0; 

    double XPrevNeg = -1.0;   

    

    double root;

    double e1;                          

    double e2; 

    double fx=0.0;

    double Fy;

    double Fxy;

    

    double Xmax = +1.0;     

    double Xmin = -1.0;     

    

    

    

    double stifstxF, stifstyF, stifcu11F, stifcu12F, stifcu22F, stifcu21F, Fcu1F, Fcu2F;



    int iter;

    for (iter = 1; iter <= iterMax; iter++) {

      

      root = sqrt(pow(-ex + ey,2) + pow(gamma,2));

      

      if (fabs(gamma) <= DBL_EPSILON) {

	if ((ey - ex)*eCommit(2) > 0) {

	  alfa[jj]=PI/2.0;

	  e1 = ex; 

	  e2 = ey; 

	}

	else {

	  alfa[jj]=0.0;

	  e1 = ey; 

	  e2 = ex; 

	}

      }

      else {

	alfa[jj]=atan((-ex + ey)/gamma + sqrt(pow((-ex + ey)/gamma,2) + 1));

	

	if (fabs(alfa[jj]-PI/2.0) <= DBL_EPSILON) {

	  e1 = ex; 

	  e2 = ey; 

	}

	else {

	  e1 = ey - gamma/2.0*tan(alfa[jj]);                          

	  e2 = ex + gamma/2.0*tan(alfa[jj]);

	}

	

      }

      

      

      fx =0.0;

      Fy = 0.0;

      Fxy = 0.0;

      

      syCommit[jj]=0.0;

      s1Commit[jj]=0.0;

      s2Commit[jj]=0.0;

      sxCommit[jj] = 0.0;

      

      double HAtot=0.0;

      double Atot=0.0;

      double Actot=0.0;

      

      stifstxF =0.0;

      stifstyF =0.0;

      stifcu11F =0.0;

      stifcu12F =0.0;

      stifcu22F =0.0;

      stifcu21F =0.0;

      Fcu1F =0.0;

      Fcu2F =0.0;

      

      

      for (int i = 2; i <= StripLoc(jj,0)+1; i++) {

	

	int fibNum=StripLoc(jj,i);

	UniaxialMaterial *theMat1 = theMaterials1[fibNum];  

	UniaxialMaterial *theMat2 = theMaterials2[fibNum];  

	

	int tag=theMat1->getTag();

	double tsy, ssy;

	double tc1, sc1, tc2, sc2;

	double tc12=0.0;

	double tc21=0.0;

	double beta1, beta2;

	

	if (tag>1000){			// to distinguish concrete & steel tag>1000 => steel

	  //	  double y = matData[2*fibNum] - yBar;		

	  double A = matData[2*fibNum+1];	

	  res = theMat1->setTrial(ey, ssy, tsy);

	  

	  fx += 0.0;

	  Fy += ssy*A;

	  Fxy += 0.0;

	  syCommit[jj] += ssy*A;

	  Atot += A; 

	  stifstyF +=tsy*A;

	  

	}

	else{			// concrete

	  //double y = matData[2*fibNum] - yBar;

	  double A = matData[2*fibNum+1];	

	  res = theMat1->setTrial(e1, sc1, tc1);

	  res = theMat2->setTrial(e2, sc2, tc2);



	  static Information theInfo;

	  double e0 = 0.0;

	  

	  const char *theData = "ec";

	  if (theMat1->getVariable(theData, theInfo) == 0)

	    e0 = theInfo.theDouble;



	  this -> beta(e0, e2, sc1, tc1, tc12, beta1);

	  this -> beta(e0, e1, sc2, tc2, tc21, beta2);

	  

	  s1Commit[jj] += sc1*A;

	  s2Commit[jj] += sc2*A;

	  Actot += A; 

	  

	  fx += ((sc1+sc2)*0.5-(sc1-sc2)*0.5*cos(2.0*alfa[jj]))*A;

	  Fy += ((sc1+sc2)*0.5+(sc1-sc2)*0.5*cos(2.0*alfa[jj]))*A;

	  Fxy += -((sc1-sc2)*0.5*sin(2.0*alfa[jj]))*A;

	  

	  stifcu11F +=tc1*A;

	  stifcu12F +=tc12*A; 

	  stifcu22F +=tc2*A;

	  stifcu21F +=tc21*A; 

	  

	  Fcu1F += sc1*A;

	  Fcu2F += sc2*A;

	  

	}

	

      }

      

      int Hloc = 0;

      for (int H = 0; H < numHFibers; H++) {				

	UniaxialMaterial *theHMat = theHMaterials[H*numHFibers + jj];

	double Hy = matHData[Hloc++];	

	double HA = matHData[Hloc++];

	double Ht, Hs;

	res = theHMat->setTrial(ex, Hs, Ht);	

	fx += Hs*HA*StripLoc(jj,1)/(L*tavg);

	HAtot += HA;

	sxCommit[jj] += Hs*HA*StripLoc(jj,1)/(L*tavg);

	stifstxF += Ht*HA*StripLoc(jj,1)/(L*tavg);

	

      }

      

      sxCommit[jj] = sxCommit[jj]/HAtot;

      syCommit[jj] = syCommit[jj]/Atot;

      s1Commit[jj] = s1Commit[jj]/Actot;

      s2Commit[jj] = s2Commit[jj]/Actot;

      

      // ex iteration...

      

      double signGamma;

      if (fabs(gamma) <= DBL_EPSILON) signGamma=1.0;

      else signGamma=fabs(gamma)/gamma;

      

      double de2 = (1.0 - (-ex + ey)*signGamma/root)*0.5;

      double de1 = (1.0 + (-ex + ey)*signGamma/root)*0.5;

      double dalfa = -gamma*0.5/pow(root,2.0);

      

      double Bgrad=stifstxF + ((stifcu11F*de1+stifcu12F*de2+stifcu22F*de2+stifcu21F*de1)/2.0-

			       (stifcu11F*de1+stifcu12F*de2-(stifcu22F*de2+stifcu21F*de1))/2.0*cos(2.0*alfa[jj])+

			       (Fcu1F-Fcu2F)*sin(2.0*alfa[jj])*dalfa);                   

      

      

      // Check the residual

      double err = fabs(fx);

      

      if (err < tol) break;

      

      if (fx>0.0) XPrevPos=ex;

      else XPrevNeg=ex;

      

      if (!(XPrevPos == -1.0) && !(XPrevNeg == -1.0)) { // assumes 1 solution available

	

	Xmax=MyMIN(MyMAX(XPrevPos,XPrevNeg),Xmax);

	Xmin=MyMAX(MyMIN(XPrevPos,XPrevNeg),Xmin);

      }

      

      if ((iter>50) && (!(XPrevPos == -1.0) && !(XPrevNeg == -1.0))) {	

	ex = (Xmax + Xmin)/2.0;

      }														// midpoint

      else {														 

	double dex = -fx/Bgrad ;									

	ex = ex + dex;

      }

    }

    

    iterOut[jj]=iter;

    exOut[jj]=ex;

    eyCommit[jj]=ey;

    e1Commit[jj]=e1;

    e2Commit[jj]=e2;

    

    double dTdEy;

    double dTdGamma;

    double dSydGamma;

    double dSydEy;

    

    if ((fabs(gamma) <= DBL_EPSILON)||(fabs(alfa[jj]) <= DBL_EPSILON)||(fabs(alfa[jj]-PI/2.0) <= DBL_EPSILON)) {

      dTdEy = 0;

      dTdGamma = (-Fcu1F + Fcu2F)/(2.*sqrt(pow(ex - ey,2)));

      dSydGamma = 0;

      dSydEy = stifcu22F - (stifcu12F*stifcu21F)/(stifcu11F + stifstxF) + stifstyF;

    }

    else {

      

      double R2=sqrt(1.0 + pow(ex - ey,2)/pow(gamma,2));

      

      dTdEy = (-((Fcu1F - Fcu2F)*(pow(ex,2) - 2*ex*ey + pow(ey,2) + pow(gamma,2))*

		 (stifcu11F + stifcu12F - 3*stifcu21F - 3*stifcu22F - 2*stifstxF)) + 

	       4*(Fcu1F - Fcu2F)*(pow(ex - ey,2) + pow(gamma,2))*(stifcu21F + stifcu22F + stifstxF)*cos(2*alfa[jj]) + 

	       (Fcu1F - Fcu2F)*(pow(ex - ey,2) + pow(gamma,2))*(stifcu11F + stifcu12F + stifcu21F + stifcu22F + 2*stifstxF)*cos(4*alfa[jj]) + 

	       2*gamma*(pow(ex - ey,2) + pow(gamma,2))*(stifcu11F + stifcu12F - stifcu21F - stifcu22F)*stifstxF*sin(2*alfa[jj]) + 

	       (-ex + ey)*gamma*R2*((-Fcu1F + Fcu2F)*(stifcu11F + stifcu12F - 3*stifcu21F - 3*stifcu22F - 2*stifstxF) + 

				    4*(Fcu1F - Fcu2F)*(stifcu21F + stifcu22F + stifstxF)*cos(2*alfa[jj]) + 

				    (Fcu1F - Fcu2F)*(stifcu11F + stifcu12F + stifcu21F + stifcu22F + 2*stifstxF)*cos(4*alfa[jj]) + 

				    2*gamma*(2*stifcu12F*stifcu21F - 2*stifcu11F*stifcu22F + (-stifcu11F + stifcu12F + stifcu21F - stifcu22F)*stifstxF)*

				    sin(2*alfa[jj])))/(2.*(-2*gamma*(pow(ex - ey,2) + pow(gamma,2))*stifcu21F*pow(cos(alfa[jj]),2) + 

							   (-ex + ey)*gamma*R2*(2*gamma*(-stifcu21F + stifcu22F)*pow(cos(alfa[jj]),2) 

										+ 8*(Fcu1F - Fcu2F)*pow(cos(alfa[jj]),3)*sin(alfa[jj]) + 

										2*gamma*(-stifcu11F + stifcu12F)*pow(sin(alfa[jj]),2)) + 

							   (-pow(-ex + ey,2) - pow(gamma,2))*(8*(-Fcu1F + Fcu2F)*pow(cos(alfa[jj]),3)*sin(alfa[jj]) + 

											      2*gamma*(2*stifstxF + stifcu22F*pow(cos(alfa[jj]),2) +

												       (stifcu11F + stifcu12F)*pow(sin(alfa[jj]),2)))));

      

      

      dTdGamma = ((-ex + ey)*(Fcu1F - Fcu2F)*(pow(ex - ey,2) + pow(gamma,2))*(-stifcu11F + 3*stifcu21F + 2*stifstxF) + 

		  4*(-ex + ey)*(Fcu1F - Fcu2F)*(pow(ex - ey,2) + pow(gamma,2))*(stifcu21F + stifstxF)*cos(2*alfa[jj]) + 

		  (-ex + ey)*(Fcu1F - Fcu2F)*(pow(ex - ey,2) + pow(gamma,2))*(stifcu11F + stifcu21F + 2*stifstxF)*cos(4*alfa[jj]) + 

		  (gamma*R2*((Fcu1F - Fcu2F)*(-((2*pow(ex - ey,2) + pow(gamma,2))*(stifcu11F - 3*stifcu21F)) + 

					      pow(gamma,2)*(stifcu12F - 3*stifcu22F) + 4*pow(ex - ey,2)*stifstxF) + 

			     4*(Fcu1F - Fcu2F)*(pow(gamma,2)*(stifcu21F - stifcu22F) + 2*pow(ex,2)*(stifcu21F + stifstxF) - 

						4*ex*ey*(stifcu21F + stifstxF) + 2*pow(ey,2)*(stifcu21F + stifstxF))*cos(2*alfa[jj]) + 

			     (Fcu1F - Fcu2F)*((2*pow(ex - ey,2) + pow(gamma,2))*(stifcu11F + stifcu21F) - pow(gamma,2)*(stifcu12F + stifcu22F) + 

					      4*pow(ex - ey,2)*stifstxF)*cos(4*alfa[jj]) + 

			     4*pow(gamma,3)*((-stifcu21F + stifcu22F)*stifstxF - stifcu12F*(stifcu21F + stifstxF) + 

					     stifcu11F*(stifcu22F + stifstxF))*sin(2*alfa[jj])))/2.) / 

	(2.*gamma*(2*gamma*(pow(ex - ey,2) + pow(gamma,2))*stifcu21F*pow(cos(alfa[jj]),2) - 

		   (-ex + ey)*gamma*R2*(2*gamma*(-stifcu21F + stifcu22F)*pow(cos(alfa[jj]),2) + 8*(Fcu1F - Fcu2F)*pow(cos(alfa[jj]),3)*sin(alfa[jj]) + 

					2*gamma*(-stifcu11F + stifcu12F)*pow(sin(alfa[jj]),2)) + 

		   (pow(ex - ey,2) + pow(gamma,2))*(8*(-Fcu1F + Fcu2F)*pow(cos(alfa[jj]),3)*sin(alfa[jj]) + 

						    2*gamma*(2*stifstxF + stifcu22F*pow(cos(alfa[jj]),2) + (stifcu11F + stifcu12F)*pow(sin(alfa[jj]),2)))));

      

      

      dSydGamma = (4*(-ex + ey)*(Fcu1F - Fcu2F)*(pow(ex - ey,2) + pow(gamma,2))*(stifcu11F + stifcu21F + stifstxF)*sin(2*alfa[jj]) + 

		   2*(-ex + ey)*(Fcu1F - Fcu2F)*(pow(ex - ey,2) + pow(gamma,2))*(stifcu11F + stifcu21F + stifstxF)*sin(4*alfa[jj]) + 

		   R2*(2*pow(gamma,4)*(-stifcu11F + stifcu12F - stifcu21F + stifcu22F)*stifstxF + 

		       2*pow(gamma,4)*(2*stifcu12F*stifcu21F - 2*stifcu11F*stifcu22F + 

				       (-stifcu11F + stifcu12F + stifcu21F - stifcu22F)*stifstxF)*cos(2*alfa[jj]) + 

		       2*(Fcu1F - Fcu2F)*gamma*((2*pow(ex - ey,2) + pow(gamma,2))*(stifcu11F + stifcu21F) - 

						pow(gamma,2)*(stifcu12F + stifcu22F) + 2*pow(ex - ey,2)*stifstxF)*sin(2*alfa[jj]) + 

		       (Fcu1F - Fcu2F)*gamma*((2*pow(ex - ey,2) + pow(gamma,2))*(stifcu11F + stifcu21F) - 

					      pow(gamma,2)*(stifcu12F + stifcu22F) + 2*pow(ex - ey,2)*stifstxF)*sin(4*alfa[jj])))/

	(2.*gamma*(2*gamma*(pow(ex - ey,2) + pow(gamma,2))*stifcu21F*pow(cos(alfa[jj]),2) - 

		   (-ex + ey)*gamma*R2*(2*gamma*(-stifcu21F + stifcu22F)*pow(cos(alfa[jj]),2) + 8*(Fcu1F - Fcu2F)*pow(cos(alfa[jj]),3)*sin(alfa[jj]) + 

					2*gamma*(-stifcu11F + stifcu12F)*pow(sin(alfa[jj]),2)) + 

		   (pow(ex - ey,2) + pow(gamma,2))*(8*(-Fcu1F + Fcu2F)*pow(cos(alfa[jj]),3)*sin(alfa[jj]) + 

						    2*gamma*(2*stifstxF + stifcu22F*pow(cos(alfa[jj]),2) + (stifcu11F + stifcu12F)*pow(sin(alfa[jj]),2)))));

      

      

      dSydEy = (-(gamma*(pow(ex - ey,2) + pow(gamma,2))*((stifcu11F + stifcu12F + stifcu21F + stifcu22F)*stifstyF + 

							 stifstxF*(stifcu11F + stifcu12F + stifcu21F + stifcu22F + 4*stifstyF))) + 

		gamma*(pow(ex - ey,2) + pow(gamma,2))*(-stifcu11F - stifcu12F + stifcu21F + stifcu22F)*(stifstxF - stifstyF)*cos(2*alfa[jj]) + 

		2*(Fcu1F - Fcu2F)*(pow(ex - ey,2) + pow(gamma,2))*(stifcu11F + stifcu12F + stifcu21F + stifcu22F + stifstxF + stifstyF)*

		sin(2*alfa[jj]) + (Fcu1F - Fcu2F)*(pow(ex - ey,2) + pow(gamma,2))*

		(stifcu11F + stifcu12F + stifcu21F + stifcu22F + stifstxF + stifstyF)*sin(4*alfa[jj]) + 

		(-ex + ey)*gamma*R2*(gamma*(stifcu11F - stifcu12F + stifcu21F - stifcu22F)*(stifstxF - stifstyF) + 

				     gamma*(-((stifcu21F - stifcu22F)*(stifstxF + stifstyF)) - stifcu12F*(4*stifcu21F + stifstxF + stifstyF) + 

					    stifcu11F*(4*stifcu22F + stifstxF + stifstyF))*cos(2*alfa[jj]) + 

				     2*(Fcu1F - Fcu2F)*(stifcu11F + stifcu12F + stifcu21F + stifcu22F + stifstxF + stifstyF)*sin(2*alfa[jj]) + 

				     (Fcu1F - Fcu2F)*(stifcu11F + stifcu12F + stifcu21F + stifcu22F + stifstxF + stifstyF)*sin(4*alfa[jj])))/

	(-2*gamma*(pow(ex - ey,2) + pow(gamma,2))*stifcu21F*pow(cos(alfa[jj]),2) + 

	 (-ex + ey)*gamma*R2*(2*gamma*(-stifcu21F + stifcu22F)*pow(cos(alfa[jj]),2) + 8*(Fcu1F - Fcu2F)*pow(cos(alfa[jj]),3)*sin(alfa[jj]) + 

			      2*gamma*(-stifcu11F + stifcu12F)*pow(sin(alfa[jj]),2)) + 

	 (-pow(-ex + ey,2) - pow(gamma,2))*(8*(-Fcu1F + Fcu2F)*pow(cos(alfa[jj]),3)*sin(alfa[jj]) + 

					    2*gamma*(2*stifstxF + stifcu22F*pow(cos(alfa[jj]),2) + (stifcu11F + stifcu12F)*pow(sin(alfa[jj]),2))));

      

    }

    

    kData[0] += dSydEy;

    kData[1] += dSydEy*StripCenterLoc(jj);

    kData[2] += dSydGamma;

    kData[3] += dSydEy*StripCenterLoc(jj);

    

    kData[4] += dSydEy*StripCenterLoc(jj)*StripCenterLoc(jj);

    kData[5] += dSydGamma*StripCenterLoc(jj);

    kData[6] += dTdEy;

    kData[7] += dTdEy*StripCenterLoc(jj);

    kData[8] += dTdGamma;

    

    sData[0] += Fy;

    sData[1] += Fy*StripCenterLoc(jj);

    sData[2] += Fxy;

    

    sy[jj] = Fy/StripLoc(jj,1); // avg stresses 

    txy[jj] = Fxy/StripLoc(jj,1);

  }



  return res;

}









const Matrix&

FiberSection2dInt::getInitialTangent(void)	

{

  kData[0] = 0.0; kData[1] = 0.0; kData[2] = 0.0; kData[3] = 0.0; kData[4] = 0.0; kData[5] = 0.0;

  kData[6] = 0.0; kData[7] = 0.0; kData[8] = 0.0;



  for (int i = 0; i < numFibers; i++) {

    

    double stifstyF, stifcu22F;

    stifstyF =0.0;

    stifcu22F =0.0;

    

    UniaxialMaterial *theMat1 = theMaterials1[i];  

    

    int tag=theMat1->getTag();

    double y = StripCenterLoc(FiberLoc(i));

    double A = matData[2*i+1];

    

    if (tag>1000){			// to distinguish concrete & steel tag>1000 => steel

      

      double tsy = theMat1->getInitialTangent();

      stifstyF =tsy*A;

    }

    else{			// concrete

      double tc1 = theMat1->getInitialTangent();

      stifcu22F =tc1*A;

    }

    

    double dTdEy = 0.0;

    double dTdGamma = stifcu22F/2.0;

    double dSydGamma = 0.0;

    double dSydEy = stifcu22F + stifstyF;

    

    kData[0] += dSydEy;

    kData[1] += dSydEy*y;

    kData[2] += dSydGamma;

    kData[3] += dSydEy*y;

    kData[4] += dSydEy*y*y;

    kData[5] += dSydGamma*y;

    kData[6] += dTdEy;

    kData[7] += dTdEy*y;

    kData[8] += dTdGamma;

  }

  

  return *ks;

}



const Matrix&

FiberSection2dInt::getSectionTangent(void)	

{

  

  return *ks;

}



const Vector&

FiberSection2dInt::getStressResultant(void)	

{

  return *s;

}



SectionForceDeformation*

FiberSection2dInt::getCopy(void)

{

  FiberSection2dInt *theCopy = new FiberSection2dInt ();

  theCopy->setTag(this->getTag());



  theCopy->numFibers = numFibers;



  if (numFibers != 0) {

    theCopy->theMaterials1 = new UniaxialMaterial *[numFibers];

	theCopy->theMaterials2 = new UniaxialMaterial *[numFibers];



    if (theCopy->theMaterials1 == 0) {

      opserr <<"FiberSection2dInt::getCopy -- failed to allocate Material pointers\n";

      exit(-1);

    }

  

    theCopy->matData = new double [numFibers*2];



    if (theCopy->matData == 0) {

      opserr << "FiberSection2dInt::getCopy -- failed to allocate double array for material data\n";

      exit(-1);

    }

			    

    for (int i = 0; i < numFibers; i++) {

      theCopy->matData[i*2] = matData[i*2];

      theCopy->matData[i*2+1] = matData[i*2+1];

      theCopy->theMaterials1[i] = theMaterials1[i]->getCopy();

	  theCopy->theMaterials2[i] = theMaterials2[i]->getCopy();

		

      if (theCopy->theMaterials1[i] == 0) {

	opserr <<"FiberSection2dInt::getCopy -- failed to get copy of a Material";

	exit(-1);

      }

    }  

  }





  theCopy->numHFibers = numHFibers;



  if (numHFibers != 0) {

    theCopy->theHMaterials = new UniaxialMaterial *[numHFibers * NStrip];



    if (theCopy->theHMaterials == 0) {

      opserr <<"FiberSection2dInt::getCopy -- failed to allocate HMaterial pointers\n";

      exit(-1);

    }

  

    theCopy->matHData = new double [numHFibers*2];



    if (theCopy->matHData == 0) {

      opserr << "FiberSection2dInt::getCopy -- failed to allocate double array for Hmaterial data\n";

      exit(-1);

    }

			    

    for (int i = 0; i < numHFibers; i++) {

      theCopy->matHData[i*2] = matHData[i*2];

      theCopy->matHData[i*2+1] = matHData[i*2+1];

      for (int jj = 0; jj < NStrip; jj++) {

	theCopy->theHMaterials[i * numHFibers + jj] = theHMaterials[i * numHFibers + jj]->getCopy();



	if (theCopy->theHMaterials[i * numHFibers + jj] == 0) {

	  opserr <<"FiberSection2dInt::getCopy -- failed to get copy of a HMaterial";

	  exit(-1);

	}

      }

    }  

  }



  theCopy->NStrip = NStrip;

  theCopy->NStrip1 = NStrip1;

  theCopy->NStrip2 = NStrip2;

  theCopy->NStrip3 = NStrip3;

  theCopy->tavg1 = tavg1;

  theCopy->tavg2 = tavg2;

  theCopy->tavg3 = tavg3;

  

for (int j = 0; j < NStrip; j++) {

  theCopy->sy[j] = sy[j];

  theCopy->txy[j] = txy[j];

  theCopy->alfa[j] = alfa[j]; 

  theCopy->alfaCommit[j] = alfaCommit[j];

  theCopy->iterOut[j] = iterOut[j];

  theCopy->iterCommit[j] = iterCommit[j];

  theCopy->exOut[j] = exOut[j];

  theCopy->exCommit[j] = exCommit[j];

  theCopy->eyCommit[j] = eyCommit[j];

  theCopy->e1Commit[j] = e1Commit[j];

  theCopy->e2Commit[j] = e2Commit[j];

  theCopy->sxCommit[j] = sxCommit[j];

  theCopy->syCommit[j] = syCommit[j];

  theCopy->s1Commit[j] = s1Commit[j];

  theCopy->s2Commit[j] = s2Commit[j];

}



  theCopy->StripCenterLoc = StripCenterLoc;

  theCopy->StripLoc = StripLoc;

  theCopy->FiberLoc = FiberLoc;



  theCopy->eCommit = eCommit;

  theCopy->e = e;

  theCopy->yBar = yBar;

  theCopy->ymin = ymin;

  theCopy->ymax = ymax;



  theCopy->kData[0] = kData[0];

  theCopy->kData[1] = kData[1];

  theCopy->kData[2] = kData[2];

  theCopy->kData[3] = kData[3];

  theCopy->kData[4] = kData[4];

  theCopy->kData[5] = kData[5];

  theCopy->kData[6] = kData[6];

  theCopy->kData[7] = kData[7];

  theCopy->kData[8] = kData[8];



  theCopy->sData[0] = sData[0];

  theCopy->sData[1] = sData[1];

  theCopy->sData[2] = sData[2];



  return theCopy;

}



const ID&

FiberSection2dInt::getType ()			

{

  return code;

}



int

FiberSection2dInt::getOrder () const	

{

  return 3;								

}



int

FiberSection2dInt::commitStateB(void)			

{

  int err = 0;



  for (int i = 0; i < numFibers; i++){

    err += theMaterials1[i]->commitState();

	err += theMaterials2[i]->commitState();

  }

  for (int H = 0; H < numHFibers; H++)

    for (int jj = 0; jj < NStrip; jj++) err += theHMaterials[H * numHFibers + jj]->commitState();



  eCommit = e;



  for (int jj = 0; jj < NStrip; jj++){

	  iterCommit[jj] = iterOut[jj];

	  alfaCommit[jj] = alfa[jj];

	  exCommit[jj] = exOut[jj];



  }



  return err;

}







int

FiberSection2dInt::revertToLastCommitB(double L)	

{

  int err = 0;



  // Last committed section deformations

  e = eCommit;

  for (int ii = 0; ii < NStrip; ii++) alfa[ii] = alfaCommit[ii];



  kData[0] = 0.0; kData[1] = 0.0; kData[2] = 0.0; kData[3] = 0.0; kData[4] = 0.0; kData[5] = 0.0;

  kData[6] = 0.0; kData[7] = 0.0; kData[8] = 0.0;



  sData[0] = 0.0; sData[1] = 0.0; sData[2] = 0.0;



  double gamma = e(2);						



for (int jj = 0; jj < NStrip; jj++) {

	double tavg;

	if (jj<NStrip1)

		tavg=tavg1;

	else {

		if (jj<NStrip2+NStrip1)

			tavg=tavg2;

		else

			tavg=tavg3;

	}



	double stifstxF, stifstyF, stifcu11F, stifcu12F, stifcu22F, stifcu21F, Fcu1F, Fcu2F;

	double Fy, Fxy;



	double ex=exCommit[jj];



	 Fy = 0.0;

	 Fxy = 0.0;



	 stifstxF =0.0;

	 stifstyF =0.0;

	 stifcu11F =0.0;

     stifcu12F =0.0;

     stifcu22F =0.0;

     stifcu21F =0.0;

     Fcu1F =0.0;

     Fcu2F =0.0;



	 double ey = e(0) + StripCenterLoc(jj)*e(1);

	 double e1;                          

	 double e2;

	 double root = sqrt(pow(-ex + ey,2) + pow(gamma,2));



	 if (fabs(alfa[jj] - PI/2) <= DBL_EPSILON) {

		 e1 = ex; 

		 e2 = ey; 

	 }

	 else {

		 e1 = ey - gamma/2*tan(alfa[jj]);                          

		 e2 = ex + gamma/2*tan(alfa[jj]);

	 }





	for (int i = 2; i <= StripLoc(jj,0)+1; i++) {

    

	  int fibNum=StripLoc(jj,i);

	  UniaxialMaterial *theMat1 = theMaterials1[fibNum];  

	  UniaxialMaterial *theMat2 = theMaterials2[fibNum];  

	  

	  int tag=theMat1->getTag();

	  

	  if (tag>1000){			// to distinguish concrete & steel tag>1000 => steel

	    double y = matData[2*fibNum] - yBar;

	    double A = matData[2*fibNum+1];	

	    err += theMat1->revertToLastCommit();

	    double tsy = theMat1->getTangent();

	    double ssy = theMat1->getStress();

	    

	    Fy += ssy*A;

	    Fxy += 0.0;

	    

	    stifstyF +=tsy*A;

	    

	  }

	  else{			// concrete

	    double y = matData[2*fibNum] - yBar;

	    double A = matData[2*fibNum+1];	

	    

	    err += theMat1->revertToLastCommit();

	    double tc1 = theMat1->getTangent();

	    double sc1 = theMat1->getStress();

	    

	    err += theMat2->revertToLastCommit();

	    double tc2 = theMat2->getTangent();

	    double sc2 = theMat2->getStress();

	    

	    double tc12=0.0;

	    double tc21=0.0;

	    double beta1, beta2;



	    double e0 = 0.0;

	    const char *theData = "ec";

	    static Information theInfo;

	    if (theMat1->getVariable(theData, theInfo) == 0)

	      e0 = theInfo.theDouble;

	    

	    this -> beta(e0, e2, sc1, tc1, tc12, beta1);

	    this -> beta(e0, e1, sc2, tc2, tc21, beta2);

	    

	    Fy += ((sc1+sc2)*0.5+(sc1-sc2)*0.5*cos(2.0*alfa[jj]))*A;

	    Fxy += -((sc1-sc2)*0.5*sin(2.0*alfa[jj]))*A;

	    

	    stifcu11F +=tc1*A;

	    stifcu12F +=tc12*A; 

	    stifcu22F +=tc2*A;

	    stifcu21F +=tc21*A; 

	    

	    Fcu1F += sc1*A;

	    Fcu2F += sc2*A;

            

	  }

	}

	int Hloc = 0;

	for (int H = 0; H < numHFibers; H++) {				

	  UniaxialMaterial *theHMat = theHMaterials[H * numHFibers + jj];

	  double Hy = matHData[Hloc++];

	  double HA = matHData[Hloc++];

	  

	  err += theHMat->revertToLastCommit();

	  double Ht = theHMat->getTangent();

	  double Hs = theHMat->getStress();

	  

	  stifstxF += Ht*HA*StripLoc(jj,1)/(L*tavg);

	  

	}

	

	

	double dTdEy;

	double dTdGamma;

	double dSydGamma;

	double dSydEy;

	

	if ((fabs(gamma) <= DBL_EPSILON)||(fabs(alfa[jj]) <= DBL_EPSILON)||(fabs(alfa[jj]-PI/2.0) <= DBL_EPSILON)) {

	  dTdEy = 0;

	  dTdGamma = (-Fcu1F + Fcu2F)/(2.*sqrt(pow(ex - ey,2)));

	  dSydGamma = 0;

	  dSydEy = stifcu22F - (stifcu12F*stifcu21F)/(stifcu11F + stifstxF) + stifstyF;

	}

	else {

	  

	  double R2=sqrt(1 + pow(ex - ey,2)/pow(gamma,2));

	  

	  dTdEy = (-((Fcu1F - Fcu2F)*(pow(ex,2) - 2*ex*ey + pow(ey,2) + pow(gamma,2))*

		     (stifcu11F + stifcu12F - 3*stifcu21F - 3*stifcu22F - 2*stifstxF)) + 

		   4*(Fcu1F - Fcu2F)*(pow(ex - ey,2) + pow(gamma,2))*(stifcu21F + stifcu22F + stifstxF)*cos(2*alfa[jj]) + 

		   (Fcu1F - Fcu2F)*(pow(ex - ey,2) + pow(gamma,2))*(stifcu11F + stifcu12F + stifcu21F + stifcu22F + 2*stifstxF)*cos(4*alfa[jj]) + 

		   2*gamma*(pow(ex - ey,2) + pow(gamma,2))*(stifcu11F + stifcu12F - stifcu21F - stifcu22F)*stifstxF*sin(2*alfa[jj]) + 

		   (-ex + ey)*gamma*R2*((-Fcu1F + Fcu2F)*(stifcu11F + stifcu12F - 3*stifcu21F - 3*stifcu22F - 2*stifstxF) + 

					4*(Fcu1F - Fcu2F)*(stifcu21F + stifcu22F + stifstxF)*cos(2*alfa[jj]) + 

					(Fcu1F - Fcu2F)*(stifcu11F + stifcu12F + stifcu21F + stifcu22F + 2*stifstxF)*cos(4*alfa[jj]) + 

					2*gamma*(2*stifcu12F*stifcu21F - 2*stifcu11F*stifcu22F + (-stifcu11F + stifcu12F + stifcu21F - stifcu22F)*stifstxF)*

					sin(2*alfa[jj])))/(2.*(-2*gamma*(pow(ex - ey,2) + pow(gamma,2))*stifcu21F*pow(cos(alfa[jj]),2) + 

							       (-ex + ey)*gamma*R2*(2*gamma*(-stifcu21F + stifcu22F)*pow(cos(alfa[jj]),2) + 8*(Fcu1F - Fcu2F)*pow(cos(alfa[jj]),3)*sin(alfa[jj]) + 

										    2*gamma*(-stifcu11F + stifcu12F)*pow(sin(alfa[jj]),2)) + 

							       (-pow(-ex + ey,2) - pow(gamma,2))*(8*(-Fcu1F + Fcu2F)*pow(cos(alfa[jj]),3)*sin(alfa[jj]) + 

												  2*gamma*(2*stifstxF + stifcu22F*pow(cos(alfa[jj]),2) + (stifcu11F + stifcu12F)*pow(sin(alfa[jj]),2)))));

	  

	  

	  dTdGamma = ((-ex + ey)*(Fcu1F - Fcu2F)*(pow(ex - ey,2) + pow(gamma,2))*(-stifcu11F + 3*stifcu21F + 2*stifstxF) + 

		      4*(-ex + ey)*(Fcu1F - Fcu2F)*(pow(ex - ey,2) + pow(gamma,2))*(stifcu21F + stifstxF)*cos(2*alfa[jj]) + 

		      (-ex + ey)*(Fcu1F - Fcu2F)*(pow(ex - ey,2) + pow(gamma,2))*(stifcu11F + stifcu21F + 2*stifstxF)*cos(4*alfa[jj]) + 

		      (gamma*R2*((Fcu1F - Fcu2F)*(-((2*pow(ex - ey,2) + pow(gamma,2))*(stifcu11F - 3*stifcu21F)) + 

						  pow(gamma,2)*(stifcu12F - 3*stifcu22F) + 4*pow(ex - ey,2)*stifstxF) + 

				 4*(Fcu1F - Fcu2F)*(pow(gamma,2)*(stifcu21F - stifcu22F) + 2*pow(ex,2)*(stifcu21F + stifstxF) - 

						    4*ex*ey*(stifcu21F + stifstxF) + 2*pow(ey,2)*(stifcu21F + stifstxF))*cos(2*alfa[jj]) + 

				 (Fcu1F - Fcu2F)*((2*pow(ex - ey,2) + pow(gamma,2))*(stifcu11F + stifcu21F) - pow(gamma,2)*(stifcu12F + stifcu22F) + 

						  4*pow(ex - ey,2)*stifstxF)*cos(4*alfa[jj]) + 

				 4*pow(gamma,3)*((-stifcu21F + stifcu22F)*stifstxF - stifcu12F*(stifcu21F + stifstxF) + 

						 stifcu11F*(stifcu22F + stifstxF))*sin(2*alfa[jj])))/2.)/

	    (2.*gamma*(2*gamma*(pow(ex - ey,2) + pow(gamma,2))*stifcu21F*pow(cos(alfa[jj]),2) - 

		       (-ex + ey)*gamma*R2*(2*gamma*(-stifcu21F + stifcu22F)*pow(cos(alfa[jj]),2) + 8*(Fcu1F - Fcu2F)*pow(cos(alfa[jj]),3)*sin(alfa[jj]) + 

					    2*gamma*(-stifcu11F + stifcu12F)*pow(sin(alfa[jj]),2)) + 

		       (pow(ex - ey,2) + pow(gamma,2))*(8*(-Fcu1F + Fcu2F)*pow(cos(alfa[jj]),3)*sin(alfa[jj]) + 

							2*gamma*(2*stifstxF + stifcu22F*pow(cos(alfa[jj]),2) + (stifcu11F + stifcu12F)*pow(sin(alfa[jj]),2)))));

	  

	  

	  dSydGamma = (4*(-ex + ey)*(Fcu1F - Fcu2F)*(pow(ex - ey,2) + pow(gamma,2))*(stifcu11F + stifcu21F + stifstxF)*sin(2*alfa[jj]) + 

		       2*(-ex + ey)*(Fcu1F - Fcu2F)*(pow(ex - ey,2) + pow(gamma,2))*(stifcu11F + stifcu21F + stifstxF)*sin(4*alfa[jj]) + 

		       R2*(2*pow(gamma,4)*(-stifcu11F + stifcu12F - stifcu21F + stifcu22F)*stifstxF + 

			   2*pow(gamma,4)*(2*stifcu12F*stifcu21F - 2*stifcu11F*stifcu22F + 

					   (-stifcu11F + stifcu12F + stifcu21F - stifcu22F)*stifstxF)*cos(2*alfa[jj]) + 

			   2*(Fcu1F - Fcu2F)*gamma*((2*pow(ex - ey,2) + pow(gamma,2))*(stifcu11F + stifcu21F) - 

						    pow(gamma,2)*(stifcu12F + stifcu22F) + 2*pow(ex - ey,2)*stifstxF)*sin(2*alfa[jj]) + 

			   (Fcu1F - Fcu2F)*gamma*((2*pow(ex - ey,2) + pow(gamma,2))*(stifcu11F + stifcu21F) - 

						  pow(gamma,2)*(stifcu12F + stifcu22F) + 2*pow(ex - ey,2)*stifstxF)*sin(4*alfa[jj])))/

	    (2.*gamma*(2*gamma*(pow(ex - ey,2) + pow(gamma,2))*stifcu21F*pow(cos(alfa[jj]),2) - 

		       (-ex + ey)*gamma*R2*(2*gamma*(-stifcu21F + stifcu22F)*pow(cos(alfa[jj]),2) + 8*(Fcu1F - Fcu2F)*pow(cos(alfa[jj]),3)*sin(alfa[jj]) + 

					    2*gamma*(-stifcu11F + stifcu12F)*pow(sin(alfa[jj]),2)) + 

		       (pow(ex - ey,2) + pow(gamma,2))*(8*(-Fcu1F + Fcu2F)*pow(cos(alfa[jj]),3)*sin(alfa[jj]) + 

							2*gamma*(2*stifstxF + stifcu22F*pow(cos(alfa[jj]),2) + (stifcu11F + stifcu12F)*pow(sin(alfa[jj]),2)))));

	  

	  

	  dSydEy = (-(gamma*(pow(ex - ey,2) + pow(gamma,2))*((stifcu11F + stifcu12F + stifcu21F + stifcu22F)*stifstyF + 

							     stifstxF*(stifcu11F + stifcu12F + stifcu21F + stifcu22F + 4*stifstyF))) + 

		    gamma*(pow(ex - ey,2) + pow(gamma,2))*(-stifcu11F - stifcu12F + stifcu21F + stifcu22F)*(stifstxF - stifstyF)*cos(2*alfa[jj]) + 

		    2*(Fcu1F - Fcu2F)*(pow(ex - ey,2) + pow(gamma,2))*(stifcu11F + stifcu12F + stifcu21F + stifcu22F + stifstxF + stifstyF)*

		    sin(2*alfa[jj]) + (Fcu1F - Fcu2F)*(pow(ex - ey,2) + pow(gamma,2))*

		    (stifcu11F + stifcu12F + stifcu21F + stifcu22F + stifstxF + stifstyF)*sin(4*alfa[jj]) + 

		    (-ex + ey)*gamma*R2*(gamma*(stifcu11F - stifcu12F + stifcu21F - stifcu22F)*(stifstxF - stifstyF) + 

					 gamma*(-((stifcu21F - stifcu22F)*(stifstxF + stifstyF)) - stifcu12F*(4*stifcu21F + stifstxF + stifstyF) + 

						stifcu11F*(4*stifcu22F + stifstxF + stifstyF))*cos(2*alfa[jj]) + 

					 2*(Fcu1F - Fcu2F)*(stifcu11F + stifcu12F + stifcu21F + stifcu22F + stifstxF + stifstyF)*sin(2*alfa[jj]) + 

					 (Fcu1F - Fcu2F)*(stifcu11F + stifcu12F + stifcu21F + stifcu22F + stifstxF + stifstyF)*sin(4*alfa[jj])))/

	    (-2*gamma*(pow(ex - ey,2) + pow(gamma,2))*stifcu21F*pow(cos(alfa[jj]),2) + 

	     (-ex + ey)*gamma*R2*(2*gamma*(-stifcu21F + stifcu22F)*pow(cos(alfa[jj]),2) + 8*(Fcu1F - Fcu2F)*pow(cos(alfa[jj]),3)*sin(alfa[jj]) + 

				  2*gamma*(-stifcu11F + stifcu12F)*pow(sin(alfa[jj]),2)) + 

	     (-pow(-ex + ey,2) - pow(gamma,2))*(8*(-Fcu1F + Fcu2F)*pow(cos(alfa[jj]),3)*sin(alfa[jj]) + 

						2*gamma*(2*stifstxF + stifcu22F*pow(cos(alfa[jj]),2) + (stifcu11F + stifcu12F)*pow(sin(alfa[jj]),2))));

	  

	  

	}

	

	kData[0] += dSydEy;

	kData[1] += dSydEy*StripCenterLoc(jj);

	kData[2] += dSydGamma;

	kData[3] += dSydEy*StripCenterLoc(jj);

    

	kData[4] += dSydEy*StripCenterLoc(jj)*StripCenterLoc(jj);

	kData[5] += dSydGamma*StripCenterLoc(jj);

	kData[6] += dTdEy;

	kData[7] += dTdEy*StripCenterLoc(jj);

	kData[8] += dTdGamma;

	

	sData[0] += Fy;

	sData[1] += Fy*StripCenterLoc(jj);

	sData[2] += Fxy;

	

	sy[jj] = Fy/StripLoc(jj,1); // avg stresses 

	txy[jj] = Fxy/StripLoc(jj,1);

	

 }

 

  return err;

}







int

FiberSection2dInt::revertToStartB(void)		

{

  // revert the fibers to start    

  int err = 0;



  kData[0] = 0.0; kData[1] = 0.0; kData[2] = 0.0; kData[3] = 0.0; kData[4] = 0.0; kData[5] = 0.0;

  kData[6] = 0.0; kData[7] = 0.0; kData[8] = 0.0;



  sData[0] = 0.0; sData[1] = 0.0;	sData[2] = 0.0;

  

	

  for (int i = 0; i < numFibers; i++) {



		double stifstyF, stifcu22F;

		stifstyF =0.0;

		stifcu22F =0.0;

		UniaxialMaterial *theMat1 = theMaterials1[i];  



		int tag=theMat1->getTag();

		double y = StripCenterLoc(FiberLoc(i));

		double A = matData[2*i+1];	



		if (tag>1000){			// to distinguish concrete & steel tag>1000 => steel

			err += theMat1->revertToStart();



			// get material stress & tangent for this strain and determine ks and fs

			double tsy = theMat1->getTangent();

			double ssy = theMat1->getStress();

			stifstyF =tsy*A;

		}

		else{			// concrete

			err += theMat1->revertToStart();

			double tc1 = theMat1->getTangent();

			stifcu22F =tc1*A;

		}



		double dTdEy = 0.0;

		double dTdGamma = stifcu22F/2.0;

		double dSydGamma = 0.0;

		double dSydEy = stifcu22F + stifstyF;



		kData[0] += dSydEy;

		kData[1] += dSydEy*y;

		kData[2] += dSydGamma;

		kData[3] += dSydEy*y;

		kData[4] += dSydEy*y*y;

		kData[5] += dSydGamma*y;

		kData[6] += dTdEy;

		kData[7] += dTdEy*y;

		kData[8] += dTdGamma;



		sData[0] += 0.0; 

		sData[1] += 0.0; 

		sData[2] += 0.0; 

	}



	for (int jj = 0; jj < NStrip; jj++) {

		int Hloc = 0;

		for (int H = 0; H < numHFibers; H++) {				

			UniaxialMaterial *theHMat = theHMaterials[H * numHFibers + jj];

			double Hy = matHData[Hloc++];	

			double HA = matHData[Hloc++];

			err += theHMat->revertToStart();

			double Ht = theHMat->getTangent();

		}

	}

  return err;

}



int

FiberSection2dInt::sendSelf(int commitTag, Channel &theChannel)	

{

  int res = 0;



  // create an id to send objects tag and numFibers, 

  //     size 3 so no conflict with matData below if just 1 fiber

  static ID data(3);

  data(0) = this->getTag();

  data(1) = numFibers;

  int dbTag = this->getDbTag();

  res += theChannel.sendID(dbTag, commitTag, data);

  if (res < 0) {

    opserr <<  "FiberSection2dInt::sendSelf - failed to send ID data\n";

    return res;

  }    



  if (numFibers != 0) {

    // create an id containingg classTag and dbTag for each material & send it

    ID materialData(2*numFibers);

    for (int i=0; i<numFibers; i++) {

      UniaxialMaterial *theMat1 = theMaterials1[i];

      UniaxialMaterial *theMat2 = theMaterials2[i];



      materialData(2*i) = theMat1->getClassTag();

      int matDbTag = theMat1->getDbTag();

      if (matDbTag == 0) {

	matDbTag = theChannel.getDbTag();

	if (matDbTag != 0)

	  theMat1->setDbTag(matDbTag);

      }

      materialData(2*i+1) = matDbTag;

    }    

    

    res += theChannel.sendID(dbTag, commitTag, materialData);

    if (res < 0) {

      opserr <<  "FiberSection2dInt::sendSelf - failed to send material data\n";

      return res;

    }    



    // send the fiber data, i.e. area and loc

    Vector fiberData(matData, 2*numFibers);

    res += theChannel.sendVector(dbTag, commitTag, fiberData);

    if (res < 0) {

      opserr <<  "FiberSection2dInt::sendSelf - failed to send material data\n";

      return res;

    }    



    // now invoke send(0 on all the materials

    for (int j=0; j<numFibers; j++){

      theMaterials1[j]->sendSelf(commitTag, theChannel);

	  theMaterials2[j]->sendSelf(commitTag, theChannel);

	}

  }

  return res;

}



int

FiberSection2dInt::recvSelf(int commitTag, Channel &theChannel,

			 FEM_ObjectBroker &theBroker)						

{

  int res = 0;

  static ID data(3);

  

  int dbTag = this->getDbTag();

  res += theChannel.recvID(dbTag, commitTag, data);

  if (res < 0) {

    opserr <<  "FiberSection2dInt::recvSelf - failed to recv ID data\n";

    return res;

  }    

  this->setTag(data(0));



  // recv data about materials objects, classTag and dbTag

  if (data(1) != 0) {

    ID materialData(2*data(1));

    res += theChannel.recvID(dbTag, commitTag, materialData);

    if (res < 0) {

      opserr <<  "FiberSection2dInt::recvSelf - failed to recv material data\n";

      return res;

    }    



    // if current arrays not of correct size, release old and resize

    if (theMaterials1 == 0 || numFibers != data(1)) {

      // delete old stuff if outa date

      if (theMaterials1 != 0) {

	for (int i=0; i<numFibers; i++){

	  delete theMaterials1[i];

	  delete theMaterials2[i];

	}

	delete [] theMaterials1;

	delete [] theMaterials2;

	if (matData != 0)

	  delete [] matData;

	matData = 0;

	theMaterials1 = 0;

	theMaterials2 = 0;

      }



      // create memory to hold material pointers and fiber data

      numFibers = data(1);

      if (numFibers != 0) {

	theMaterials1 = new UniaxialMaterial *[numFibers];

	theMaterials2 = new UniaxialMaterial *[numFibers];



	if (theMaterials1 == 0) {

	  opserr <<"FiberSection2dInt::recvSelf -- failed to allocate Material pointers\n";

	  exit(-1);

	}

	

	for (int j=0; j<numFibers; j++){

	  theMaterials1[j] = 0;

	  theMaterials2[j] = 0;

	}

	matData = new double [numFibers*2];



	if (matData == 0) {

	  opserr <<"FiberSection2dInt::recvSelf  -- failed to allocate double array for material data\n";

	  exit(-1);

	}

      }

    }



    Vector fiberData(matData, 2*numFibers);

    res += theChannel.recvVector(dbTag, commitTag, fiberData);

    if (res < 0) {

      opserr <<  "FiberSection2dInt::recvSelf - failed to recv material data\n";

      return res;

    }    



    int i;

    for (i=0; i<numFibers; i++) {

      int classTag = materialData(2*i);

      int dbTag = materialData(2*i+1);



      // if material pointed to is blank or not of corrcet type, 

      // release old and create a new one

      if (theMaterials1[i] == 0){

	theMaterials1[i] = theBroker.getNewUniaxialMaterial(classTag);

	theMaterials2[i] = theBroker.getNewUniaxialMaterial(classTag);



	  }

      else if (theMaterials1[i]->getClassTag() != classTag) {

	delete theMaterials1[i];

	theMaterials1[i] = theBroker.getNewUniaxialMaterial(classTag);    

	delete theMaterials2[i];

	theMaterials2[i] = theBroker.getNewUniaxialMaterial(classTag);      

	

      }



      if (theMaterials1[i] == 0) {

	opserr <<"FiberSection2dInt::recvSelf -- failed to allocate double array for material data\n";

	exit(-1);

      }



      theMaterials1[i]->setDbTag(dbTag);

      res += theMaterials1[i]->recvSelf(commitTag, theChannel, theBroker);

	  theMaterials2[i]->setDbTag(dbTag);

      res += theMaterials2[i]->recvSelf(commitTag, theChannel, theBroker);

    }



    double Qz = 0.0;

    double A  = 0.0;

    double yLoc, Area;



    // Recompute centroid

    for (i = 0; i < numFibers; i++) {

      yLoc = -matData[2*i];

      Area = matData[2*i+1];

      A  += Area;

      Qz += yLoc*Area;

    }

    

    yBar = -Qz/A;

  }    



  return res;

}



void

FiberSection2dInt::Print(OPS_Stream &s, int flag)

{

  s << "\nFiberSection2d, tag: " << this->getTag() << endln;

  s << "\tSection code: " << code;

  s << "\tNumber of Fibers: " << numFibers << endln;

  s << "\tCentroid: " << -yBar << endln;



  if (flag == 1) {

    int loc = 0;

    for (int i = 0; i < numFibers; i++) {

      s << "\nLocation (y) = (" << -matData[loc++] << ")";

      s << "\nArea = " << matData[loc++] << endln;

      theMaterials1[i]->Print(s, flag);

	  theMaterials2[i]->Print(s, flag);

    }

  }

}



Response*
FiberSection2dInt::setResponse(const char **argv, int argc, OPS_Stream &output)	
{
  const ID &type = this->getType();
  int typeSize = this->getOrder();
  
  Response *theResponse =0;
  
  output.tag("SectionOutput");
  output.attr("secType", this->getClassType());
  output.attr("secTag", this->getTag());

  // deformations
  if (strcmp(argv[0],"deformations") == 0 || strcmp(argv[0],"deformation") == 0) {
    for (int i=0; i<typeSize; i++) {
      int code = type(i);
      switch (code){
      case SECTION_RESPONSE_MZ:
	output.tag("ResponseType","kappaZ");
	break;
      case SECTION_RESPONSE_P:
	output.tag("ResponseType","eps");
	break;
      case SECTION_RESPONSE_VY:
	output.tag("ResponseType","gammaY");
	break;
      case SECTION_RESPONSE_MY:
	output.tag("ResponseType","kappaY");
	break;
      case SECTION_RESPONSE_VZ:
	output.tag("ResponseType","gammaZ");
	break;
      case SECTION_RESPONSE_T:
	output.tag("ResponseType","theta");
	break;
      default:
	output.tag("ResponseType","Unknown");
      }
    }
    theResponse =  new MaterialResponse(this, 1, this->getSectionDeformation());
    return theResponse;
  
  // forces
  } else if (strcmp(argv[0],"forces") == 0 || strcmp(argv[0],"force") == 0) {
    for (int i=0; i<typeSize; i++) {
      int code = type(i);
      switch (code){
      case SECTION_RESPONSE_MZ:
	output.tag("ResponseType","Mz");
	break;
      case SECTION_RESPONSE_P:
	output.tag("ResponseType","P");
	break;
      case SECTION_RESPONSE_VY:
	output.tag("ResponseType","Vy");
	break;
      case SECTION_RESPONSE_MY:
	output.tag("ResponseType","My");
	break;
      case SECTION_RESPONSE_VZ:
	output.tag("ResponseType","Vz");
	break;
      case SECTION_RESPONSE_T:
	output.tag("ResponseType","T");
	break;
      default:
	output.tag("ResponseType","Unknown");
      }
    }
    theResponse =  new MaterialResponse(this, 2, this->getStressResultant());
    return theResponse;
  
  // force and deformation
  } else if (strcmp(argv[0],"forceAndDeformation") == 0) { 
    for (int j=0; j<typeSize; j++) {
      int code = type(j);
      switch (code){
      case SECTION_RESPONSE_MZ:
	output.tag("ResponseType","kappaZ");
	break;
      case SECTION_RESPONSE_P:
	output.tag("ResponseType","eps");
	break;
      case SECTION_RESPONSE_VY:
	output.tag("ResponseType","gammaY");
	break;
      case SECTION_RESPONSE_MY:
	output.tag("ResponseType","kappaY");
	break;
      case SECTION_RESPONSE_VZ:
	output.tag("ResponseType","gammaZ");
	break;
      case SECTION_RESPONSE_T:
	output.tag("ResponseType","theta");
	break;
      default:
	output.tag("ResponseType","Unknown");
      }
    }
    for (int i=0; i<typeSize; i++) {
      int code = type(i);
      switch (code){
      case SECTION_RESPONSE_MZ:
	output.tag("ResponseType","Mz");
	break;
      case SECTION_RESPONSE_P:
	output.tag("ResponseType","P");
	break;
      case SECTION_RESPONSE_VY:
	output.tag("ResponseType","Vy");
	break;
      case SECTION_RESPONSE_MY:
	output.tag("ResponseType","My");
	break;
      case SECTION_RESPONSE_VZ:
	output.tag("ResponseType","Vz");
	break;
      case SECTION_RESPONSE_T:
	output.tag("ResponseType","T");
	break;
      default:
	output.tag("ResponseType","Unknown");
      }
    }

    theResponse =  new MaterialResponse(this, 4, Vector(2*this->getOrder()));
    return theResponse;
  }  

  // strip sigma y
  else if (strcmp(argv[0],"sigmaY") == 0)
    return new MaterialResponse(this, 105, this->getSigmaY());

  // strip Tau
  else if (strcmp(argv[0],"tau") == 0)
    return new MaterialResponse(this, 106, this->getTau());

  // strip Alpha 
  else if (strcmp(argv[0],"alpha") == 0)
    return new MaterialResponse(this, 107, this->getAlpha());

  // strip Alpha iter
  else if (strcmp(argv[0],"iter") == 0)
    return new MaterialResponse(this, 108, this->getIter());

  // strip eX
  else if (strcmp(argv[0],"eX") == 0)
    return new MaterialResponse(this, 109, this->getEX());

  // strip ey
  else if (strcmp(argv[0],"eY") == 0)
    return new MaterialResponse(this, 110, this->getEY());

  // strip e1
  else if (strcmp(argv[0],"e1") == 0)
    return new MaterialResponse(this, 111, this->getE1());

  // strip e2 
  else if (strcmp(argv[0],"e2") == 0)
    return new MaterialResponse(this, 112, this->getE2());

  // strip sX
  else if (strcmp(argv[0],"sX") == 0)
    return new MaterialResponse(this, 113, this->getSX());

  // strip sy
  else if (strcmp(argv[0],"sY") == 0)
    return new MaterialResponse(this, 114, this->getSY());

  // strip s1
  else if (strcmp(argv[0],"s1") == 0)
    return new MaterialResponse(this, 115, this->getS1());

  // strip s2                                                                                     
  else if (strcmp(argv[0],"s2") == 0)
    return new MaterialResponse(this, 116, this->getS2());

  // Check if fiber response is requested
  else if ((strcmp(argv[0],"fiber") == 0) || (strcmp(argv[0],"fiber1") == 0)) {	
    int key = numFibers;
    int passarg = 2;

    if (argc <= 2)          
      return 0;
    if (argc <= 3) {		  
      key = atoi(argv[1]);
      
      if (key < numFibers)
	return theMaterials1[key]->setResponse(&argv[passarg],argc-passarg,output);
      else
         return 0;
	}

    if (argc > 4) {         // find fiber closest to coord. with mat tag
      int matTag = atoi(argv[3]);
      double yCoord = atof(argv[1]);
      double closestDist;
      double ySearch, dy;
      double distance;
      int j;
      // Find first fiber with specified material tag
      for (j = 0; j < numFibers; j++) {
	if (matTag == theMaterials1[j]->getTag()) {
	  ySearch = -matData[2*j];
	  dy = ySearch-yCoord;
	  closestDist = fabs(dy);
	  key = j;
	  break;
	}
      }
      // Search the remaining fibers
      for ( ; j < numFibers; j++) {
	if (matTag == theMaterials1[j]->getTag()) {
	  ySearch = -matData[2*j];
	  dy = ySearch-yCoord;
	  distance = fabs(dy);
	  if (distance < closestDist) {
	    closestDist = distance;
	    key = j;
	  }
	}
      }
      passarg = 4;
    }

    else {                  // fiber near-to coordinate specified
      double yCoord = atof(argv[1]);
      double closestDist;
      double ySearch, dy;
      double distance;
      ySearch = -matData[0];
      dy = ySearch-yCoord;
      closestDist = fabs(dy);
      key = 0;
      for (int j = 1; j < numFibers; j++) {
	ySearch = -matData[2*j];
	dy = ySearch-yCoord;
	distance = fabs(dy);
	if (distance < closestDist) {
	  closestDist = distance;
	  key = j;
	}
      }
      passarg = 3;
    }
    
    if (key < numFibers)
      return theMaterials1[key]->setResponse(&argv[passarg],argc-passarg,output);
    else
      return 0;
  }

  

  else if (strcmp(argv[0],"fiber2") == 0) {
    int key = numFibers;
    int passarg = 2;
    
    if (argc <= 2)          
      return 0;

	if (argc <= 3) {		  
      key = atoi(argv[1]);
      if (key < numFibers)
         return theMaterials2[key]->setResponse(&argv[passarg],argc-passarg,output);
      else
         return 0;
	}

    if (argc > 4) {         // find fiber closest to coord. with mat tag
      int matTag = atoi(argv[3]);
      double yCoord = atof(argv[1]);
      double closestDist;
      double ySearch, dy;
      double distance;
      int j;
      // Find first fiber with specified material tag
      for (j = 0; j < numFibers; j++) {
	if (matTag == theMaterials2[j]->getTag()) {
	  ySearch = -matData[2*j];
	  dy = ySearch-yCoord;
	  closestDist = fabs(dy);
	  key = j;
	  break;
	}
      }
      // Search the remaining fibers
      for ( ; j < numFibers; j++) {
	if (matTag == theMaterials2[j]->getTag()) {
	  ySearch = -matData[2*j];
	  dy = ySearch-yCoord;
	  distance = fabs(dy);
	  if (distance < closestDist) {
	    closestDist = distance;
	    key = j;
	  }
	}
      }
      passarg = 4;
    }

    else {                  // fiber near-to coordinate specified
      double yCoord = atof(argv[1]);
      double closestDist;
      double ySearch, dy;
      double distance;
      ySearch = -matData[0];
      dy = ySearch-yCoord;
      closestDist = fabs(dy);
      key = 0;
      for (int j = 1; j < numFibers; j++) {
	ySearch = -matData[2*j];
	dy = ySearch-yCoord;
	distance = fabs(dy);
	if (distance < closestDist) {
	  closestDist = distance;
	  key = j;
	}
      }
      passarg = 3;
    }
    
    if (key < numFibers)
      return theMaterials2[key]->setResponse(&argv[passarg],argc-passarg,output);
    else
      return 0;
  }

  


  // Check if fiber response is requested
  else if (strcmp(argv[0],"Hfiber") == 0) {
    int HFibOut=atoi(argv[1])-1;
	int key = numHFibers;
    int passarg = 3;
    
    if (argc <= 3)          
      return 0;

	if (argc <= 4) {		  
      key = atoi(argv[2]);
      if (key < numHFibers)
         return theHMaterials[key*numHFibers + HFibOut]->setResponse(&argv[passarg],argc-passarg,output);
      else
         return 0;
	}

    if (argc > 5) {         // find fiber closest to coord. with mat tag
      int matTag = atoi(argv[4]);
      double yCoord = atof(argv[2]);
      double closestDist;
      double ySearch, dy;
      double distance;
      int j;
      // Find first fiber with specified material tag
      for (j = 0; j < numHFibers; j++) {
	if (matTag == theHMaterials[j * numHFibers + HFibOut]->getTag()) {
	  ySearch = -matHData[2*j];
	  dy = ySearch-yCoord;
	  closestDist = fabs(dy);
	  key = j;
	  break;
	}
      }
      // Search the remaining fibers
      for ( ; j < numHFibers; j++) {
	if (matTag == theHMaterials[j * numHFibers + HFibOut]->getTag()) {
	  ySearch = -matHData[2*j];
	  dy = ySearch-yCoord;
	  distance = fabs(dy);
	  if (distance < closestDist) {
	    closestDist = distance;
	    key = j;
	  }
	}
      }
      passarg = 5;
    }

    else {                  // fiber near-to coordinate specified
      double yCoord = atof(argv[2]);
      double closestDist;
      double ySearch, dy;
      double distance;
      ySearch = -matHData[0];
      dy = ySearch-yCoord;
      closestDist = fabs(dy);
      key = 0;
      for (int j = 1; j < numHFibers; j++) {
	ySearch = -matHData[2*j];
	dy = ySearch-yCoord;
	distance = fabs(dy);
	if (distance < closestDist) {
	  closestDist = distance;
	  key = j;
	}
      }
      passarg = 4;
    }
    
    if (key < numHFibers)
      return theHMaterials[key * numHFibers + HFibOut]->setResponse(&argv[passarg],argc-passarg,output);
    else
      return 0;
  }
 
  else
    return 0;
}





int 
FiberSection2dInt::getResponse(int responseID, Information &sectInfo)
{
  switch (responseID) {
    
  case 1:
    return sectInfo.setVector(this->getSectionDeformation());
    
  case 2:
    return sectInfo.setVector(this->getStressResultant());
    
  case 3:
    return sectInfo.setMatrix(this->getSectionTangent());
    
  case 4: {
    Vector &theVec = *(sectInfo.theVector);
    const Vector &e = this->getSectionDeformation();
    const Vector &s = this->getStressResultant();
    int order = this->getOrder();
    for (int i = 0; i < order; i++) {
      theVec(i) = e(i);
      theVec(i+order) = s(i);
    }

    return sectInfo.setVector(theVec);
  }

    
  case 105:
    return sectInfo.setVector(this->getSigmaY());
    
  case 106:
    return sectInfo.setVector(this->getTau());
    
  case 107:
    return sectInfo.setVector(this->getAlpha());
    
  case 108:
    return sectInfo.setVector(this->getIter());
    
  case 109:
    return sectInfo.setVector(this->getEX());
    
  case 110:
    return sectInfo.setVector(this->getEY());
    
  case 111:
    return sectInfo.setVector(this->getE1());
    
  case 112:
    return sectInfo.setVector(this->getE2());
    
  case 113:
    return sectInfo.setVector(this->getSX());
    
  case 114:
    return sectInfo.setVector(this->getSY());
    
  case 115:
    return sectInfo.setVector(this->getS1());
    
  case 116:
    return sectInfo.setVector(this->getS2());
    
  default:
    return -1;
  }
}







// AddingSensitivity:BEGIN ////////////////////////////////////
const Vector &
FiberSection2dInt::getSectionDeformationSensitivity(int gradNumber)
{
	static Vector dummy(2);
	return dummy;
}


const Vector &
FiberSection2dInt::getStressResultantSensitivity(int gradNumber, bool conditional)
{
	static Vector dummy(2);	
	return dummy;
}


const Matrix &
FiberSection2dInt::getSectionTangentSensitivity(int gradNumber)
{
	static Matrix something(2,2);
	something.Zero();
	return something;
}

int
FiberSection2dInt::commitSensitivity(const Vector& defSens, int gradNumber, int numGrads)
{
  return 0;
}

// AddingSensitivity:END ///////////////////////////////////
