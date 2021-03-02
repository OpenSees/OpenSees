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

// $Revision: 1.2 $
// $Date: 2013-07-31 00:00:00 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/KikuchiBearing/KikuchiBearing.h,v $

// Written: Ken Ishii
// Created: Jan 2013
// Modified: Feb 17, 2015
//
// KikuchiBearing model
//
// Description: This file contains the implementation of the KikuchiBearing class.
//              This file contains the function to parse the TCL input
// element KikuchiBearing eleTag? iNode? jNode?
//                        -shape shape? -size size? totalRubber? <-totalHeight totalHeight?>
//                        -nMSS nMSS? -matMSS matMSSTag? <-limDisp limDisp?>
//                        -nMNS nMNS? -matMNS matMNSTag? <-lambda lambda?>
//                        <-orient <x1? x2? x3?> yp1? yp2? yp3?> <-mass m?>
//                        <-noPDInput> <-noTilt> <-adjustPDOutput ci? cj?> <-doBalance limFo? limFi? nIter?>

// This element works in 3-dim 6-dof model.
// The multiple normal springs are arranged parallel to local-x axis.
// The multiple shear springs are distributed in local-yz plane.


#include <ID.h>
#include <Vector.h>

#include <KikuchiBearing.h>

#include <Domain.h>
#include <Node.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Renderer.h>
#include <Information.h>
#include <ElementResponse.h>
#include <UniaxialMaterial.h>

#include <float.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <elementAPI.h>


// Bessel function, Copyright(C) 1996 Takuya OOURA
#ifdef _WIN32
extern "C" double dbesi0(double);
extern "C" double dbesi1(double);
#else
extern double dbesi0(double);
extern double dbesi1(double);
#endif


static bool errDetected(bool ifNoError, const char *msg){
  if (ifNoError){
    opserr << "" << endln;
    opserr << "========================================" << endln;
    opserr << "KikuchiBearing element : input error detected" << endln;
    opserr << "------------------------------" << endln;
  }
  opserr << "  " << msg << endln;
  return false;
};


void* OPS_KikuchiBearing()
{

    //3-dim, 6dof
    int ndm = OPS_GetNDM();
    int ndf = OPS_GetNDF();

    if (ndm != 3 || ndf != 6) {
    	opserr << "ndm=" << ndm << ", ndf=" << ndf << endln;
    	opserr << "WARNING KikuchiBearing command only works when ndm is 3 and ndf is 6" << endln;
    	return 0;
    }

    //arguments (necessary)
    int eleTag;
    int iNode;
    int jNode;

    //arguments (necessary, input with -???)
    int shape;
    double size;
    double totalRubber;
    int nMSS;
    int matMSSTag;
    UniaxialMaterial *matMSS;
    int nMNS;
    int matMNSTag;
    UniaxialMaterial *matMNS;

    //arguments (optional, input with -???)
    double totalHeight = -1.0; //default: Norm(I->J)
    double limDisp = -1.0; //default: INF
    double lambda = -1.0; //default: INF
    Vector oriX(0); //default: local-x Vec(I->J)
    Vector oriYp(3); oriYp(0) = 0.0; oriYp(1) = 1.0; oriYp(2) = 0.0; //default: global-Y
    double mass = 0.0;
    bool ifPDInput = true;
    bool ifTilt = true;
    double adjCi = 0.5;
    double adjCj = 0.5;
    bool ifBalance = false;
    double limFo = -1.0; //default: INF
    double limFi = -1.0; //default: INF
    int nIter = 1;

    // input comfirmation
    int recvShape  = 0;
    int recvSize   = 0;
    int recvHeight = 0;
    int recvNMSS  = 0;
    int recvMatMSS  = 0;
    int recvLimDisp = 0;
    int recvNMNS  = 0;
    int recvMatMNS  = 0;
    int recvLambda = 0;
    int recvOrient = 0;
    int recvMass   = 0;
    int recvIfPD = 0;
    int recvIfTl = 0;
    int recvAdj = 0;
    int recvBal = 0;

    //
    Element *theElement = 0;


    //error flag
    bool ifNoError = true;



    if (OPS_GetNumRemainingInputArgs() < 3)  { //element KikuchiBearing eleTag? iNode? jNode?

	ifNoError = errDetected(ifNoError,"insufficient arguments");

    } else {

    	//argv[2~4]
    	int numdata = 1;
    	if (OPS_GetIntInput(&numdata, &eleTag) < 0) {
    	    ifNoError = errDetected(ifNoError,"invalid eleTag");
    	}

    	if (OPS_GetIntInput(&numdata, &iNode) < 0) {
    	    ifNoError = errDetected(ifNoError,"invalid iNode");
    	}

    	if (OPS_GetIntInput(&numdata, &jNode) < 0) {
    	    ifNoError = errDetected(ifNoError,"invalid jNode");
    	}
    
	//argv[5~]
    	while (OPS_GetNumRemainingInputArgs() > 0) {
      
	    double value;
     	    const char* flag = OPS_GetString();

	    if (strcmp(flag,"-shape")==0 && OPS_GetNumRemainingInputArgs()>0) { // -shape shape?
    		const char* shapeflag = OPS_GetString();
    		if (strcmp(shapeflag,"round") == 0) {
    		    shape = 1; //round
    		} else if (strcmp(shapeflag,"square") == 0) {
    		    shape = 2; //square
    		} else {
    		    ifNoError = errDetected(ifNoError,"invalid shape (\"round\" or \"square\" are available)");
    		}
	
    		recvShape++ ;


    	    } else if (strcmp(flag,"-size")==0 && OPS_GetNumRemainingInputArgs()>1) { // -size size? totalRubber?

    		numdata = 1;
    		if (OPS_GetDoubleInput(&numdata, &size)<0 || size<=0) {
    		    ifNoError = errDetected(ifNoError,"invalid size");
    		}

    		if (OPS_GetDoubleInput(&numdata, &totalRubber)<0 || totalRubber<=0) {
    		    ifNoError = errDetected(ifNoError,"invalid totalRubber");
    		}

    		recvSize++ ;
	
    	    } else if (strcmp(flag,"-totalHeight")==0 && OPS_GetNumRemainingInputArgs()>0) { // -totalHeight totalHeight?

    		numdata = 1;
    		if (OPS_GetDoubleInput(&numdata, &totalHeight)<0 || totalHeight<=0) {
    		    ifNoError = errDetected(ifNoError,"invalid totalHeight");
    		}

    		recvHeight++ ;

    	    } else if (strcmp(flag,"-nMSS")==0 && OPS_GetNumRemainingInputArgs()>0) { // -nMSS nMSS?

    		numdata = 1;
    		if (OPS_GetIntInput(&numdata, &nMSS)<0 || nMSS<=0) {
    		    ifNoError = errDetected(ifNoError,"invalid nMSS");
    		}

    		recvNMSS++ ;


    	    } else if (strcmp(flag,"-matMSS")==0 && OPS_GetNumRemainingInputArgs()>0) { // -matMSS matMSSTag?

    		numdata = 1;
    		if (OPS_GetIntInput(&numdata, &matMSSTag)<0) {
    		    ifNoError = errDetected(ifNoError,"invalid matMSSTag");
    		}

    		matMSS = OPS_getUniaxialMaterial(matMSSTag);
    		if (matMSS == 0)  {
    		    ifNoError = errDetected(ifNoError,"material for MSS model not found");
    		}

    		recvMatMSS++ ;


    	    } else if (strcmp(flag,"-limDisp")==0 && OPS_GetNumRemainingInputArgs()>0) { // <-limDisp limDisp?>

    		numdata = 1;
    		if (OPS_GetDoubleInput(&numdata, &limDisp)<0 || limDisp<0) {
    		    ifNoError = errDetected(ifNoError,"invalid limDisp");
    		}

    		recvLimDisp++ ;


    	    } else if (strcmp(flag,"-nMNS")==0 && OPS_GetNumRemainingInputArgs()>0) { // -nMNS nMNS?

    		numdata = 1;
    		if (OPS_GetIntInput(&numdata, &nMNS)<0 || nMNS<=0) {
    		    ifNoError = errDetected(ifNoError,"invalid nMNS");
    		}

    		recvNMNS++ ;


    	    } else if (strcmp(flag,"-matMNS")==0 && OPS_GetNumRemainingInputArgs()>0) { // -matMNS matMNSTag?
    		numdata = 1;
    		if (OPS_GetIntInput(&numdata, &matMNSTag)<0) {
    		    ifNoError = errDetected(ifNoError,"invalid matMNSTag");
    		}

    		matMNS = OPS_getUniaxialMaterial(matMNSTag);
    		if (matMNS == 0)  {
    		    ifNoError = errDetected(ifNoError,"material for MNS model not found");
    		}

    		recvMatMNS++ ;


    	    } else if (strcmp(flag,"-lambda")==0 && OPS_GetNumRemainingInputArgs()>0) { // <-lambda lambda?>
    		numdata = 1;
    		if (OPS_GetDoubleInput(&numdata, &lambda)<0 || lambda<0) {
    		    ifNoError = errDetected(ifNoError,"invalid lambda");
    		}

    		recvLambda++ ;


    	    } else if (strcmp(flag,"-orient")==0 && OPS_GetNumRemainingInputArgs()>=6) { // <-orient x1? x2? x3? yp1? yp2? yp3?>

    		oriX.resize(3);

    		for (int j=1; j<=3; j++) {
    		    numdata = 1;
    		    if (OPS_GetDoubleInput(&numdata, &value)<0) {
    			ifNoError = errDetected(ifNoError,"invalid orient");
    		    } else {
    			oriX(j-1) = value;
    		    }
    		}

    		for (int j=1; j<=3; j++) {
    		    numdata = 1;
    		    if (OPS_GetDoubleInput(&numdata, &value)<0) {
    			ifNoError = errDetected(ifNoError,"invalid orient");
    		    } else {
    			oriYp(j-1) = value;
    		    }
    		}

    		recvOrient++ ;

    	    } else if (strcmp(flag,"-orient")==0 && OPS_GetNumRemainingInputArgs()>=3) { // <-orient yp1? yp2? yp3?>

    		for (int j=1; j<=3; j++) {
    		    numdata = 1;
    		    if (OPS_GetDoubleInput(&numdata, &value)<0) {
    			ifNoError = errDetected(ifNoError,"invalid orient");
    		    } else {
    			oriYp(j-1) = value;
    		    }
    		}

    		recvOrient++ ;

    	    } else if (strcmp(flag,"-mass")==0 && OPS_GetNumRemainingInputArgs()>0) { // <-mass mass?>

    		numdata = 1;
    		if (OPS_GetDoubleInput(&numdata, &mass)<0 || mass<=0) {
    		    ifNoError = errDetected(ifNoError,"invalid mass");
    		}

    		recvMass++ ;

    	    } else if (strcmp(flag,"-noPDInput")==0) { // <-noPDInput>
	
    		ifPDInput = false;

    		recvIfPD++ ;

    	    } else if (strcmp(flag,"-noTilt")==0) { // <-noTilt>
	
    		ifTilt = false;

    		recvIfTl++ ;
	
    	    } else if (strcmp(flag,"-adjustPDOutput")==0 && OPS_GetNumRemainingInputArgs()>1) { // -adjustPDOutput ci? cj?

    		numdata = 1;
    		if (OPS_GetDoubleInput(&numdata, &adjCi)<0 || adjCi<=0) {
    		    ifNoError = errDetected(ifNoError,"invalid ci");
    		}

    		if (OPS_GetDoubleInput(&numdata, &adjCj)<0 || adjCj<=0) {
    		    ifNoError = errDetected(ifNoError,"invalid cj");
    		}

    		recvAdj++ ;

    	    } else if (strcmp(flag,"-doBalance")==0 && OPS_GetNumRemainingInputArgs()>2) { // -doBalance limFo? limFi? nIter?

    		numdata = 1;
    		if (OPS_GetDoubleInput(&numdata, &limFo)<0 || limFo<=0) {
    		    ifNoError = errDetected(ifNoError,"invalid limFo");
    		}

    		if (OPS_GetDoubleInput(&numdata, &limFi)<0 || limFi<=0) {
    		    ifNoError = errDetected(ifNoError,"invalid limFi");
    		}

    		if (OPS_GetIntInput(&numdata, &nIter)<0 || nIter<=0) {
    		    ifNoError = errDetected(ifNoError,"invalid nIter");
    		}

    		ifBalance = true;

    		recvBal++ ;


    	    } else { //invalid option
	
    		ifNoError = errDetected(ifNoError,"invalid optional arguments");
    		break;

    	    }

    	}
    
    } //end input
  

    // input cofirmation
    // necessary arguments
    if (recvShape != 1)  {
    	char buf[100];
    	sprintf(buf,"wrong number of -shape inputs (got %d inputs, but want 1 input)",recvShape);
    	ifNoError = errDetected(ifNoError,buf);
    }

    if (recvSize != 1)  {
    	char buf[100];
    	sprintf(buf,"wrong number of -size inputs (got %d inputs, but want 1 input)",recvSize);
    	ifNoError = errDetected(ifNoError,buf);
    }

    if (recvNMSS != 1)  {
    	char buf[100];
    	sprintf(buf,"wrong number of -NMSS inputs (got %d inputs, but want 1 input)",recvNMSS);
    	ifNoError = errDetected(ifNoError,buf);
    }

    if (recvMatMSS != 1)  {
    	char buf[100];
    	sprintf(buf,"wrong number of -matMSS inputs (got %d inputs, but want 1 input)",recvMatMSS);
    	ifNoError = errDetected(ifNoError,buf);
    }

    if (recvNMNS != 1)  {
    	char buf[100];
    	sprintf(buf,"wrong number of -NMNS inputs (got %d inputs, but want 1 input)",recvNMNS);
    	ifNoError = errDetected(ifNoError,buf);
    }

    if (recvMatMNS != 1)  {
    	char buf[100];
    	sprintf(buf,"wrong number of -matMNS inputs (got %d inputs, but want 1 input)",recvMatMNS);
    	ifNoError = errDetected(ifNoError,buf);
    }


    //optional arguments
    if (recvHeight >= 2)  {
    	char buf[100];
    	sprintf(buf,"wrong number of -totalHeight inputs (got %d inputs, but want 1 input)",recvHeight);
    	ifNoError = errDetected(ifNoError,buf);
    }

    if (recvLimDisp >= 2)  {
    	char buf[100];
    	sprintf(buf,"wrong number of -limDisp inputs (got %d inputs, but want 1 input)",recvLimDisp);
    	ifNoError = errDetected(ifNoError,buf);
    }

    if (recvLambda >= 2)  {
    	char buf[100];
    	sprintf(buf,"wrong number of -lambda inputs (got %d inputs, but want 1 input)",recvLambda);
    	ifNoError = errDetected(ifNoError,buf);
    }

    if (recvOrient >= 2)  {
    	char buf[100];
    	sprintf(buf,"wrong number of -ori inputs (got %d inputs, but want 1 input)",recvOrient);
    	ifNoError = errDetected(ifNoError,buf);
    }

    if (recvMass >= 2)  {
    	char buf[100];
    	sprintf(buf,"wrong number of -mass inputs (got %d inputs, but want 1 input)",recvMass);
    	ifNoError = errDetected(ifNoError,buf);
    }

    if (recvIfPD >= 2)  {
    	char buf[100];
    	sprintf(buf,"wrong number of -noPDInput inputs (got %d inputs, but want 1 input)",recvIfPD);
    	ifNoError = errDetected(ifNoError,buf);
    }

    if (recvIfTl >= 2)  {
    	char buf[100];
    	sprintf(buf,"wrong number of -noTilt inputs (got %d inputs, but want 1 input)",recvIfTl);
    	ifNoError = errDetected(ifNoError,buf);
    }

    if (recvAdj >= 2)  {
    	char buf[100];
    	sprintf(buf,"wrong number of -adjustPDOutput inputs (got %d inputs, but want 1 input)",recvAdj);
    	ifNoError = errDetected(ifNoError,buf);
    }

    if (recvBal >= 2)  {
    	char buf[100];
    	sprintf(buf,"wrong number of -doBalance inputs (got %d inputs, but want 1 input)",recvBal);
    	ifNoError = errDetected(ifNoError,buf);
    }

  
    //if error detected
    if (!ifNoError) {
    	opserr << "------------------------------" << endln;
    	//input:
    	//printCommand(argc, argv);
    	//want:
    	opserr << "Want: element KikuchiBearing eleTag? iNode? jNode?\n";
    	opserr << "                             -shape shape? -size size? totalRubber? <-totalHeight totalHeight?>\n";
    	opserr << "                             -nMSS nMSS? -matMSS matMSSTag? <-lim limDisp?>\n";
    	opserr << "                             -nMNS nMNS? -matMNS matMNSTag? <-lambda lambda?>\n";
    	opserr << "                             <-orient <x1? x2? x3?> yp1? yp2? yp3?> <-mass m?>\n";
    	opserr << "                             <-noPDInput> <-noTilt> <-adjustPDOutput ci? cj?> <-doBalance limFo? limFi? nIter?>\n";
    	opserr << "========================================" << endln;
    	opserr << "" << endln;
    	return 0;
    }

    // now create the KikuchiBearing
    theElement = new KikuchiBearing(eleTag, iNode, jNode, shape, size, totalRubber, totalHeight, nMSS, matMSS, limDisp, nMNS, matMNS, lambda, oriYp, oriX, mass, ifPDInput, ifTilt, adjCi, adjCj, ifBalance, limFo, limFi, nIter);
    return theElement;

}


// initialize the class wide variables
Matrix KikuchiBearing::theMatrix(12,12);
Vector KikuchiBearing::theVector(12);
Vector KikuchiBearing::theLoad(12);

Vector KikuchiBearing::commitDij18(18);
Vector KikuchiBearing::trialDij18(18);

Vector KikuchiBearing::commitFij(12);
Vector KikuchiBearing::trialFij(12);

Matrix KikuchiBearing::Kij18(18,18);
Matrix KikuchiBearing::Kij18_11(12,12);
Matrix KikuchiBearing::Kij18_12(12,6);
Matrix KikuchiBearing::Kij18_21(6,12);
Matrix KikuchiBearing::Kij18_22(6,6);
Matrix KikuchiBearing::invKij18_22(6,6);
Matrix KikuchiBearing::Kij(12,12);

Vector KikuchiBearing::Fij(12);
Vector KikuchiBearing::Fmn(6);


Vector KikuchiBearing::stfCpnt(19);
Vector KikuchiBearing::frcCpnt(12);
Vector KikuchiBearing::dspCpnt(9);



KikuchiBearing::KikuchiBearing(int Tag, int Nd1, int Nd2,
			       int Shape, double Size, double TotalRubber, double TotalHeight,
			       int NMSS, UniaxialMaterial *MatMSS, double LimDisp,
			       int NMNS, UniaxialMaterial *MatMNS, double Lambda,
			       const Vector OriYp, const Vector OriX, double Mass,
			       bool IfPDInput, bool IfTilt,
			       double AdjCi, double AdjCj,
			       bool IfBalance, double LimFo, double LimFi, int NIter)
  : Element(Tag, ELE_TAG_KikuchiBearing),
    connectedExternalNodes(2),
    shape(Shape),size(Size),totalRubber(TotalRubber), totalHeight(TotalHeight),
    nMSS(NMSS),limDisp(LimDisp),
    nMNS(NMNS),lambda(Lambda),
    oriX(OriX), oriYp(OriYp), mass(Mass),
    ifPDInput(IfPDInput), ifTilt(IfTilt),
    adjCi(AdjCi), adjCj(AdjCj),
    ifBalance(IfBalance), limFo(LimFo), limFi(LimFi), nIter(NIter),
    Tgl(12,12), Tlb(6,12),
    basicDisp(6), localDisp(12), basicForce(6),
    localIncrDisp(12),incrDispij(12),incrDispmn(6),localForceij(12)
{
  // ensure the connectedExternalNode ID is of correct size & set values
  if (connectedExternalNodes.Size() != 2)  {
    opserr << "KikuchiBearing::setUp() - element: "
	   << this->getTag() << " failed to create an ID of size 2\n";
  }
  
  connectedExternalNodes(0) = Nd1;
  connectedExternalNodes(1) = Nd2;
  
  // set node pointers to NULL
  for (int i=0; i<2; i++)
    theNodes[i] = 0;
  
  // check material input

  //MNS
  if (MatMNS == 0)  {
    opserr << "KikuchiBearing::KikuchiBearing() - "
	   << "null uniaxial material pointer passed.\n";
    exit(-1);
  }

  theINodeMNSMaterials = new UniaxialMaterial* [nMNS*nMNS]; // nMNS^2 materials, MNS(i-Node)
  theJNodeMNSMaterials = new UniaxialMaterial* [nMNS*nMNS]; // nMNS^2 materials, MNS(j-Node)

  //MSS
  if (MatMSS == 0)  {
    opserr << "KikuchiBearing::KikuchiBearing() - "
	   << "null uniaxial material pointer passed.\n";
    exit(-1);
  }

  theMidMSSMaterials = new UniaxialMaterial* [nMSS]; // nMSS materials, MSS(mid-height)



  // get copies of the uniaxial materials
  //MNS
  for (int i=0; i<nMNS*nMNS; i++)  {
    theINodeMNSMaterials[i] = MatMNS->getCopy();
    if (theINodeMNSMaterials[i] == 0) {
      opserr << "KikuchiBearing::KikuchiBearing() - "
 	     << "failed to copy uniaxial material.\n";
      exit(-1);
    }

    theJNodeMNSMaterials[i] = MatMNS->getCopy();
    if (theJNodeMNSMaterials[i] == 0) {
      opserr << "KikuchiBearing::KikuchiBearing() - "
 	     << "failed to copy uniaxial material.\n";
      exit(-1);
    }
  }
  
  //MSS
  for (int i=0; i<nMSS; i++)  {
    theMidMSSMaterials[i] = MatMSS->getCopy();

    if (theMidMSSMaterials[i] == 0) {
      opserr << "KikuchiBearing::KikuchiBearing() - "
 	     << "failed to copy uniaxial material.\n";
      exit(-1);
    }
  }


  //material to calculate Feq and Seq
  dmyMSSMaterial = MatMSS->getCopy();
  if (dmyMSSMaterial == 0) {
    opserr << "KikuchiBearing::KikuchiBearing() - "
	   << "failed to copy uniaxial material.\n";
    exit(-1);
  }


  //MNS
  //position of centroid, distribution factor
  posLy = new double [nMNS*nMNS]; //local-y position
  posLz = new double [nMNS*nMNS]; //local-z position
  distFct = new double [nMNS*nMNS]; //distribution factor
  if (shape == 1){//round shape

    incA = (M_PI*size*size)/(4.0*nMNS*nMNS); //area of each normal spring

    int p = nMNS%2;
    int k = -1; //index
  
    for(int i=1; i<=((nMNS+p)/2); i++) {
      //circle
      if(i==1 && p==1) {
	k++; //k=0
	double r = 0.0;
	posLy[k] = 0.0;
	posLz[k] = 0.0;
	
	//distribution factor
	if(lambda < 0) { //uniform
	  distFct[k] = 1.0;
	} else if(lambda == 0) { //parabolic
	  distFct[k] = 2.0;
	} else {
	  distFct[k] = (1.0-1.0/dbesi0(lambda));
	}

	continue;
      }

      //circular sector
      double tht = M_PI/(4.0*(2*i-1-p));
      double r1 = ((2.0*i-2-p)/(2.0*nMNS))*size;
      double r2 = ((2.0*i-p)/(2.0*nMNS))*size;
      double r = (2.0/3.0)*(sin(tht)/tht)*((r1*r1+r1*r2+r2*r2)/(r1+r2));
      for(int j=1; j<=(4*(2*i-1-p)); j++){
	k++; //k=(2*i-2-p)^2+j
	double alp = (2*j-1)*tht;
	posLy[k] = r*cos(alp);
	posLz[k] = r*sin(alp);

	//distribution factor
	if(lambda < 0) {
	  distFct[k] = 1.0; //uniform
	} else if(lambda == 0) {
	  distFct[k] = 2.0*(1.0-r/(size/2)*r/(size/2)); //parabolic
	} else {
	  distFct[k] = (1.0-dbesi0(lambda*r/(size/2))/dbesi0(lambda));
	}
      }
    }

  } else { //square shape

    incA = (size*size)/(nMNS*nMNS); //area of each normal spring

    int k = -1; //index
 
    for(int i=0; i<nMNS; i++) {
      for(int j=0; j<nMNS; j++) {

	k++;

	posLy[k] = (2.0*i-nMNS+1)/nMNS * (size/2);
	posLz[k] = (2.0*j-nMNS+1)/nMNS * (size/2);

	double ry = posLy[k]/(size/2);
	double rz = posLz[k]/(size/2);
	
	//distribution factor
	if(lambda < 0) { //uniform
	  distFct[k] = 1.0;
	} else {
	  distFct[k] = 0.0;
	  for (int m=1; m<=19; m+=2) {
	    double alpha = sqrt(lambda*lambda + (m*M_PI/2.0)*(m*M_PI/2.0));
	    distFct[k] += sin(m*M_PI/2.0) * 4/(m*M_PI) * (lambda*lambda)/(alpha*alpha) * (1.0 - cosh(alpha*rz)/cosh(alpha)) * cos(m*M_PI/2.0*ry);
	  }
	}
      }
    }
    
  }

  //ave(distFct) = 1.0
  double aveDistFct = 0.0;
  for(int k=0; k<nMNS*nMNS; k++) {
    aveDistFct += distFct[k];//sumDistFct
  }
  aveDistFct = aveDistFct / (nMNS*nMNS);
  for(int k=0; k<nMNS*nMNS; k++) {
    distFct[k] = distFct[k] / aveDistFct;  
  }

  
  //commitStrain for MNS
  commitStrnIMns = new double [nMNS*nMNS];
  commitStrnJMns = new double [nMNS*nMNS];
  for (int i=0; i<nMNS*nMNS; i++)  {
    commitStrnIMns[i] = 0.0;
    commitStrnJMns[i] = 0.0;
  }

  
  //MSS
  //arrangement of each spring
  cosTht = new double [nMSS];
  sinTht = new double [nMSS];
  
  for (int i=0; i<nMSS; i++)  {
    cosTht[i] = cos(M_PI*i/nMSS);
    sinTht[i] = sin(M_PI*i/nMSS);
  }

  //commitDisp for MSS
  commitDspMss = new double [nMSS];
  for (int i=0; i<nMSS; i++)  {
    commitDspMss[i] = 0.0;
  }


  //stiff springs in mid-height
  double sect_aa = 0.0; //area
  double sect_ii = 0.0; //moment of inertia of area
  double sect_pp = 0.0; //polar moment of inertia of area
  
  if (shape==1) { //round
    sect_aa = (M_PI*size*size)/(4.0);
    sect_ii = (M_PI*size*size*size*size)/(64.0);
    sect_pp = (sect_ii)*(2.0);
  } else { //square
    sect_aa = (size*size);
    sect_ii = (size*size*size*size)/(12.0);
    sect_pp = (sect_ii)*(2.0);
  }

  stfMidX  = ((theINodeMNSMaterials[0]->getInitialTangent()) * sect_aa/(totalRubber)) * 1e4;
  stfMidRY = ((theINodeMNSMaterials[0]->getInitialTangent()) * sect_ii/(totalRubber)) * 1e4;
  stfMidRZ = ((theINodeMNSMaterials[0]->getInitialTangent()) * sect_ii/(totalRubber)) * 1e4;
  stfMidRX = ((theMidMSSMaterials[0]->getInitialTangent()) * sect_pp/(totalRubber)) * 1e4;

  // initialize variables
  this->revertToStart();
  
}


KikuchiBearing::KikuchiBearing()
  : Element(0, ELE_TAG_KikuchiBearing),
    connectedExternalNodes(2),
    shape(0),size(0.0),totalRubber(0.0), totalHeight(0.0),
    nMSS(0),limDisp(0.0),
    nMNS(0),lambda(0.0),
    oriX(0), oriYp(0), mass(0.0),
    ifPDInput(false), ifTilt(false),
    adjCi(0.0), adjCj(0.0),
    ifBalance(false), limFo(0.0), limFi(0.0), nIter(0),
    Tgl(12,12), Tlb(6,12),
    basicDisp(6), localDisp(12), basicForce(6),
    localIncrDisp(12),incrDispij(12),incrDispmn(6),localForceij(12)
{	


  // ensure the connectedExternalNode ID is of correct size & set values
  if (connectedExternalNodes.Size() != 2)  {
    opserr << "KikuchiBearing::KikuchiBearing() - "
	   <<  "failed to create an ID of size 2\n";
    exit(-1);
  }
  
  // set node pointers to NULL
  for (int i=0; i<2; i++)
    theNodes[i] = 0;

  // set material array pointers to NULL
  theINodeMNSMaterials = 0;
  theJNodeMNSMaterials = 0;
  theMidMSSMaterials = 0;
  dmyMSSMaterial = 0;

}


KikuchiBearing::~KikuchiBearing()
{
  // invoke the destructor on any objects created by the object
  // that the object still holds a pointer to


  // MSS
  if (theMidMSSMaterials != 0) {
    for (int i=0; i<nMSS; i++)
      if (theMidMSSMaterials[i] != 0)
	delete theMidMSSMaterials[i];
    delete [] theMidMSSMaterials;
  }

  if (cosTht != 0)
    delete [] cosTht;

  if (sinTht != 0)
    delete [] sinTht;

  if (dmyMSSMaterial != 0)
    delete dmyMSSMaterial;

  if (commitDspMss != 0)
    delete [] commitDspMss;

  // MNS
  if (theINodeMNSMaterials != 0) {
    for (int i=0; i<(nMNS*nMNS); i++)
      if (theINodeMNSMaterials[i] != 0)
	delete theINodeMNSMaterials[i];
    delete [] theINodeMNSMaterials;
  }

  if (theJNodeMNSMaterials != 0) {
    for (int i=0; i<(nMNS*nMNS); i++)
      if (theJNodeMNSMaterials[i] != 0)
	delete theJNodeMNSMaterials[i];
    delete [] theJNodeMNSMaterials;
  }

  if (posLy != 0)
    delete [] posLy;

  if (posLz != 0)
    delete [] posLz;

  if (distFct != 0)
    delete [] distFct;

  if (commitStrnIMns != 0)
    delete [] commitStrnIMns;

  if (commitStrnJMns != 0)
    delete [] commitStrnJMns;


}


int KikuchiBearing::getNumExternalNodes() const
{
  return 2;
}


const ID& KikuchiBearing::getExternalNodes() 
{
  return connectedExternalNodes;
}


Node** KikuchiBearing::getNodePtrs() 
{
  return theNodes;
}


int KikuchiBearing::getNumDOF() 
{
  return 12;
}


void KikuchiBearing::setDomain(Domain *theDomain)
{
  // check Domain is not null - invoked when object removed from a domain
  if (!theDomain)  {
    theNodes[0] = 0;
    theNodes[1] = 0;
    
    return;
  }
  
  // first set the node pointers
  theNodes[0] = theDomain->getNode(connectedExternalNodes(0));
  theNodes[1] = theDomain->getNode(connectedExternalNodes(1));	
  
  // if can't find both - send a warning message
  if (!theNodes[0] || !theNodes[1])  {
    if (!theNodes[0])  {
      opserr << "WARNING KikuchiBearing::setDomain() - Nd1: " 
	     << connectedExternalNodes(0) << " does not exist in the model for ";
    } else  {
      opserr << "WARNING KikuchiBearing::setDomain() - Nd2: " 
	     << connectedExternalNodes(1) << " does not exist in the model for ";
    }
    opserr << "KikuchiBearing ele: " << this->getTag() << endln;
    
    return;
  }
  
  // now determine the number of dof and the dimension    
  int dofNd1 = theNodes[0]->getNumberDOF();
  int dofNd2 = theNodes[1]->getNumberDOF();	
  
  // if differing dof at the ends - print a warning message
  if (dofNd1 != 6)  {
    opserr << "KikuchiBearing::setDomain() - node 1: "
	   << connectedExternalNodes(0) << " has incorrect number of DOF (not 6)\n";
    return;
  }
  if (dofNd2 != 6)  {
    opserr << "KikuchiBearing::setDomain() - node 2: "
	   << connectedExternalNodes(1) << " has incorrect number of DOF (not 6)\n";
    return;
  }
  
  // call the base class method
  this->DomainComponent::setDomain(theDomain);
  
  // set up the transformation matrix for orientation
  this->setUp();
}   	 


int KikuchiBearing::commitState()
{
  int errCode = 0;

  //get rid of internal unbalanced force
  int iite = 0;

  while(ifBalance){
    //number of iterations
    if(iite>=nIter){
      opserr << "KikuchiBearing::KikuchiBearing() - "
	     << "inner iteration failed (lmtI) \n";
      break;
    }
    
    //refer to finite displacement in the element (trial)
    subRefFntDisp(false);
    
    //calculate stiffness components
    subCalcStfCpnt();
    
    //calculate force components
    subCalcFrcCpnt();

    //make K18 matrix (full)
    subMakeKij18();

    //make submatrices
    subSubmatKij18();
    
    //calculate Fmn (and Fij)
    subMakeFijFmn();

    //incrDispij is 0
    //calculate incrDispmn
    Kij18_22.Invert(invKij18_22);
    incrDispmn = (-1.0) * (invKij18_22 *  Fmn );

    //move internal nodes
    for (int j=0; j<=5; j++)  {
      trialDij18(j+12) += incrDispmn(j);
    }
    
    //setTrialStrain for materials
    incrDispij.Zero();
    subSetMaterialStrains(false);

    //calculate stiffness components and force components again
    subCalcStfCpnt();
    subCalcFrcCpnt();

    //make K18 matrix (full)
    subMakeKij18();
    
    //make submatrices
    subSubmatKij18();
    
    //calculate Fij and Fmn
    subMakeFijFmn();
    
    //calculate equivalent internal force vector
    Kij18_22.Invert(invKij18_22); // inv_Kij18_22 = inv(Kij18_22)
    theVector = (-1.0) * (Kij18_12 * invKij18_22 * Fmn) + Fij;

 
    //external unbalanced force
    for (int i=0; i<=11; i++)  {
      if(fabs(localForceij(i)-theVector(i))>limFo){
	opserr << "KikuchiBearing::KikuchiBearing() - "
	       << "inner iteration failed (lmtO) \n";
      }
    }

    //internal unbalanced force
    bool ifContinue = false;
    for (int i=0; i<=5; i++)  {
      if(fabs(Fmn(i))>limFi){
	ifContinue = true;
      }
    }

    if(ifContinue){
      iite++;
      opserr << "inner iteration\n";
      continue;
    } else {
      opserr << "inner iteration done\n";
      break;
    }

  }
  


  // commit material models

  // MSS
  for (int i=0; i<nMSS; i++) {
    errCode += theMidMSSMaterials[i]->commitState();

    commitDspMss[i] = theMidMSSMaterials[i]->getStrain();
  }


  // MNS
  for (int i=0; i<(nMNS*nMNS); i++) {
    errCode += theINodeMNSMaterials[i]->commitState();
    errCode += theJNodeMNSMaterials[i]->commitState();

    commitStrnIMns[i] = theINodeMNSMaterials[i]->getStrain();
    commitStrnJMns[i] = theJNodeMNSMaterials[i]->getStrain();
  }


  // stiff springs in mid-height
  commitDspMidX  = trialDspMidX;
  commitDspMidRY = trialDspMidRY;
  commitDspMidRZ = trialDspMidRZ;
  commitDspMidRX = trialDspMidRX;
  commitFrcMidX  = trialFrcMidX;
  commitFrcMidRY = trialFrcMidRY;
  commitFrcMidRZ = trialFrcMidRZ;
  commitFrcMidRX = trialFrcMidRX;


  //Dij, Fij
  commitDij18 = trialDij18;
  commitFij = trialFij;
  
  
  //equivalent coefficient for MSS
  subCalcMSSFeqSeq();

  return errCode;

}


int KikuchiBearing::revertToLastCommit()
{
  int errCode = 0;
  
  // revert material models
  for (int i=0; i<(nMNS*nMNS); i++)
    errCode += theINodeMNSMaterials[i]->revertToLastCommit();

  for (int i=0; i<(nMNS*nMNS); i++)
    errCode += theJNodeMNSMaterials[i]->revertToLastCommit();

  for (int i=0; i<nMSS; i++)
    errCode += theMidMSSMaterials[i]->revertToLastCommit();


  return errCode;
}


int KikuchiBearing::revertToStart()
{   

  int errCode=0;

  // trial variables
  basicDisp.Zero();
  basicForce.Zero();

  // Dij, Fij
  commitDij18.Zero();
  trialDij18.Zero();

  commitFij.Zero();
  trialFij.Zero();

  dspCpnt.Zero();
  
  // revert material models
  for (int i=0; i<(nMNS*nMNS); i++) {
    errCode += theINodeMNSMaterials[i]->revertToStart();
    errCode += theJNodeMNSMaterials[i]->revertToStart();
    
    commitStrnIMns[i] = 0.0;
    commitStrnJMns[i] = 0.0;
  }

  for (int i=0; i<nMSS; i++) {
    errCode += theMidMSSMaterials[i]->revertToStart();

    commitDspMss[i] = 0.0;
  }

  dmyMSSMaterial->revertToStart();
  

  commitDspMidX  = 0.0;
  commitDspMidRY = 0.0;
  commitDspMidRZ = 0.0;
  commitDspMidRX = 0.0;
  trialDspMidX  = 0.0;
  trialDspMidRY = 0.0;
  trialDspMidRZ = 0.0;
  trialDspMidRX = 0.0;


  //equivalent coefficient for MSS
  subCalcMSSFeqSeq();
  
  //initialStiff
  subCalcStfCpntInit();
  subMakeKij18();
  subReductKij();

  return errCode;
}


int KikuchiBearing::update()
{

  // get global trial displacements and velocities
  const Vector &dsp1 = theNodes[0]->getTrialDisp();
  const Vector &dsp2 = theNodes[1]->getTrialDisp();
  
  static Vector globalDisp(12);
  for (int i=0; i<6; i++)  {
    globalDisp(i)   = dsp1(i);
    globalDisp(i+6) = dsp2(i);
  }

  
  // transform response from the global to the local system
  localDisp    = Tgl*globalDisp;
  
  // transform response from the local to the basic system
  basicDisp    = Tlb*localDisp;


  //-------------------------------------

  //---0 incremental displacement of external nodes
  // get globalIncrDisp (global-X,Y,Z,RX,RY,RZ)
  const Vector &globalIncrDispINode = theNodes[0]->getIncrDisp();
  const Vector &globalIncrDispJNode = theNodes[1]->getIncrDisp();

  static Vector globalIncrDisp(12);
  for (int i=0; i<6; i++)  {
    globalIncrDisp(i)   = globalIncrDispINode(i);
    globalIncrDisp(i+6) = globalIncrDispJNode(i);
  }
  
  // transform (global-X,Y,Z,RX,RY,RZ) -> (local-x,y,z,rx,ry,rz)
  localIncrDisp = Tgl*globalIncrDisp;

  //---1

  //refer to finite displacement in the element (commit)
  subRefFntDisp();

  //calculate stiffness components
  subCalcStfCpnt();

  //calculate force components
  subCalcFrcCpnt();

  //make K18 matrix (full)
  subMakeKij18();

  //make submatrices
  subSubmatKij18();

  //calculate Fij and Fmn
  subMakeFijFmn();

  //incremental displacement of external nodes
  incrDispij = localIncrDisp;

  //incremental displacement of internal nodes
  Kij18_22.Invert(invKij18_22); // inv_Kij18_22 = inv(Kij18_22)
  incrDispmn = (-1.0) * (invKij18_22 * ( Kij18_21 * incrDispij + Fmn ));


  //---2

  //setTrialStrain for materials
  subSetMaterialStrains();

  //calculate stiffness components and force components again
  subCalcStfCpnt();
  subCalcFrcCpnt();

  //make K18 matrix
  subMakeKij18();

  //make submatrices
  subSubmatKij18();

  //calculate Fij and Fmn
  subMakeFijFmn();

  //calculate equivalent internal force vector
  Kij18_22.Invert(invKij18_22); // inv_Kij18_22 = inv(Kij18_22)
  localForceij = (-1.0) * (Kij18_12 * invKij18_22 * Fmn) + Fij;


  if (ifPDInput) {
    //P-Delta moment adjustment for reaction force
    double hdsb = ( commitDij18(7)+incrDispij(7) ) - ( commitDij18(1)+incrDispij(1) );// local-y relative deformation between i&j
    double hdsc = ( commitDij18(8)+incrDispij(8) ) - ( commitDij18(2)+incrDispij(2) );// local-z relative deformation between i&j
    
    double fnrm  = -localForceij(6);// axial force (j-node)
    
    double pdfbi =  fnrm*hdsb*adjCi;
    double pdfbj =  fnrm*hdsb*adjCj;
    double pdfci =  fnrm*hdsc*adjCi;
    double pdfcj =  fnrm*hdsc*adjCj;
    
    localForceij(4)  +=  -pdfbi  ;// local-ry (i-node)
    localForceij(5)  +=   pdfci  ;// local-rz (i)
    localForceij(10) +=  -pdfbj  ;// local-ry (j)
    localForceij(11) +=   pdfcj  ;// local-rz (j)
  }
  

  //reduct K18 matrix
  subReductKij();


  //---3
  for (int j=0; j<12; j++)  {
    trialDij18(j) = commitDij18(j) + incrDispij(j);
  }
  for (int j=0; j<6; j++)  {
    trialDij18(j+12) = commitDij18(j+12) + incrDispmn(j);
  }

  trialFij = localForceij;

  //---4
  basicForce = 0.5*Tlb*localForceij;

  return 0;
}


const Matrix& KikuchiBearing::getTangentStiff()
{
  // zero the matrix
  theMatrix.Zero();

  // transform from basic to local system
  static Matrix localStiff(12,12);
  localStiff = Kij;

  // transform from local to global system
  theMatrix.addMatrixTripleProduct(0.0, Tgl, localStiff, 1.0);

  return theMatrix;

}


const Matrix& KikuchiBearing::getInitialStiff()
{

  // zero the matrix
  theMatrix.Zero();
  
  // transform from basic to local system
  static Matrix localStiff(12,12);
  localStiff = Kij;
  
  // transform from local to global system
  theMatrix.addMatrixTripleProduct(0.0, Tgl, localStiff, 1.0);

  return theMatrix;
}


const Matrix& KikuchiBearing::getMass()
{
  // zero the matrix
  theMatrix.Zero();
  
  // check for quick return
  if (mass == 0.0)  {
    return theMatrix;
  }    
  
  double m = 0.5*mass;
  for (int i = 0; i < 3; i++)  {
    theMatrix(i,i)     = m;
    theMatrix(i+6,i+6) = m;
  }
  
  return theMatrix; 
}


void KikuchiBearing::zeroLoad()
{
  theLoad.Zero();
}


int KikuchiBearing::addLoad(ElementalLoad *theLoad, double loadFactor)
{  
  opserr <<"KikuchiBearing::addLoad() - "
	 << "load type unknown for element: "
	 << this->getTag() << endln;
  
  return -1;
}


int KikuchiBearing::addInertiaLoadToUnbalance(const Vector &accel)
{
  // check for quick return
  if (mass == 0.0)  {
    return 0;
  }
  
  // get R * accel from the nodes
  const Vector &Raccel1 = theNodes[0]->getRV(accel);
  const Vector &Raccel2 = theNodes[1]->getRV(accel);
  
  if (6 != Raccel1.Size() || 6 != Raccel2.Size())  {
    opserr << "KikuchiBearing::addInertiaLoadToUnbalance() - "
	   << "matrix and vector sizes are incompatible\n";
    return -1;
  }
  
  // want to add ( - fact * M R * accel ) to unbalance
  // take advantage of lumped mass matrix
  double m = 0.5*mass;
  for (int i = 0; i < 3; i++)  {
    theLoad(i)   -= m * Raccel1(i);
    theLoad(i+6) -= m * Raccel2(i);
  }
  
  return 0;
}


const Vector& KikuchiBearing::getResistingForce()
{
  // zero the residual
  theVector.Zero();
  
  // determine resisting forces in local system
  static Vector localForce(12);
  localForce = localForceij;
  
  // determine resisting forces in global system
  theVector = Tgl^localForce;

  // subtract external load
  theVector.addVector(1.0, theLoad, -1.0);

  
  return theVector;
}


const Vector& KikuchiBearing::getResistingForceIncInertia()
{	
  theVector = this->getResistingForce();
  
  // add the damping forces if rayleigh damping
  if (alphaM != 0.0 || betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0)
    theVector += this->getRayleighDampingForces();
  
  // now include the mass portion
  if (mass != 0.0)  {
    const Vector &accel1 = theNodes[0]->getTrialAccel();
    const Vector &accel2 = theNodes[1]->getTrialAccel();    
    
    double m = 0.5*mass;
    for (int i = 0; i < 3; i++)  {
      theVector(i)   += m * accel1(i);
      theVector(i+6) += m * accel2(i);
    }
  }
  
  return theVector;
}


int KikuchiBearing::sendSelf(int commitTag, Channel &sChannel)
{
  return -1;
}


int KikuchiBearing::recvSelf(int commitTag, Channel &rChannel,
				  FEM_ObjectBroker &theBroker)
{
  return -1;
}


int KikuchiBearing::displaySelf(Renderer &theViewer,
				int displayMode, float fact,
				const char **modes, int numMode)
{
    static Vector v1(3);
    static Vector v2(3);

    theNodes[0]->getDisplayCrds(v1, fact, displayMode);
    theNodes[1]->getDisplayCrds(v2, fact, displayMode);

    return theViewer.drawLine(v1, v2, 1.0, 1.0, this->getTag());
}


void KikuchiBearing::Print(OPS_Stream &s, int flag)
{
  
  if (flag == 0)  {
    // print everything
    s << "Element: " << this->getTag(); 
    s << "  type: KikuchiBearing  iNode: " << connectedExternalNodes(0);
    s << "  jNode: " << connectedExternalNodes(1) << endln;
    s << "  Material : " << theINodeMNSMaterials[0]->getTag() << endln;
    s << "  mass: " << mass << endln;
    // determine resisting forces in global system
    s << "  resisting force: " << this->getResistingForce() << endln;
  } else if (flag == 1)  {
    // does nothing
  }
}


Response* KikuchiBearing::setResponse(const char **argv, int argc,
					   OPS_Stream &output)
{
  Response *theResponse = 0;
  
  output.tag("ElementOutput");
  output.attr("eleType","KikuchiBearing");
  output.attr("eleTag",this->getTag());
  output.attr("node1",connectedExternalNodes[0]);
  output.attr("node2",connectedExternalNodes[1]);
  
  // global forces
  if (strcmp(argv[0],"force") == 0 || strcmp(argv[0],"forces") == 0 ||
      strcmp(argv[0],"globalForce") == 0 || strcmp(argv[0],"globalForces") == 0)
    {
      output.tag("ResponseType","Px_1");
      output.tag("ResponseType","Py_1");
      output.tag("ResponseType","Pz_1");
      output.tag("ResponseType","Mx_1");
      output.tag("ResponseType","My_1");
      output.tag("ResponseType","Mz_1");
      output.tag("ResponseType","Px_2");
      output.tag("ResponseType","Py_2");
      output.tag("ResponseType","Pz_2");
      output.tag("ResponseType","Mx_2");
      output.tag("ResponseType","My_2");
      output.tag("ResponseType","Mz_2");
      
      theResponse = new ElementResponse(this, 1, theVector);
    }
  // local forces
  else if (strcmp(argv[0],"localForce") == 0 || strcmp(argv[0],"localForces") == 0)
    {
      output.tag("ResponseType","N_ 1");
      output.tag("ResponseType","Vy_1");
      output.tag("ResponseType","Vz_1");
      output.tag("ResponseType","T_1");
      output.tag("ResponseType","My_1");
      output.tag("ResponseType","Tz_1");
      output.tag("ResponseType","N_2");
      output.tag("ResponseType","Py_2");
      output.tag("ResponseType","Pz_2");
      output.tag("ResponseType","T_2");
      output.tag("ResponseType","My_2");
      output.tag("ResponseType","Mz_2");
      
      theResponse = new ElementResponse(this, 2, theVector);
    }
  // basic forces
  else if (strcmp(argv[0],"basicForce") == 0 || strcmp(argv[0],"basicForces") == 0)
    {
      output.tag("ResponseType","qb1");
      output.tag("ResponseType","qb2");
      output.tag("ResponseType","qb3");
      output.tag("ResponseType","qb4");
      output.tag("ResponseType","qb5");
      output.tag("ResponseType","qb6");
      
      theResponse = new ElementResponse(this, 3, Vector(6));
    }
  // local displacements
  else if (strcmp(argv[0],"localDisplacement") == 0 ||
	   strcmp(argv[0],"localDisplacements") == 0)
    {
      output.tag("ResponseType","ux_1");
      output.tag("ResponseType","uy_1");
      output.tag("ResponseType","uz_1");
      output.tag("ResponseType","rx_1");
      output.tag("ResponseType","ry_1");
      output.tag("ResponseType","rz_1");
      output.tag("ResponseType","ux_2");
      output.tag("ResponseType","uy_2");
      output.tag("ResponseType","uz_2");
      output.tag("ResponseType","rx_2");
      output.tag("ResponseType","ry_2");
      output.tag("ResponseType","rz_2");
      
      theResponse = new ElementResponse(this, 4, theVector);
    }
  // basic displacements
  else if (strcmp(argv[0],"deformation") == 0 || strcmp(argv[0],"deformations") == 0 || 
	   strcmp(argv[0],"basicDeformation") == 0 || strcmp(argv[0],"basicDeformations") == 0 ||
	   strcmp(argv[0],"basicDisplacement") == 0 || strcmp(argv[0],"basicDisplacements") == 0)
    {
      output.tag("ResponseType","ub1");
      output.tag("ResponseType","ub2");
      output.tag("ResponseType","ub3");
      output.tag("ResponseType","ub4");
      output.tag("ResponseType","ub5");
      output.tag("ResponseType","ub6");
      
      theResponse = new ElementResponse(this, 5, Vector(6));
    }
  
  output.endTag(); // ElementOutput
  
  return theResponse;
}


int KikuchiBearing::getResponse(int responseID, Information &eleInfo)
{
  switch (responseID)  {
  case 1:  // global forces
    return eleInfo.setVector(this->getResistingForce());
    
  case 2:  // local forces
    return eleInfo.setVector(localForceij);
    
  case 3:  // basic forces
    return eleInfo.setVector(basicForce);
    
  case 4:  // local displacements
    return eleInfo.setVector(localDisp);
    
  case 5:  // basic displacements
    return eleInfo.setVector(basicDisp);
    
  default:
    return -1;
  }
}


// set up the transformation matrix for orientation
void KikuchiBearing::setUp()
{ 

  const Vector &end1Crd = theNodes[0]->getCrds();
  const Vector &end2Crd = theNodes[1]->getCrds();	
  Vector oriXp = end2Crd - end1Crd;
  if (totalHeight < 0)  totalHeight = oriXp.Norm(); //total height of bearing, default: Norm(I->J)

  if (totalHeight > DBL_EPSILON)  {
    if (oriX.Size() == 0)  {
      oriX.resize(3);
      oriX = oriXp;
    } else  {
      opserr << "WARNING KikuchiBearing::setUp() - " 
	     << "element: " << this->getTag() << endln
	     << "ignoring nodes and using specified "
	     << "local x vector to determine orientation\n";
    }
  }
  // check that vectors for orientation are of correct size
  if (oriX.Size() != 3 || oriYp.Size() != 3)  {
    opserr << "KikuchiBearing::setUp() - "
	   << "element: " << this->getTag() << endln
	   << "incorrect dimension of orientation vectors\n";
    exit(-1);
  }

    
  // establish orientation of element for the transformation matrix
  // z = x cross yp
  Vector oriZ(3);
  oriZ(0) = oriX(1)*oriYp(2) - oriX(2)*oriYp(1);
  oriZ(1) = oriX(2)*oriYp(0) - oriX(0)*oriYp(2);
  oriZ(2) = oriX(0)*oriYp(1) - oriX(1)*oriYp(0);
  
  // y = z cross x
  Vector oriY(3);
  oriY(0) = oriZ(1)*oriX(2) - oriZ(2)*oriX(1);
  oriY(1) = oriZ(2)*oriX(0) - oriZ(0)*oriX(2);
  oriY(2) = oriZ(0)*oriX(1) - oriZ(1)*oriX(0);
    
  // compute length(norm) of vectors
  double xn = oriX.Norm();
  double yn = oriY.Norm();
  double zn = oriZ.Norm();
  
  // check valid x and y vectors, i.e. not parallel and of zero length
  if (xn == 0 || yn == 0 || zn == 0)  {
    opserr << "KikuchiBearing::setUp() - "
	   << "element: " << this->getTag() << endln
	   << "invalid orientation vectors\n";
    exit(-1);
  }
  
  // create transformation matrix from global to local system
  Tgl.Zero();
  Tgl(0,0) = Tgl(3,3) = Tgl(6,6) = Tgl(9,9)   = oriX(0)/xn;
  Tgl(0,1) = Tgl(3,4) = Tgl(6,7) = Tgl(9,10)  = oriX(1)/xn;
  Tgl(0,2) = Tgl(3,5) = Tgl(6,8) = Tgl(9,11)  = oriX(2)/xn;
  Tgl(1,0) = Tgl(4,3) = Tgl(7,6) = Tgl(10,9)  = oriY(0)/yn;
  Tgl(1,1) = Tgl(4,4) = Tgl(7,7) = Tgl(10,10) = oriY(1)/yn;
  Tgl(1,2) = Tgl(4,5) = Tgl(7,8) = Tgl(10,11) = oriY(2)/yn;
  Tgl(2,0) = Tgl(5,3) = Tgl(8,6) = Tgl(11,9)  = oriZ(0)/zn;
  Tgl(2,1) = Tgl(5,4) = Tgl(8,7) = Tgl(11,10) = oriZ(1)/zn;
  Tgl(2,2) = Tgl(5,5) = Tgl(8,8) = Tgl(11,11) = oriZ(2)/zn;
  
  // create transformation matrix from local to basic system (linear)
  Tlb.Zero();
  Tlb(0,0) = Tlb(1,1) = Tlb(2,2) = Tlb(3,3) = Tlb(4,4) = Tlb(5,5) = -1.0;
  Tlb(0,6) = Tlb(1,7) = Tlb(2,8) = Tlb(3,9) = Tlb(4,10) = Tlb(5,11) = 1.0;

}

// subroutine

//refer to finite displacement in the element
void KikuchiBearing::subRefFntDisp(bool ifCommit)
{

  //dspCpnt
  double tym  = 0.0 ;//(0)
  double tzm  = 0.0 ;//(1)
  double tyn  = 0.0 ;//(2)
  double tzn  = 0.0 ;//(3)
  double dsy  = 0.0 ;//(4)
  double dsz  = 0.0 ;//(5)
  double hf   = 0.0 ;//(6)
  double dsy2 = 0.0 ;//(7)
  double dsz2 = 0.0 ;//(8)

  //commitDij or trialDij
  double uiy = 0.0;
  double uiz = 0.0;
  double ujy = 0.0;
  double ujz = 0.0;
    


  if (ifCommit) {
    //displacement of old step
    uiy = commitDij18(1);
    uiz = commitDij18(2);
    ujy = commitDij18(7);
    ujz = commitDij18(8);
  } else {
    //displacement of new step
    uiy = trialDij18(1);
    uiz = trialDij18(2);
    ujy = trialDij18(7);
    ujz = trialDij18(8);
  }
  

  //consider tilt of rigid link (or not)
  if (ifTilt) {
    if (ifCommit) {
      tym = commitDij18(13);
      tzm = commitDij18(14);
      tyn = commitDij18(16);
      tzn = commitDij18(17);
    } else {
      tym = trialDij18(13);
      tzm = trialDij18(14);
      tyn = trialDij18(16);
      tzn = trialDij18(17);
    }
  } else {
    tym = 0.0;
    tzm = 0.0;
    tyn = 0.0;
    tzn = 0.0;
  }
  

  //consider P-Delta moment (or not)
  if (ifPDInput) {
    dsy = (ujy-uiy) - (tzm+tzn)*totalHeight/2.0;
    dsz = (ujz-uiz) + (tym+tyn)*totalHeight/2.0;
  } else {
    dsy = 0.0;
    dsz = 0.0;
  }
  

  hf = totalHeight/2.0;
  dsy2 = dsy/2.0;
  dsz2 = dsz/2.0;


  //dspCpnt output
  dspCpnt(0) = tym;
  dspCpnt(1) = tzm;
  dspCpnt(2) = tyn;
  dspCpnt(3) = tzn;
  dspCpnt(4) = dsy;
  dspCpnt(5) = dsz;

  dspCpnt(6) = hf;
  dspCpnt(7) = dsy2;
  dspCpnt(8) = dsz2;


  return;
}

//setTrialStrain for materials
void KikuchiBearing::subSetMaterialStrains(bool ifCommit)
{
  //dspCpnt
  double tym  = dspCpnt(0);
  double tzm  = dspCpnt(1);
  double tyn  = dspCpnt(2);
  double tzn  = dspCpnt(3);
  double dsy  = dspCpnt(4);
  double dsz  = dspCpnt(5);
  double hf   = dspCpnt(6);
  double dsy2 = dspCpnt(7);
  double dsz2 = dspCpnt(8);
  
  //incremental displacement of m'-node
  double dxmd = incrDispmn(0) + tzm*incrDispij(1) - tym*incrDispij(2) + dsz2*incrDispmn(1) - dsy2*incrDispmn(2) ;
  double dymd = incrDispij(1) +  hf*incrDispmn(2);
  double dzmd = incrDispij(2) -  hf*incrDispmn(1);
  double txmd = incrDispij(3);
  double tymd = incrDispmn(1);
  double tzmd = incrDispmn(2);

  //incremental displacement of n'-node
  double dxnd = incrDispmn(3) + tzn*incrDispij(7) - tyn*incrDispij(8) - dsz2*incrDispmn(4) + dsy2*incrDispmn(5) ;
  double dynd = incrDispij(7) -  hf*incrDispmn(5);
  double dznd = incrDispij(8) +  hf*incrDispmn(4);
  double txnd = incrDispij(9);
  double tynd = incrDispmn(4);
  double tznd = incrDispmn(5);
  
  //relative incremental deformation between  m'&n'
  double dsxx = dxnd - dxmd;
  double dsyy = dynd - dymd;
  double dszz = dznd - dzmd;
  double drxx = txnd - txmd;
  double dryy = tynd - tymd;
  double drzz = tznd - tzmd;


  //MSS (force-deformation)
  for (int i=0; i<nMSS; i++)  {
    double tmpCommitDisp  = (ifCommit) ? commitDspMss[i] : theMidMSSMaterials[i]->getStrain() ;
    double tmpIncrDisp = dsyy * cosTht[i] + dszz * sinTht[i];
    double tmpDisp = tmpCommitDisp + tmpIncrDisp;

    theMidMSSMaterials[i]->setTrialStrain(tmpDisp);
  }
  

  //i-MNS (stress-strain, length=totalRubber/2)
  for (int i=0; i<(nMNS*nMNS); i++)  {
    double tmpCommitStrain = (ifCommit) ? commitStrnIMns[i] : theINodeMNSMaterials[i]->getStrain();
    double tmpIncrDisp = incrDispmn(0) + posLz[i]*incrDispmn(1) - posLy[i]*incrDispmn(2) - incrDispij(0) - posLz[i]*incrDispij(4) + posLy[i]*incrDispij(5) ;
    double tmpStrain = tmpCommitStrain +  tmpIncrDisp/(totalRubber/2.0);
    theINodeMNSMaterials[i]->setTrialStrain(tmpStrain);
  }


  //j-MNS (stress-strain, length=totalRubber/2)
  for (int i=0; i<(nMNS*nMNS); i++)  {
    double tmpCommitStrain = (ifCommit) ? commitStrnJMns[i] : theJNodeMNSMaterials[i]->getStrain();
    double tmpIncrDisp = incrDispij(6) + posLz[i]*incrDispij(10) - posLy[i]*incrDispij(11) - incrDispmn(3) - posLz[i]*incrDispmn(4) + posLy[i]*incrDispmn(5);
    double tmpStrain = tmpCommitStrain +  tmpIncrDisp/(totalRubber/2.0);
    theJNodeMNSMaterials[i]->setTrialStrain(tmpStrain);
  }
  
  //stiff springs in mid-height
  trialDspMidX  = (ifCommit) ? (commitDspMidX + dsxx) : (trialDspMidX + dsxx);
  trialDspMidRY = (ifCommit) ? (commitDspMidRY + dryy) : (trialDspMidRY + dryy);
  trialDspMidRZ = (ifCommit) ? (commitDspMidRZ + drzz) : (trialDspMidRZ + drzz);
  trialDspMidRX = (ifCommit) ? (commitDspMidRX + drxx) : (trialDspMidRX + drxx);

  return;
}


//calculate force components
void KikuchiBearing::subCalcFrcCpnt()
{

  //frcCpnt
  double fsy = 0.0 ;//(0)MSS Sum(fs*cos)
  double fsz = 0.0 ;//(1)    Sum(fs*sin)
  double fn  = 0.0 ;//(2)mid fn
  double fi  = 0.0 ;//(3)i-MNS Sum(fi)
  double fiy = 0.0 ;//(4)      Sum(fi*yi)
  double fiz = 0.0 ;//(5)      Sum(fi*xi)
  double fj  = 0.0 ;//(6)j-MNS Sum(fj)
  double fjy = 0.0 ;//(7)      Sum(fj*yj)
  double fjz = 0.0 ;//(8)      Sum(fj*xj)
  double fry = 0.0 ;//(9)mid my
  double frz = 0.0 ;//(10)   mz
  double frx = 0.0 ;//(11)   mx


  //MSS (force-deformation)
  for (int i=0; i<nMSS; i++)  {
    double tmpForce = theMidMSSMaterials[i]->getStress();
    fsy += tmpForce * cosTht[i];
    fsz += tmpForce * sinTht[i];
  }
  //equivalent coefficient for MSS
  fsy *= mssFeq;
  fsz *= mssFeq;


  //i-MNS (stress-strain, area=incA)
  for (int i=0; i<(nMNS*nMNS); i++)  {
    double tmpForce = theINodeMNSMaterials[i]->getStress()*incA*distFct[i];
    fi += tmpForce ;
    fiy += tmpForce * posLy[i];
    fiz += tmpForce * posLz[i];
  }

  //j-MNS (stress-strain, area=incA)
  for (int i=0; i<(nMNS*nMNS); i++)  {
    double tmpForce = theJNodeMNSMaterials[i]->getStress()*incA*distFct[i];
    fj += tmpForce ;
    fjy += tmpForce * posLy[i];
    fjz += tmpForce * posLz[i];
  }

  //stiff springs in mid-height
  trialFrcMidX  = stfMidX  * trialDspMidX;
  trialFrcMidRY = stfMidRY * trialDspMidRY;
  trialFrcMidRZ = stfMidRZ * trialDspMidRZ;
  trialFrcMidRX = stfMidRX * trialDspMidRX;

  fn  = trialFrcMidX;
  fry = trialFrcMidRY;
  frz = trialFrcMidRZ;
  frx = trialFrcMidRX;

  //frcCpnt output
  frcCpnt(0)  = fsy ;
  frcCpnt(1)  = fsz ;
  frcCpnt(2)  = fn  ;
  frcCpnt(3)  = fi  ;
  frcCpnt(4)  = fiy ;
  frcCpnt(5)  = fiz ;
  frcCpnt(6)  = fj  ;
  frcCpnt(7)  = fjy ;
  frcCpnt(8)  = fjz ;
  frcCpnt(9)  = fry ;
  frcCpnt(10) = frz ;
  frcCpnt(11) = frx ;
     
  return;
}



//calculate stiffness components
void KikuchiBearing::subCalcStfCpnt(){subCalcStfCpnt_main(false); return;}
void KikuchiBearing::subCalcStfCpntInit(){subCalcStfCpnt_main(true); return;}
void KikuchiBearing::subCalcStfCpnt_main(bool ifInit)
{

  //stfCpnt
  double kscc = 0.0 ;//(0)MSS Sum(ks*cos*cos)
  double kscs = 0.0 ;//(1)    Sum(ks*cos*sin)
  double ksss = 0.0 ;//(2)    Sum(ks*sin*sin)
  double kn   = 0.0 ;//(3)mid x(stiff)
  double ki   = 0.0 ;//(4)i-MNS Sum(ki)
  double kiy  = 0.0 ;//(5)      Sum(ki*yi)
  double kiz  = 0.0 ;//(6)      Sum(ki*zi)
  double kiyy = 0.0 ;//(7)      Sum(ki*yi*yi)
  double kiyz = 0.0 ;//(8)      Sum(ki*yi*zi)
  double kizz = 0.0 ;//(9)      Sum(ki*zi*zi)
  double kj   = 0.0 ;//(10)j-MNS Sum(kj)
  double kjy  = 0.0 ;//(11)      Sum(kj*yj)
  double kjz  = 0.0 ;//(12)      Sum(kj*zj)
  double kjyy = 0.0 ;//(13)      Sum(kj*yi*yj)
  double kjyz = 0.0 ;//(14)      Sum(kj*yj*zj)
  double kjzz = 0.0 ;//(15)      Sum(kj*zj*zj)
  double kry  = 0.0 ;//(16)mid ry(stiff)
  double krz  = 0.0 ;//(17)    rz(stiff)
  double krx  = 0.0 ;//(18)    rx(stiff)


  //MSS (force-deformation)
  for (int i=0; i<nMSS; i++)  {
    double tmpStiff = (ifInit) ? theMidMSSMaterials[i]->getInitialTangent() : theMidMSSMaterials[i]->getTangent();
    kscc += tmpStiff * cosTht[i] * cosTht[i];
    kscs += tmpStiff * cosTht[i] * sinTht[i];
    ksss += tmpStiff * sinTht[i] * sinTht[i];
  }
  //equivalent coefficient of MSS
  kscc *= mssSeq;
  kscs *= mssSeq;
  ksss *= mssSeq;
  

  //i-MNS (stress-strain, area=incA, length=totalRubber/2)
  for (int i=0; i<(nMNS*nMNS); i++)  {
    double tmpTangent = (ifInit) ? theINodeMNSMaterials[i]->getInitialTangent() : theINodeMNSMaterials[i]->getTangent();
    double tmpStiff = tmpTangent*incA*distFct[i]/(totalRubber/2.0);
    ki   += tmpStiff ;
    kiy  += tmpStiff * posLy[i] ;
    kiz  += tmpStiff * posLz[i] ;
    kiyy += tmpStiff * posLy[i] * posLy[i] ;
    kiyz += tmpStiff * posLy[i] * posLz[i] ;
    kizz += tmpStiff * posLz[i] * posLz[i] ;
  }


  //j-MNS (stress-strain, area=incA, length=totalRubber/2)
  for (int i=0; i<(nMNS*nMNS); i++)  {
    double tmpTangent = (ifInit) ? theJNodeMNSMaterials[i]->getInitialTangent() : theJNodeMNSMaterials[i]->getTangent();
    double tmpStiff = tmpTangent*incA*distFct[i]/(totalRubber/2.0);
    kj   += tmpStiff ;
    kjy  += tmpStiff * posLy[i] ;
    kjz  += tmpStiff * posLz[i] ;
    kjyy += tmpStiff * posLy[i] * posLy[i] ;
    kjyz += tmpStiff * posLy[i] * posLz[i] ;
    kjzz += tmpStiff * posLz[i] * posLz[i] ;
  }

  
  //stiff springs in mid-height
  kn  = stfMidX;
  kry = stfMidRY;
  krz = stfMidRZ;
  krx = stfMidRX;
  

  //stfCpnt output
  stfCpnt(0)  = kscc;
  stfCpnt(1)  = kscs;
  stfCpnt(2)  = ksss;
  stfCpnt(3)  = kn  ;
  stfCpnt(4)  = ki  ;
  stfCpnt(5)  = kiy ;
  stfCpnt(6)  = kiz ;
  stfCpnt(7)  = kiyy;
  stfCpnt(8)  = kiyz;
  stfCpnt(9)  = kizz;
  stfCpnt(10) = kj  ;
  stfCpnt(11) = kjy ;
  stfCpnt(12) = kjz ;
  stfCpnt(13) = kjyy;
  stfCpnt(14) = kjyz;
  stfCpnt(15) = kjzz;
  stfCpnt(16) = kry ;
  stfCpnt(17) = krz ;
  stfCpnt(18) = krx ;

  return;
}


//make K18 matrix (full)
void KikuchiBearing::subMakeKij18()
{
  //localLocalStiff
  static Matrix Kim(6,6)   ;//i-m
  static Matrix Kjn(6,6)   ;//j-n
  static Matrix Kmn(12,12) ;//m-n

  //stfCpnt
  double kscc = stfCpnt(0) ;//MSS Sum(ks*cos*cos)
  double kscs = stfCpnt(1) ;//    Sum(ks*cos*sin)
  double ksss = stfCpnt(2) ;//    Sum(ks*sin*sin)
  double kn   = stfCpnt(3) ;//mid x(stiff)
  double ki   = stfCpnt(4) ;//i-MNS Sum(ki)
  double kiy  = stfCpnt(5) ;//      Sum(ki*yi)
  double kiz  = stfCpnt(6) ;//      Sum(ki*zi)
  double kiyy = stfCpnt(7) ;//      Sum(ki*yi*yi)
  double kiyz = stfCpnt(8) ;//      Sum(ki*yi*zi)
  double kizz = stfCpnt(9) ;//      Sum(ki*zi*zi)
  double kj   = stfCpnt(10);//j-MNS Sum(kj)
  double kjy  = stfCpnt(11);//      Sum(kj*yj)
  double kjz  = stfCpnt(12);//      Sum(kj*zj)
  double kjyy = stfCpnt(13);//      Sum(kj*yi*yj)
  double kjyz = stfCpnt(14);//      Sum(kj*yj*zj)
  double kjzz = stfCpnt(15);//      Sum(kj*zj*zj)
  double kry  = stfCpnt(16);//mid ry (stiff)
  double krz  = stfCpnt(17);//    rz (stiff)
  double krx  = stfCpnt(18);//    rx (stiff)

  //dspCpnt
  double tym  = dspCpnt(0);
  double tzm  = dspCpnt(1);
  double tyn  = dspCpnt(2);
  double tzn  = dspCpnt(3);
  double dsy  = dspCpnt(4);
  double dsz  = dspCpnt(5);
  double hf   = dspCpnt(6);
  double dsy2 = dspCpnt(7);
  double dsz2 = dspCpnt(8);

  //Kmn
  Kmn( 0, 0) =  kn;
  Kmn( 0, 1) =  kn*tzm;
  Kmn( 0, 2) = -kn*tym;
  Kmn( 0, 3) =  0.0;
  Kmn( 0, 4) =  kn*dsz2;
  Kmn( 0, 5) = -kn*dsy2;
  Kmn( 0, 6) = -kn;
  Kmn( 0, 7) = -kn*tzn;
  Kmn( 0, 8) =  kn*tyn;
  Kmn( 0, 9) =  0.0;
  Kmn( 0,10) =  kn*dsz2;
  Kmn( 0,11) = -kn*dsy2;

  Kmn( 1, 1) =  kscc    + kn*tzm*tzm;
  Kmn( 1, 2) =  kscs    - kn*tym*tzm;
  Kmn( 1, 3) =  0.0;
  Kmn( 1, 4) = -kscs*hf + kn*dsz2*tzm;
  Kmn( 1, 5) =  kscc*hf - kn*dsy2*tzm;
  Kmn( 1, 6) =          - kn*tzm;
  Kmn( 1, 7) = -kscc    - kn*tzm*tzn;
  Kmn( 1, 8) = -kscs    + kn*tyn*tzm;
  Kmn( 1, 9) =  0.0;
  Kmn( 1,10) = -kscs*hf + kn*dsz2*tzm;
  Kmn( 1,11) =  kscc*hf - kn*dsy2*tzm;

  Kmn( 2, 2) =  ksss    + kn*tym*tym;
  Kmn( 2, 3) =  0.0;
  Kmn( 2, 4) = -ksss*hf - kn*dsz2*tym;
  Kmn( 2, 5) =  kscs*hf + kn*dsy2*tym;
  Kmn( 2, 6) =            kn*tym;
  Kmn( 2, 7) = -kscs    + kn*tym*tzn;
  Kmn( 2, 8) = -ksss    - kn*tym*tyn;
  Kmn( 2, 9) =  0.0;
  Kmn( 2,10) = -ksss*hf - kn*dsz2*tym;
  Kmn( 2,11) =  kscs*hf + kn*dsy2*tym;

  Kmn( 3, 3) =  krx;
  Kmn( 3, 4) =  0.0;
  Kmn( 3, 5) =  0.0;
  Kmn( 3, 6) =  0.0;
  Kmn( 3, 7) =  0.0;
  Kmn( 3, 8) =  0.0;
  Kmn( 3, 9) = -krx;
  Kmn( 3,10) =  0.0;
  Kmn( 3,11) =  0.0;

  Kmn( 4, 4) =  ksss*hf*hf + kn*dsz2*dsz2 + kry;
  Kmn( 4, 5) = -kscs*hf*hf - kn*dsz2*dsy2;
  Kmn( 4, 6) =             - kn*dsz2;
  Kmn( 4, 7) =  kscs*hf    - kn*dsz2*tzn;
  Kmn( 4, 8) =  ksss*hf    + kn*dsz2*tyn;
  Kmn( 4, 9) =  0.0;
  Kmn( 4,10) =  ksss*hf*hf + kn*dsz2*dsz2 - kry;
  Kmn( 4,11) = -kscs*hf*hf - kn*dsz2*dsy2;

  Kmn( 5, 5) =  kscc*hf*hf + kn*dsy2*dsy2 + krz;
  Kmn( 5, 6) =               kn*dsy2;
  Kmn( 5, 7) = -kscc*hf    + kn*dsy2*tzn;
  Kmn( 5, 8) = -kscs*hf    - kn*dsy2*tyn;
  Kmn( 5, 9) =  0.0;
  Kmn( 5,10) = -kscs*hf*hf - kn*dsz2*dsy2;
  Kmn( 5,11) =  kscc*hf*hf + kn*dsy2*dsy2 - krz;

  Kmn( 6, 6) =  kn;
  Kmn( 6, 7) =  kn*tzn;
  Kmn( 6, 8) = -kn*tyn;
  Kmn( 6, 9) =  0.0;
  Kmn( 6,10) =             - kn*dsz2;
  Kmn( 6,11) =               kn*dsy2;

  Kmn( 7, 7) =  kscc       + kn*tzn*tzn;
  Kmn( 7, 8) =  kscs       - kn*tyn*tzn;
  Kmn( 7, 9) =  0.0;
  Kmn( 7,10) =  kscs*hf    - kn*dsz2*tzn;
  Kmn( 7,11) = -kscc*hf    + kn*dsy2*tzn;

  Kmn( 8, 8) =  ksss       + kn*tyn*tyn;
  Kmn( 8, 9) =  0.0;
  Kmn( 8,10) =  ksss*hf    + kn*dsz2*tyn;
  Kmn( 8,11) = -kscs*hf    - kn*dsy2*tyn;

  Kmn( 9, 9) =  krx;
  Kmn( 9,10) =  0.0;
  Kmn( 9,11) =  0.0;

  Kmn(10,10) =  ksss*hf*hf + kn*dsz2*dsz2 + kry;
  Kmn(10,11) = -kscs*hf*hf - kn*dsz2*dsy2;

  Kmn(11,11) =  kscc*hf*hf + kn*dsy2*dsy2 + krz;

  for (int i=0; i<=10; i++)  { //symmetric
    for (int j=i+1; j<=11; j++)  {
      Kmn(j,i) = Kmn(i,j);
    }
  }

  //Kim
  Kim(0,0) =  ki;
  Kim(0,1) =  kiz;
  Kim(0,2) = -kiy;
  Kim(0,3) = -ki;
  Kim(0,4) = -kiz;
  Kim(0,5) =  kiy;

  Kim(1,1) =  kizz;
  Kim(1,2) = -kiyz;
  Kim(1,3) = -kiz;
  Kim(1,4) = -kizz;
  Kim(1,5) =  kiyz;

  Kim(2,2) =  kiyy;
  Kim(2,3) =  kiy;
  Kim(2,4) =  kiyz;
  Kim(2,5) = -kiyy;

  Kim(3,3) =  ki;
  Kim(3,4) =  kiz;
  Kim(3,5) = -kiy;

  Kim(4,4) =  kizz;
  Kim(4,5) = -kiyz;

  Kim(5,5) =  kiyy;

  for (int i=0; i<=4; i++)  { //symmetric
    for (int j=i+1; j<=5; j++)  {
      Kim(j,i) = Kim(i,j);
    }
  }

  //Kjn
  Kjn(0,0) =  kj;
  Kjn(0,1) =  kjz;
  Kjn(0,2) = -kjy;
  Kjn(0,3) = -kj;
  Kjn(0,4) = -kjz;
  Kjn(0,5) =  kjy;

  Kjn(1,1) =  kjzz;
  Kjn(1,2) = -kjyz;
  Kjn(1,3) = -kjz;
  Kjn(1,4) = -kjzz;
  Kjn(1,5) =  kjyz;

  Kjn(2,2) =  kjyy;
  Kjn(2,3) =  kjy;
  Kjn(2,4) =  kjyz;
  Kjn(2,5) = -kjyy;

  Kjn(3,3) =  kj;
  Kjn(3,4) =  kjz;
  Kjn(3,5) = -kjy;

  Kjn(4,4) =  kjzz;
  Kjn(4,5) = -kjyz;

  Kjn(5,5) =  kjyy;

  for (int i=0; i<=4; i++)  { //symmetric
    for (int j=i+1; j<=5; j++)  {
      Kjn(j,i) = Kjn(i,j);
    }
  }


  //Kij18
  Kij18( 0, 0) =              Kim( 0, 0);
  Kij18( 0, 1) =              0.0;
  Kij18( 0, 2) =              0.0;
  Kij18( 0, 3) =              0.0;
  Kij18( 0, 4) =              Kim( 0, 1);
  Kij18( 0, 5) =              Kim( 0, 2);
  Kij18( 0, 6) =              0.0;
  Kij18( 0, 7) =              0.0;
  Kij18( 0, 8) =              0.0;
  Kij18( 0, 9) =              0.0;
  Kij18( 0,10) =              0.0;
  Kij18( 0,11) =              0.0;
  Kij18( 0,12) =              Kim( 0, 3);
  Kij18( 0,13) =              Kim( 0, 4);
  Kij18( 0,14) =              Kim( 0, 5);
  Kij18( 0,15) =              0.0;
  Kij18( 0,16) =              0.0;
  Kij18( 0,17) =              0.0;

  Kij18( 1, 1) = Kmn( 1, 1);
  Kij18( 1, 2) = Kmn( 1, 2);
  Kij18( 1, 3) = Kmn( 1, 3);
  Kij18( 1, 4) = 0.0;
  Kij18( 1, 5) = 0.0;
  Kij18( 1, 6) = 0.0;
  Kij18( 1, 7) = Kmn( 1, 7);
  Kij18( 1, 8) = Kmn( 1, 8);
  Kij18( 1, 9) = Kmn( 1, 9);
  Kij18( 1,10) = 0.0;
  Kij18( 1,11) = 0.0;
  Kij18( 1,12) = Kmn( 1, 0);
  Kij18( 1,13) = Kmn( 1, 4);
  Kij18( 1,14) = Kmn( 1, 5);
  Kij18( 1,15) = Kmn( 1, 6);
  Kij18( 1,16) = Kmn( 1,10);
  Kij18( 1,17) = Kmn( 1,11);

  Kij18( 2, 2) = Kmn( 2, 2);
  Kij18( 2, 3) = Kmn( 2, 3);
  Kij18( 2, 4) = 0.0;
  Kij18( 2, 5) = 0.0;
  Kij18( 2, 6) = 0.0;
  Kij18( 2, 7) = Kmn( 2, 7);
  Kij18( 2, 8) = Kmn( 2, 8);
  Kij18( 2, 9) = Kmn( 2, 9);
  Kij18( 2,10) = 0.0;
  Kij18( 2,11) = 0.0;
  Kij18( 2,12) = Kmn( 2, 0);
  Kij18( 2,13) = Kmn( 2, 4);
  Kij18( 2,14) = Kmn( 2, 5);
  Kij18( 2,15) = Kmn( 2, 6);
  Kij18( 2,16) = Kmn( 2,10);
  Kij18( 2,17) = Kmn( 2,11);

  Kij18( 3, 3) = Kmn( 3, 3);
  Kij18( 3, 4) = 0.0;
  Kij18( 3, 5) = 0.0;
  Kij18( 3, 6) = 0.0;
  Kij18( 3, 7) = Kmn( 3, 7);
  Kij18( 3, 8) = Kmn( 3, 8);
  Kij18( 3, 9) = Kmn( 3, 9);
  Kij18( 3,10) = 0.0;
  Kij18( 3,11) = 0.0;
  Kij18( 3,12) = Kmn( 3, 0);
  Kij18( 3,13) = Kmn( 3, 4);
  Kij18( 3,14) = Kmn( 3, 5);
  Kij18( 3,15) = Kmn( 3, 6);
  Kij18( 3,16) = Kmn( 3,10);
  Kij18( 3,17) = Kmn( 3,11);

  Kij18( 4, 4) =              Kim( 1, 1);
  Kij18( 4, 5) =              Kim( 1, 2);
  Kij18( 4, 6) =              0.0;
  Kij18( 4, 7) =              0.0;
  Kij18( 4, 8) =              0.0;
  Kij18( 4, 9) =              0.0;
  Kij18( 4,10) =              0.0;
  Kij18( 4,11) =              0.0;
  Kij18( 4,12) =              Kim( 1, 3);
  Kij18( 4,13) =              Kim( 1, 4);
  Kij18( 4,14) =              Kim( 1, 5);
  Kij18( 4,15) =              0.0;
  Kij18( 4,16) =              0.0;
  Kij18( 4,17) =              0.0;

  Kij18( 5, 5) =              Kim( 2, 2);
  Kij18( 5, 6) =              0.0;
  Kij18( 5, 7) =              0.0;
  Kij18( 5, 8) =              0.0;
  Kij18( 5, 9) =              0.0;
  Kij18( 5,10) =              0.0;
  Kij18( 5,11) =              0.0;
  Kij18( 5,12) =              Kim( 2, 3);
  Kij18( 5,13) =              Kim( 2, 4);
  Kij18( 5,14) =              Kim( 2, 5);
  Kij18( 5,15) =              0.0;
  Kij18( 5,16) =              0.0;
  Kij18( 5,17) =              0.0;
  
  Kij18( 6, 6) =              Kjn( 0, 0);
  Kij18( 6, 7) =              0.0;
  Kij18( 6, 8) =              0.0;
  Kij18( 6, 9) =              0.0;
  Kij18( 6,10) =              Kjn( 0, 1);
  Kij18( 6,11) =              Kjn( 0, 2);
  Kij18( 6,12) =              0.0;
  Kij18( 6,13) =              0.0;
  Kij18( 6,14) =              0.0;
  Kij18( 6,15) =              Kjn( 0, 3);
  Kij18( 6,16) =              Kjn( 0, 4);
  Kij18( 6,17) =              Kjn( 0, 5);

  Kij18( 7, 7) = Kmn( 7, 7);
  Kij18( 7, 8) = Kmn( 7, 8);
  Kij18( 7, 9) = Kmn( 7, 9);
  Kij18( 7,10) = 0.0;
  Kij18( 7,11) = 0.0;
  Kij18( 7,12) = Kmn( 7, 0);
  Kij18( 7,13) = Kmn( 7, 4);
  Kij18( 7,14) = Kmn( 7, 5);
  Kij18( 7,15) = Kmn( 7, 6);
  Kij18( 7,16) = Kmn( 7,10);
  Kij18( 7,17) = Kmn( 7,11);

  Kij18( 8, 8) = Kmn( 8, 8);
  Kij18( 8, 9) = Kmn( 8, 9);
  Kij18( 8,10) = 0.0;
  Kij18( 8,11) = 0.0;
  Kij18( 8,12) = Kmn( 8, 0);
  Kij18( 8,13) = Kmn( 8, 4);
  Kij18( 8,14) = Kmn( 8, 5);
  Kij18( 8,15) = Kmn( 8, 6);
  Kij18( 8,16) = Kmn( 8,10);
  Kij18( 8,17) = Kmn( 8,11);

  Kij18( 9, 9) = Kmn( 9, 9);
  Kij18( 9,10) = 0.0;
  Kij18( 9,11) = 0.0;
  Kij18( 9,12) = Kmn( 9, 0);
  Kij18( 9,13) = Kmn( 9, 4);
  Kij18( 9,14) = Kmn( 9, 5);
  Kij18( 9,15) = Kmn( 9, 6);
  Kij18( 9,16) = Kmn( 9,10);
  Kij18( 9,17) = Kmn( 9,11);

  Kij18(10,10) =              Kjn( 1, 1);
  Kij18(10,11) =              Kjn( 1, 2);
  Kij18(10,12) =              0.0;
  Kij18(10,13) =              0.0;
  Kij18(10,14) =              0.0;
  Kij18(10,15) =              Kjn( 1, 3);
  Kij18(10,16) =              Kjn( 1, 4);
  Kij18(10,17) =              Kjn( 1, 5);

  Kij18(11,11) =              Kjn( 2, 2);
  Kij18(11,12) =              0.0;
  Kij18(11,13) =              0.0;
  Kij18(11,14) =              0.0;
  Kij18(11,15) =              Kjn( 2, 3);
  Kij18(11,16) =              Kjn( 2, 4);
  Kij18(11,17) =              Kjn( 2, 5);

  Kij18(12,12) = Kmn( 0, 0) + Kim( 3, 3);
  Kij18(12,13) = Kmn( 0, 4) + Kim( 3, 4);
  Kij18(12,14) = Kmn( 0, 5) + Kim( 3, 5);
  Kij18(12,15) = Kmn( 0, 6);
  Kij18(12,16) = Kmn( 0,10);
  Kij18(12,17) = Kmn( 0,11);

  Kij18(13,13) = Kmn( 4, 4) + Kim( 4, 4);
  Kij18(13,14) = Kmn( 4, 5) + Kim( 4, 5);
  Kij18(13,15) = Kmn( 4, 6);
  Kij18(13,16) = Kmn( 4,10);
  Kij18(13,17) = Kmn( 4,11);

  Kij18(14,14) = Kmn( 5, 5) + Kim( 5, 5);
  Kij18(14,15) = Kmn( 5, 6);
  Kij18(14,16) = Kmn( 5,10);
  Kij18(14,17) = Kmn( 5,11);

  Kij18(15,15) = Kmn( 6, 6) + Kjn( 3, 3);
  Kij18(15,16) = Kmn( 6,10) + Kjn( 3, 4);
  Kij18(15,17) = Kmn( 6,11) + Kjn( 3, 5);

  Kij18(16,16) = Kmn(10,10) + Kjn( 4, 4);
  Kij18(16,17) = Kmn(10,11) + Kjn( 4, 5);

  Kij18(17,17) = Kmn(11,11) + Kjn( 5, 5);

  for (int i=0; i<=16; i++)  { //symmetric
    for (int j=i+1; j<=17; j++)  {
      Kij18(j,i) = Kij18(i,j);
    }
  }

  return;
}


//make submatrices
void KikuchiBearing::subSubmatKij18()
{

  for (int i=0; i<=11; i++)  {
    for (int j=0; j<=11; j++)  {
      Kij18_11(i,j) = Kij18(i,j);
    }
  }

  for (int i=0; i<=11; i++)  {
    for (int j=0; j<=5; j++)  {
      Kij18_12(i,j) = Kij18(i,j+12);
    }
  }

  for (int i=0; i<=5; i++)  {
    for (int j=0; j<=11; j++)  {
      Kij18_21(i,j) = Kij18(i+12,j);
    }
  }

  for (int i=0; i<=5; i++)  {
    for (int j=0; j<=5; j++)  {
      Kij18_22(i,j) = Kij18(i+12,j+12);
    }
  }
  
  return;
}

//reduct K18 matrix
void KikuchiBearing::subReductKij()
{
  //reduction
  subSubmatKij18(); //make submatrices
  Kij18_22.Invert(invKij18_22); // invKij18_22 = inv(Kij18_22)
  Kij = Kij18_11 - Kij18_12*invKij18_22*Kij18_21;

  return;
}

//calculate Fij and Fmn
void KikuchiBearing::subMakeFijFmn()
{

  //localLocalForce
  static Vector fmn0(12)  ;//m-n
  static Vector fim0(6)  ;//i-m
  static Vector fjn0(6)  ;//j-n
  
  //frcCpnt
  double fsy = frcCpnt(0);//MSS Sum(fs*cos)
  double fsz = frcCpnt(1);//    Sum(fs*sin)
  double fn  = frcCpnt(2);//mid fn
  double fi  = frcCpnt(3);//i-MNS Sum(fi)
  double fiy = frcCpnt(4);//      Sum(fi*yi)
  double fiz = frcCpnt(5);//      Sum(fi*xi)
  double fj  = frcCpnt(6);//j-MNS Sum(fj)
  double fjy = frcCpnt(7);//      Sum(fj*yj)
  double fjz = frcCpnt(8);//      Sum(fj*xj)
  double fry = frcCpnt(9);//mid my
  double frz = frcCpnt(10);//   mz
  double frx = frcCpnt(11);//   mx

  //dspCpnt
  double tym  = dspCpnt(0);
  double tzm  = dspCpnt(1);
  double tyn  = dspCpnt(2);
  double tzn  = dspCpnt(3);
  double dsy  = dspCpnt(4);
  double dsz  = dspCpnt(5);
  double hf   = dspCpnt(6);
  double dsy2 = dspCpnt(7);
  double dsz2 = dspCpnt(8);


  //fmn0
  fmn0( 0) = -fn                     ;// fxm
  fmn0( 1) = -fn*tzm  - fsy          ;// fym
  fmn0( 2) =  fn*tym  - fsz          ;// fzm
  fmn0( 3) = -frx                    ;// mxm
  fmn0( 4) = -fn*dsz2 + fsz*hf - fry ;// mym
  fmn0( 5) =  fn*dsy2 - fsy*hf - frz ;// mzm
  fmn0( 6) =  fn                     ;// fxn
  fmn0( 7) =  fn*tzn  + fsy          ;// fyn
  fmn0( 8) = -fn*tyn  + fsz          ;// fzn
  fmn0( 9) =  frx                    ;// mxn
  fmn0(10) = -fn*dsz2 + fsz*hf + fry ;// myn
  fmn0(11) =  fn*dsy2 - fsy*hf + frz ;// mzn

  //fim0
  fim0(0) = -fi    ;// fxi
  fim0(1) = -fiz   ;// myi
  fim0(2) =  fiy   ;// mzi
  fim0(3) =  fi    ;// fxm
  fim0(4) =  fiz   ;// mym
  fim0(5) = -fiy   ;// mzm
  
  //fjn0
  fjn0(0) =  fj    ;// fxj
  fjn0(1) =  fjz   ;// myj
  fjn0(2) = -fjy   ;// mzj
  fjn0(3) = -fj    ;// fxn
  fjn0(4) = -fjz   ;// myn
  fjn0(5) =  fjy   ;// mzn

  //internal nodes
  Fmn(0) =  fim0(3) + fmn0( 0)              ;// fxm
  Fmn(1) =  fim0(4) + fmn0( 4)              ;// mym
  Fmn(2) =  fim0(5) + fmn0( 5)              ;// mzm
  Fmn(3) =            fmn0( 6) + fjn0(3)    ;// fxn
  Fmn(4) =            fmn0(10) + fjn0(4)    ;// myn
  Fmn(5) =            fmn0(11) + fjn0(5)    ;// mzn
  
  //external nodes
  Fij( 0) =  fim0(0)                        ;// fxi
  Fij( 1) =           fmn0( 1)              ;// fyi(=fym)
  Fij( 2) =           fmn0( 2)              ;// fzi(=fzm)
  Fij( 3) =           fmn0( 3)              ;// mxi(=mxm)
  Fij( 4) =  fim0(1)                        ;// myi
  Fij( 5) =  fim0(2)                        ;// mzi
  Fij( 6) =                      fjn0(0)    ;// fxj
  Fij( 7) =           fmn0( 7)              ;// fyj(=fyn)
  Fij( 8) =           fmn0( 8)              ;// fzj(=fzn)
  Fij( 9) =           fmn0( 9)              ;// mxj(=mxn)
  Fij(10) =                      fjn0(1)    ;// myj
  Fij(11) =                      fjn0(2)    ;// mzj
 
  return;
}

//calculate Feq and Seq
void KikuchiBearing::subCalcMSSFeqSeq()
{

    //equivalent coefficient for force and stiffness
    if (limDisp >= 0) {
      double uRef, fRef, sRef;//deformation, force, stiffness
      double uCmp, fSum, sSum;
      double refDisp;
      
      //reference disp
      refDisp = sqrt(basicDisp(1)*basicDisp(1)+basicDisp(2)*basicDisp(2));
      uRef = (refDisp>limDisp) ? refDisp : limDisp;

      //material to calculate Feq and Seq
      dmyMSSMaterial->setTrialStrain(uRef,0);
      fRef = dmyMSSMaterial->getStress();
      sRef = dmyMSSMaterial->getTangent();
      
      //total force, total stiffness
      fSum = 0.0;
      sSum = 0.0;
      for (int i=0; i<nMSS; i++)  {
        uCmp = uRef * cosTht[i];
        dmyMSSMaterial->setTrialStrain(uCmp,0);
    	fSum += dmyMSSMaterial->getStress() * cosTht[i];
    	sSum += dmyMSSMaterial->getTangent() * cosTht[i] * cosTht[i];
      }

      //Feq, Seq
      mssFeq = fRef/fSum;
      mssSeq = sRef/sSum;
      
    } else {

      mssFeq = 1.0;
      mssSeq = 1.0;
      
    }
    
}


