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
                                                                                                                                                                                           
// Written: vavgen  (Evangelos Avgenakis)
// Some preexisting code from basic OpenSees elements was used
// Created: July 2020
// Revision: A
//
// Description: This element can be used to describe the rocking motion of 2d deformable bodies,
// 		either elastic or inelastic, under static or dynamic loading. Apart from the deformability
// 		along the length of the member, the element is able to account for deformability near the
// 		contact area, where nonlinear stress distributions develop and sections do not remain plane.
// 		Furthermore, the element is able to account for constraints along the length of the rocking
// 		member imposed by other structural members, as well as sliding and upthrow.
//
// References:
// 	1. Avgenakis E. and Psycharis I.N. (2017) “Modeling of Rocking Elastic Flexible Bodies under Static
//  		Loading Considering the Nonlinear Stress Distribution at Their Base.” Journal of Structural
//		Engineering 143(7): 04017051.
//	2. Avgenakis, E. and Psycharis, I. N. (2019) “Determination of the nonlinear displacement distribution
//		of the semi-infinite strip–Application to deformable rocking bodies.” International Journal
//		of Solids and Structures, 170, 22-37.
//	3. Avgenakis E. and Psycharis I.N. (2020) “Modeling of inelastic rocking bodies under cyclic loading.”
//		Journal of Engineering Mechanics 146(4): 04020020.
// 	4. Avgenakis E. and Psycharis I.N. (2020) “An integrated macroelement formulation for the dynamic
//		response of inelastic deformable rocking bodies.” Earthquake Engineering and Structural Dynamics.

#include "RockingBC.h"
#include <ElementalLoad.h>
#include <elementAPI.h>

#include <Domain.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>

#include <CrdTransf.h>
#include <Information.h>
#include <Parameter.h>
#include <ElementResponse.h>
#include <Renderer.h>

#include <elementAPI.h>
#include <G3Globals.h>

#include <Node.h>
#include <Message.h>
#include <UniaxialMaterial.h>

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <fstream>

#ifdef _USRDLL
#define ELE_TAG_RockingBC                 199
#define RockingBC_dll 1
#include "Eigen/Dense"
#include <windows.h>
#define OPS_Export extern "C" _declspec(dllexport)
#elif _MACOSX
#define OPS_Export extern "C" __attribute__((visibility("default")))
#else
#define OPS_Export extern "C"
#define RockingBC_dll 0
#endif

//OPS_Export void* OPS_RockingBC() //CHANGE FOR DLL CREATION
void* OPS_RockingBC()
{
    if(OPS_GetNumRemainingInputArgs() < 10) {
	opserr<<"Insufficient mandatory arguments: eleTag,iNode,jNode,Nw,E,nu,sy,B,w,mu; Optional arguments: convlim,maxtries,af,aflim,convlimmult,usecomstiff,useshear,blevery\n";
	return 0;
    }

    int ndm = OPS_GetNDM();
    int ndf = OPS_GetNDF();
    if(ndm != 2 || ndf != 3) {
	opserr<<"ndm must be 2 and ndf must be 3\n";
	return 0;
    }

    // inputs: 
    int iData[4];
    int numData = 4;
    if(OPS_GetIntInput(&numData,&iData[0]) < 0) return 0;

    double data[6];
	numData = 6;
    if(OPS_GetDoubleInput(&numData,&data[0]) < 0) return 0;

	// options
	double convlim = 1.0e-14;
	int maxtries = 100;
	double af = 1.0;
	double aflim = 0.4;
	double convlimmult = 1.0;
	int usecomstiff = 0; // Me -1 xrisimopoiei panta kai to Wcommit
	int useshear = 0;
	int blevery = 1;
	double NlimN = 0.1;
	double NlimT = 10.0;
	double Dtlim = 1.0e-8;
	int errorifNexceeds = 0;
	int useUelNM = 1;

	numData = 1;
	while (OPS_GetNumRemainingInputArgs() > 0) {
		std::string type = OPS_GetString();
		if (type == "-convlim") {
			if (OPS_GetNumRemainingInputArgs() > 0) {
				if (OPS_GetDoubleInput(&numData, &convlim) < 0) return 0;
			}
		}
		else if (type == "-maxtries") {
			if (OPS_GetNumRemainingInputArgs() > 0) {
				if (OPS_GetIntInput(&numData, &maxtries) < 0) return 0;
			}
		}
		else if (type == "-af") {
			if (OPS_GetNumRemainingInputArgs() > 0) {
				if (OPS_GetDoubleInput(&numData, &af) < 0) return 0;
			}
		}
		else if (type == "-aflim") {
			if (OPS_GetNumRemainingInputArgs() > 0) {
				if (OPS_GetDoubleInput(&numData, &aflim) < 0) return 0;
			}
		}
		else if (type == "-convlimmult") {
			if (OPS_GetNumRemainingInputArgs() > 0) {
				if (OPS_GetDoubleInput(&numData, &convlimmult) < 0) return 0;
			}
		}
		else if (type == "-usecomstiff") {
			if (OPS_GetNumRemainingInputArgs() > 0) {
				if (OPS_GetIntInput(&numData, &usecomstiff) < 0) return 0;
			}
		}
		else if (type == "-useshear") {
			if (OPS_GetNumRemainingInputArgs() > 0) {
				if (OPS_GetIntInput(&numData, &useshear) < 0) return 0;
			}
		}
		else if (type == "-blevery") {
			if (OPS_GetNumRemainingInputArgs() > 0) {
				if (OPS_GetIntInput(&numData, &blevery) < 0) return 0;
			}
		}
		else if (type == "-NlimN") {
			if (OPS_GetNumRemainingInputArgs() > 0) {
				if (OPS_GetDoubleInput(&numData, &NlimN) < 0) return 0;
			}
		}
		else if (type == "-NlimT") {
			if (OPS_GetNumRemainingInputArgs() > 0) {
				if (OPS_GetDoubleInput(&numData, &NlimT) < 0) return 0;
			}
		}
		else if (type == "-Dtlim") {
			if (OPS_GetNumRemainingInputArgs() > 0) {
				if (OPS_GetDoubleInput(&numData, &Dtlim) < 0) return 0;
			}
		}
		else if (type == "-errorifNexceeds") {
			if (OPS_GetNumRemainingInputArgs() > 0) {
				if (OPS_GetIntInput(&numData, &errorifNexceeds) < 0) return 0;
			}
		}
		else if (type == "-useUelNM") {
			if (OPS_GetNumRemainingInputArgs() > 0) {
				if (OPS_GetIntInput(&numData, &useUelNM) < 0) return 0;
			}
		}
	}

	if (aflim > af) { aflim = af; }

    return new RockingBC(iData[0], iData[1], iData[2], iData[3], data[0],data[1],data[2], data[3],
		data[4], data[5], convlim, maxtries, af, aflim,
		convlimmult, usecomstiff, useshear, blevery, NlimN, NlimT, Dtlim, errorifNexceeds, useUelNM);
}

RockingBC::RockingBC()
  :Element(0,ELE_TAG_RockingBC), 
  E(0.0), nu(0.0), sy(0.0), B(0.0), w(0.0), alpha(0.0), mu(0.0), Fe(6),
	connectedExternalNodes(2),
	cosTheta(0), sinTheta(0),
	nodeIPtr(0), nodeJPtr(0), L(0), ue(6), uecommit(6), uepr(6), Nw(2),
	convlim(0), maxtries(0), af(0), aflim(0), convlimmult(0), usecomstiff(0), useshear(0), blevery(0),
	NlimN(0), NlimT(0), Dtlim(0), errorifNexceeds(0), useUelNM(0)
{
  // set node pointers to NULL
  for (int i=0; i<2; i++)
    theNodes[i] = 0;     

  // Parameters
  EI = 0.0;
  EA = 0.0;
  G = 0.0;
  GA = 0.0;
}

RockingBC::RockingBC(int tag, int Nd1, int Nd2, int nw,
				 double e, double Nu, double Sy, double bb, double ww, double Mu,
				double Convlim, int Maxtries, double Af, double Aflim,
	double Convlimmult, int Usecomstiff, int Useshear, int Blevery,
	double NLimN, double NLimT, double DTlim, int ErrorifNexceeds, int UseUelNM)
  :Element(tag,ELE_TAG_RockingBC), 
  E(e), nu(Nu), sy(Sy), B(bb), w(ww), alpha(1.2), Fe(6),
  connectedExternalNodes(2),
  cosTheta(0), sinTheta(0),
  nodeIPtr(0), nodeJPtr(0), L(0), ue(6), uecommit(6), uepr(6), Nw(nw), mu(Mu),
  convlim(Convlim), maxtries(Maxtries), af(Af), aflim(Aflim),
	convlimmult(Convlimmult), usecomstiff(Usecomstiff), useshear(Useshear), blevery(Blevery),
	NlimN(NLimN), NlimT(NLimT), Dtlim(DTlim), errorifNexceeds(ErrorifNexceeds),
	useUelNM(UseUelNM)
{
  connectedExternalNodes(0) = Nd1;
  connectedExternalNodes(1) = Nd2;
   
  tagtag = tag;
 
  //theCoordTransf = coordTransf.getCopy2d();
  //  
  //if (!theCoordTransf) {
  //  opserr << "RockingBC::RockingBC -- failed to get copy of coordinate transformation\n";
  //  exit(01);
  //}

  // set node pointers to NULL
  theNodes[0] = 0;
  theNodes[1] = 0;

  // Parameters
  b = B / 2.0;
  I = w*B*B*B / 12.0;
  A = w*B;
  EI = E*I;
  EA = E*A;
  G = E / 2.0 /(1.+nu);
  GA = G*A;

  if (Nw < 2) { opserr << "Nw must be larger than 2"; exit(-1); }
  if (useshear) { noshear = false; }

  Yw = Vector(Nw);
  double yexp = 1.0;

  for (int i = 0; i != Nw/2+1; i++)
  {
	  Yw(i) = -1.0 + pow(2.0*i/(Nw-1),yexp);
  }
  for (int i = Nw / 2 - 1; i != -1; i--)
  {
	  Yw(Nw-i-1) = 1.0 - pow(2.0*i/(Nw - 1), yexp);
  }
  //std::cout << Yw << std::endl;

  UNM_calc(Yw, UN, UM);

  W = Vector(Nw);
  Winit = Vector(Nw);
  Wpr = Vector(Nw);
  Wcommit = Vector(Nw);

  Youter(0) = -1.;
  Youter(1) = 1.;
  
  ey = sy / E;

  for (size_t i = 0; i != Nw - 1; i++) {
	  Vec vvv1{ Yw[i],Yw[i + 1] };
	  Vec vvv2{ 0.0, 0.0 };
	  Vecint vvv3{ 0 };
	  Matrix mmm4{};
	  Ysi.push_back(vvv1);
	  Ysi_com.push_back(vvv1);
	  Yupi.push_back(vvv1);
	  Yupi_com.push_back(vvv1);
	  Si.push_back(vvv2);
	  Si_com.push_back(vvv2);
	  Upi.push_back(vvv2);
	  Upi_com.push_back(vvv2);
	  Ua_pos.push_back(vvv2);
	  Ys_cats.push_back(vvv3);
	  dYsi_dW.push_back(mmm4);
	  dSi_dW.push_back(mmm4);
  }

  Ec = Vector(Nw);
  El = Vector(Nw);
  
  Nints = Vector(Nw-1);
  Mints = Vector(Nw-1);
  dNints_dW = Matrix(Nw-1,Nw);
  dMints_dW = Matrix(Nw-1,Nw);

  Ys = Yw;
  Ys_com = Yw;
  Yup = Yw;
  Yup_com = Yw;
  Up = Vector(Nw);
  Up_com = Vector(Nw);
  S = Vector(Nw);
  S_com = Vector(Nw);
  Ks = Vector(Nw - 1);
  Ks_com = Vector(Nw - 1);
  Kup = Vector(Nw - 1);
  Kup_com = Vector(Nw - 1);

  Ud = Vector(Nw);
  dUd_dW = Matrix(Nw, Nw);
  dUd_due = Matrix(Nw, 6);
  dW_due = Matrix(Nw, 6);
  dW_due_com = Matrix(Nw, 6);
  dW_due_pr = Matrix(Nw, 6);

  dun_dW = Matrix(3, Nw);
  dues_dW = Matrix(6, Nw);
  dsL_dW = Vector(Nw);

  Upl = Vector(Nw);
  Ua = Vector(Nw);
  dUa_dW = Matrix(Nw, Nw);
  dYouter_dW = Matrix(2, Nw);
  dN_dW = Vector(Nw);
  dM_dW = Vector(Nw);
  dQ_dW = Vector(Nw);
  dt_dW = Vector(Nw);
  dut_dW = Matrix(2, Nw);
  durf_dW = Matrix(2, Nw);
  Uel = Vector(Nw);
  dUel_dW = Matrix(Nw,Nw);

  dFn2_dW = Matrix(2, Nw);
  dFn_dW = Matrix(3, Nw);
  dFntot_dW = Matrix(3, Nw);
  dFnD_dW = Matrix(3, Nw);
  dFes_dW = Matrix(6, Nw);
  dFeV_dW = Matrix(6, Nw);

  dlim1_dW = Vector(Nw);
  dlim2_dW = Vector(Nw);

  due5_due(5) = 1.;

  DW = Vector(Nw);
  DWcommit = Vector(Nw);
  dUd_dW = Matrix(Nw, Nw);
  dUd_due = Matrix(Nw, 6);

  Uel_com = Vector(Nw);
  UB = Matrix(Nw, 0);
  dUB_dR = Matrix(Nw, 0);
  UB_R.clear();

  Im1 = Vector(Nw);
  Jm1 = Vector(Nw);
  Im1_calc(Yw, Im1);
  Jm1_calc(Yw, Jm1);

  dgQ_dW = Vector(Nw);
  dth2_dW = Vector(Nw);
  
}

RockingBC::~RockingBC()
{
}

int
RockingBC::getNumExternalNodes(void) const
{
    return 2;
}

const ID &
RockingBC::getExternalNodes(void) 
{
    return connectedExternalNodes;
}

Node **
RockingBC::getNodePtrs(void) 
{
  return theNodes;
}

int
RockingBC::getNumDOF(void)
{
    return 6;
}

void
RockingBC::setDomain(Domain *theDomain)
{
  if (theDomain == 0) {
    opserr << "RockingBC::setDomain -- Domain is null\n";
    exit(-1);
  }
    
    theNodes[0] = theDomain->getNode(connectedExternalNodes(0));
    theNodes[1] = theDomain->getNode(connectedExternalNodes(1));    
    
    if (theNodes[0] == 0) {
      opserr << "RockingBC::setDomain -- Node 1: " << connectedExternalNodes(0) << " does not exist\n";
      exit(-1);
    }
			      
    if (theNodes[1] == 0) {
      opserr << "RockingBC::setDomain -- Node 2: " << connectedExternalNodes(1) << " does not exist\n";
      exit(-1);
    }

    int dofNd1 = theNodes[0]->getNumberDOF();
    int dofNd2 = theNodes[1]->getNumberDOF();    
    
    if (dofNd1 != 3) {
      opserr << "RockingBC::setDomain -- Node 1: " << connectedExternalNodes(0) 
	     << " has incorrect number of DOF\n";
      exit(-1);
    }
    
    if (dofNd2 != 3) {
      opserr << "RockingBC::setDomain -- Node 2: " << connectedExternalNodes(1) 
	     << " has incorrect number of DOF\n";
      exit(-1);
    }
	
    this->DomainComponent::setDomain(theDomain);
    
 //   if (theCoordTransf->initialize(theNodes[0], theNodes[1]) != 0) {
	//opserr << "RockingBC::setDomain -- Error initializing coordinate transformation\n";
	//exit(-1);
 //   }
	this->initialize(theNodes[0], theNodes[1]);
    
    double L = this->getInitialLength();

    if (L == 0.0) {
      opserr << "RockingBC::setDomain -- Element has zero length\n";
      exit(-1);
    }

}

int
RockingBC::initialize(Node *nodeIPointer, Node *nodeJPointer)
{
	int error;

	nodeIPtr = nodeIPointer;
	nodeJPtr = nodeJPointer;

	if ((!nodeIPtr) || (!nodeJPtr))
	{
		opserr << "\nRockingBC::initialize";
		opserr << "\ninvalid pointers to the element nodes\n";
		return -1;
	}

	// get element length and orientation
	if ((error = this->compElemtLengthAndOrient()))
		return error;

	fr_calc();

	k1 = 1.0 / fr(2, 2);
	k2 = fr(2, 1) / fr(2, 2);
	frr(0, 0) = fr(0, 0);
	frr(1, 1) = fr(1, 1) - fr(1, 2) * fr(2, 1) / fr(2, 2);

	dw1_due = Vector(6); dw1_due(2) = 1.0;
	dr_due = Vector(6); dr_due(4) = 1.0 / L; dr_due(1) = -1.0 / L;
	dw2_due = Vector(6); dw2_due(5) = 1.0;

	durth_due = Matrix(2, 6);
	durth_dW = Matrix(2, Nw);

	CC = Matrix(Nw, 2);
	for (int i = 0; i != Nw; i++) {
		CC(i, 0) = 1.0;
		CC(i, 1) = Yw(i);
	}

	BB(0, 0) = 1. / b;
	BB(1, 1) = 1.;

	CB = CC*BB;

	dutar_due = Matrix(2, 6);
	dutar_dW = Matrix(2, Nw);

	TF1(0, 0) = 1;
	TF1(1, 1) = b;
	TF1(2, 1) = -b;
	TF1(2, 2) = -L;
	TF = b*w*TF1*E;

	return 0;
}

int
RockingBC::compElemtLengthAndOrient(void)
{
	// element projection
	static Vector dx(2);

	dx = nodeJPtr->getCrds() - nodeIPtr->getCrds();

	// calculate the element length
	L = dx.Norm();

	if (L == 0.0)
	{
		opserr << "\nRockingBC::compElemtLengthAndOrien: 0 length\n";
		return -2;
	}

	// calculate the element local x axis components (direction cosines)
	// wrt to the global coordinates 
	cosTheta = dx(0) / L;
	sinTheta = dx(1) / L;

	return 0;
}

int
RockingBC::commitState()
{

  int retVal = 0;
  // call element commitState to do any base class stuff
  if ((retVal = this->Element::commitState()) != 0) {
    opserr << "RockingBC::commitState () - failed in base class";
  }    
  //retVal += theCoordTransf->commitState();
  
  uecommit = ue;
  kecommit = ke;
  Fecommit = Fe;
  sLcommit = sL;
  DWcommit = W - Wcommit;
  Dtcommit = Dt;

  curtime = this->DomainComponent::getDomain()->getCurrentTime();
  committedtime = curtime;

  Wcommit = W;
  Fn_com = Fn;
  FnVec_com = FnVec;
  FnD_com = FnD;
  Uel_com = Uel;

  dW_due_com = dW_due;

  if (useUelNM) {
	  Ysi_com = Ysi;
	  Si_com = Si;
	  Yupi_com = Yupi;
	  Upi_com = Upi;
  }
  else {
	  Ys_com = Ys;
	  S_com = S;
	  Ks_com = Ks;
	  Yup_com = Yup;
	  Up_com = Up;
	  Kup_com = Kup;
  }

  UB = UBnew;
  dUB_dR = dUBnew_dR;
  UB_R = UBnew_R;

  if (slidmode_com != newslidmode) {
	  std::cout << "Changed sliding mode into " << newslidmode << std::endl;
  }
  slidmode_com = newslidmode;
  
  comcount++;

  //Bilinearization
  if (useUelNM && blevery > 0) {
	  if (comcount%blevery == 0) {
		  for (size_t i = 0; i != Ysi_com.size(); i++) {

			  int_bilin(Ys_cats[i], Ysi_com[i], Si_com[i], Yupi_com[i], Upi_com[i], Ua_pos[i], ey, ysi_new, si_new, yupi_new, upi_new);
			  Ysi_com[i] = ysi_new;
			  Si_com[i] = si_new;
			  Yupi_com[i] = yupi_new;
			  Upi_com[i] = upi_new;
		  }
	  }
  }

  if (isdynamic) {
	  dyncount += 1;
	}
  triesfromcommitstate = 0;

  return retVal;
}

int
RockingBC::revertToLastCommit()
{
	//theCoordTransf->revertToLastCommit();
	ue = uecommit;
	ke = kecommit;
	Fe = Fecommit;
	sL = sLcommit;
	W = Wcommit;
	Fn = Fn_com;
	FnD = FnD_com;

	if (useUelNM) {
		Ysi = Ysi_com;
		Si = Si_com;
		Upi = Upi_com;
		Yupi = Yupi_com;
	}
	else {
		Up = Up_com;
		Yup = Yup_com;
		Kup = Kup_com;
		Ys = Ys_com;
		S = S_com;
		Ks = Ks_com;
		Uel = Uel_com;
	}

	slidmode = slidmode_com;

	curtime = committedtime;

	dW_due = dW_due_com;

	hasreverted = true;

	return 0;
}

int
RockingBC::revertToStart()
{
	//theCoordTransf->revertToStart();
	ue.Zero();
	this->getInitialStiff();
	return 0;
}

bool
RockingBC::is_analysis_dynamic(void)
{
	const Vector &velI = nodeIPtr->getTrialVel();
	const Vector &velJ = nodeJPtr->getTrialVel();
	const Vector &accelI = nodeIPtr->getTrialAccel();
	const Vector &accelJ = nodeJPtr->getTrialAccel();

	static Vector vg = Vector(6);
	static Vector ag = Vector(6);
	for (int i = 0; i < 3; i++) {
		vg(i) = velI(i);
		vg(i + 3) = velJ(i);
		ag(i) = accelI(i);
		ag(i + 3) = accelJ(i);
	}

	bool isdyn{};
	if (vg.Norm() + ag.Norm() == 0) {
		isdyn = false;
	}
	else {
		isdyn = true;
	}

	return isdyn;
}

int
RockingBC::update(void)
{

	triesfromcommitstate += 1;

	kepr = ke;
	Fepr = Fe;
	uepr = ue;
	sLpr = sL;
	Wpr = W;

	dW_due_pr = dW_due;

	//theCoordTransf->update();
	// get global displacements 
	const Vector& dispI = nodeIPtr->getTrialDisp();
	const Vector& dispJ = nodeJPtr->getTrialDisp();

	static Vector ug(6);
	for (int i = 0; i < 3; i++) {
		ug(i) = dispI(i);
		ug(i + 3) = dispJ(i);
	}

	// transform global end displacements to local coordinates

	ue(0) = cosTheta * ug(0) + sinTheta * ug(1);
	ue(1) = cosTheta * ug(1) - sinTheta * ug(0);
	ue(2) = ug(2);
	ue(3) = cosTheta * ug(3) + sinTheta * ug(4);
	ue(4) = cosTheta * ug(4) - sinTheta * ug(3);
	ue(5) = ug(5);

	// determine displacements in the basic system eliminating rigid body modes 
	//this->transfLocalDisplsToBasic(ul);

	// compute the transformation matrix from local to the basic system
	//this->compTransfMatrixBasicLocalDisps(Tbl_d);
	//this->compTransfMatrixBasicLocalForces(Tbl_f);

	if (hasreverted) {
		hasreverted = false;
		triesfromcommitstate = 0;
		return 0;
	} else {
		return this->state_determination();
	}

}

double
RockingBC::getDt(void) {
	curtime = this->DomainComponent::getDomain()->getCurrentTime();
	//std::cout << "Dt= " << curtime - committedtime << std::endl;
	return curtime - committedtime;
}

void
RockingBC::compTransfMatrixLocalGlobal(Matrix &Tlg)
{
	// setup transformation matrix from global to local coordinates
	Tlg.Zero();

	Tlg(0, 0) = Tlg(3, 3) = cosTheta;
	Tlg(0, 1) = Tlg(3, 4) = sinTheta;
	Tlg(1, 0) = Tlg(4, 3) = -sinTheta;
	Tlg(1, 1) = Tlg(4, 4) = cosTheta;
	Tlg(2, 2) = Tlg(5, 5) = 1.0;
}

const Vector &
RockingBC::getLocalTrialDisp(void)
{
	return ue;
}

const Vector &
RockingBC::getLocalIncrDeltaDisp(void)
{
	// dub = ub - ubpr;
	due = ue;
	due.addVector(1.0, uepr, -1.0);

	return due;
}

const Vector &
RockingBC::getLocalIncrDisp(void)
{
	// Dub = ub - ubcommit;
	Due = ue;
	Due.addVector(1.0, uecommit, -1.0);

	return Due;
}


double
RockingBC::getInitialLength(void)
{
	return L;
}

const Vector &
RockingBC::getGlobalResistingForce(const Vector &Fe)
{

	// transform resisting forces from the basic system to local coordinates
	//this->compTransfMatrixBasicLocalForces(Tbl_f);
	//static Vector pl(6);
	//pl.addMatrixTransposeVector(0.0, Tbl_f, pb, 1.0);    // pl = Tbl ^ pb;

	pg(0) = cosTheta*Fe[0] - sinTheta*Fe[1];
	pg(1) = sinTheta*Fe[0] + cosTheta*Fe[1];

	pg(3) = cosTheta*Fe[3] - sinTheta*Fe[4];
	pg(4) = sinTheta*Fe[3] + cosTheta*Fe[4];

	pg(2) = Fe[2];
	pg(5) = Fe[5];

	return pg;
}

const Matrix &
RockingBC::getGlobalStiffMatrix(const Matrix &ke)
{
	// transform tangent stiffness matrix from the basic system to local coordinates
	//static Matrix kl(6, 6);
	//this->compTransfMatrixBasicLocalForces(Tbl_f);
	//kl.addMatrixTripleProduct(0.0, Tbl_f, kb, 1.0);      // kl = Tbl ^ kb * Tbl;

													   // add geometric stiffness matrix
	//kl.addMatrix(1.0, this->getGeomStiffMatrix(pb), 1.0);

	// transform tangent  stiffness matrix from local to global coordinates

	// kg.addMatrixTripleProduct(0.0, Tlg, kl, 1.0);
	double s2, c2, cs;

	s2 = sinTheta*sinTheta;
	c2 = cosTheta*cosTheta;
	cs = sinTheta*cosTheta;

	double k11, k12, k13, k21, k22, k23, k31, k32, k33;

	k11 = ke(0, 0);    k12 = ke(0, 1);    k13 = ke(0, 2);
	k21 = ke(1, 0);    k22 = ke(1, 1);    k23 = ke(1, 2);
	k31 = ke(2, 0);    k32 = ke(2, 1);    k33 = ke(2, 2);

	kg(0, 0) = c2*k11 + s2*k22 - cs*(k21 + k12);
	kg(1, 0) = c2*k21 - s2*k12 + cs*(k11 - k22);
	kg(2, 0) = cosTheta*k31 - sinTheta*k32;

	kg(0, 1) = c2*k12 - s2*k21 + cs*(k11 - k22);
	kg(1, 1) = c2*k22 + s2*k11 + cs*(k21 + k12);
	kg(2, 1) = sinTheta*k31 + cosTheta*k32;

	kg(0, 2) = cosTheta*k13 - sinTheta*k23;
	kg(1, 2) = sinTheta*k13 + cosTheta*k23;
	kg(2, 2) = k33;

	k11 = ke(0, 3);    k12 = ke(0, 4);    k13 = ke(0, 5);
	k21 = ke(1, 3);    k22 = ke(1, 4);    k23 = ke(1, 5);
	k31 = ke(2, 3);    k32 = ke(2, 4);    k33 = ke(2, 5);

	kg(0, 3) = c2*k11 + s2*k22 - cs*(k21 + k12);
	kg(1, 3) = c2*k21 - s2*k12 + cs*(k11 - k22);
	kg(2, 3) = cosTheta*k31 - sinTheta*k32;

	kg(0, 4) = c2*k12 - s2*k21 + cs*(k11 - k22);
	kg(1, 4) = c2*k22 + s2*k11 + cs*(k21 + k12);
	kg(2, 4) = sinTheta*k31 + cosTheta*k32;

	kg(0, 5) = cosTheta*k13 - sinTheta*k23;
	kg(1, 5) = sinTheta*k13 + cosTheta*k23;
	kg(2, 5) = k33;

	k11 = ke(3, 0);    k12 = ke(3, 1);    k13 = ke(3, 2);
	k21 = ke(4, 0);    k22 = ke(4, 1);    k23 = ke(4, 2);
	k31 = ke(5, 0);    k32 = ke(5, 1);    k33 = ke(5, 2);

	kg(3, 0) = c2*k11 + s2*k22 - cs*(k21 + k12);
	kg(4, 0) = c2*k21 - s2*k12 + cs*(k11 - k22);
	kg(5, 0) = cosTheta*k31 - sinTheta*k32;

	kg(3, 1) = c2*k12 - s2*k21 + cs*(k11 - k22);
	kg(4, 1) = c2*k22 + s2*k11 + cs*(k21 + k12);
	kg(5, 1) = sinTheta*k31 + cosTheta*k32;

	kg(3, 2) = cosTheta*k13 - sinTheta*k23;
	kg(4, 2) = sinTheta*k13 + cosTheta*k23;
	kg(5, 2) = k33;

	k11 = ke(3, 3);    k12 = ke(3, 4);    k13 = ke(3, 5);
	k21 = ke(4, 3);    k22 = ke(4, 4);    k23 = ke(4, 5);
	k31 = ke(5, 3);    k32 = ke(5, 4);    k33 = ke(5, 5);

	kg(3, 3) = c2*k11 + s2*k22 - cs*(k21 + k12);
	kg(4, 3) = c2*k21 - s2*k12 + cs*(k11 - k22);
	kg(5, 3) = cosTheta*k31 - sinTheta*k32;

	kg(3, 4) = c2*k12 - s2*k21 + cs*(k11 - k22);
	kg(4, 4) = c2*k22 + s2*k11 + cs*(k21 + k12);
	kg(5, 4) = sinTheta*k31 + cosTheta*k32;

	kg(3, 5) = cosTheta*k13 - sinTheta*k23;
	kg(4, 5) = sinTheta*k13 + cosTheta*k23;
	kg(5, 5) = k33;

	return kg;
}

const Matrix &
RockingBC::getInitialGlobalStiffMatrix(const Matrix &kb)
{
	// transform tangent stiffness matrix from the basic system to local coordinates
	static Matrix kl(6, 6);
	static Matrix T(3, 6);

	T(0, 0) = -1.0;
	T(1, 0) = 0;
	T(2, 0) = 0;

	T(0, 1) = 0;
	T(1, 1) = 1 / L;
	T(2, 1) = 1 / L;

	T(0, 2) = 0;
	T(1, 2) = 0;
	T(2, 2) = 1;

	T(0, 3) = 1;
	T(1, 3) = 0;
	T(2, 3) = 0;

	T(0, 4) = 0;
	T(1, 4) = -1 / L;
	T(2, 4) = -1 / L;

	T(0, 5) = 0;
	T(1, 5) = 1;
	T(2, 5) = 0;

	kl.addMatrixTripleProduct(0.0, T, kb, 1.0);      // kl = Tbl ^ kb * Tbl;

													 // add geometric stiffness matrix
													 // kl.addMatrix(1.0, this->getGeomStiffMatrix(pb), 1.0);

													 // transform tangent  stiffness matrix from local to global coordinates

													 // kg.addMatrixTripleProduct(0.0, Tlg, kl, 1.0);
	double s2, c2, cs;

	s2 = sinTheta*sinTheta;
	c2 = cosTheta*cosTheta;
	cs = sinTheta*cosTheta;

	double k11, k12, k13, k21, k22, k23, k31, k32, k33;

	k11 = kl(0, 0);    k12 = kl(0, 1);    k13 = kl(0, 2);
	k21 = kl(1, 0);    k22 = kl(1, 1);    k23 = kl(1, 2);
	k31 = kl(2, 0);    k32 = kl(2, 1);    k33 = kl(2, 2);

	kg(0, 0) = c2*k11 + s2*k22 - cs*(k21 + k12);
	kg(1, 0) = c2*k21 - s2*k12 + cs*(k11 - k22);
	kg(2, 0) = cosTheta*k31 - sinTheta*k32;

	kg(0, 1) = c2*k12 - s2*k21 + cs*(k11 - k22);
	kg(1, 1) = c2*k22 + s2*k11 + cs*(k21 + k12);
	kg(2, 1) = sinTheta*k31 + cosTheta*k32;

	kg(0, 2) = cosTheta*k13 - sinTheta*k23;
	kg(1, 2) = sinTheta*k13 + cosTheta*k23;
	kg(2, 2) = k33;

	k11 = kl(0, 3);    k12 = kl(0, 4);    k13 = kl(0, 5);
	k21 = kl(1, 3);    k22 = kl(1, 4);    k23 = kl(1, 5);
	k31 = kl(2, 3);    k32 = kl(2, 4);    k33 = kl(2, 5);

	kg(0, 3) = c2*k11 + s2*k22 - cs*(k21 + k12);
	kg(1, 3) = c2*k21 - s2*k12 + cs*(k11 - k22);
	kg(2, 3) = cosTheta*k31 - sinTheta*k32;

	kg(0, 4) = c2*k12 - s2*k21 + cs*(k11 - k22);
	kg(1, 4) = c2*k22 + s2*k11 + cs*(k21 + k12);
	kg(2, 4) = sinTheta*k31 + cosTheta*k32;

	kg(0, 5) = cosTheta*k13 - sinTheta*k23;
	kg(1, 5) = sinTheta*k13 + cosTheta*k23;
	kg(2, 5) = k33;

	k11 = kl(3, 0);    k12 = kl(3, 1);    k13 = kl(3, 2);
	k21 = kl(4, 0);    k22 = kl(4, 1);    k23 = kl(4, 2);
	k31 = kl(5, 0);    k32 = kl(5, 1);    k33 = kl(5, 2);

	kg(3, 0) = c2*k11 + s2*k22 - cs*(k21 + k12);
	kg(4, 0) = c2*k21 - s2*k12 + cs*(k11 - k22);
	kg(5, 0) = cosTheta*k31 - sinTheta*k32;

	kg(3, 1) = c2*k12 - s2*k21 + cs*(k11 - k22);
	kg(4, 1) = c2*k22 + s2*k11 + cs*(k21 + k12);
	kg(5, 1) = sinTheta*k31 + cosTheta*k32;

	kg(3, 2) = cosTheta*k13 - sinTheta*k23;
	kg(4, 2) = sinTheta*k13 + cosTheta*k23;
	kg(5, 2) = k33;

	k11 = kl(3, 3);    k12 = kl(3, 4);    k13 = kl(3, 5);
	k21 = kl(4, 3);    k22 = kl(4, 4);    k23 = kl(4, 5);
	k31 = kl(5, 3);    k32 = kl(5, 4);    k33 = kl(5, 5);

	kg(3, 3) = c2*k11 + s2*k22 - cs*(k21 + k12);
	kg(4, 3) = c2*k21 - s2*k12 + cs*(k11 - k22);
	kg(5, 3) = cosTheta*k31 - sinTheta*k32;

	kg(3, 4) = c2*k12 - s2*k21 + cs*(k11 - k22);
	kg(4, 4) = c2*k22 + s2*k11 + cs*(k21 + k12);
	kg(5, 4) = sinTheta*k31 + cosTheta*k32;

	kg(3, 5) = cosTheta*k13 - sinTheta*k23;
	kg(4, 5) = sinTheta*k13 + cosTheta*k23;
	kg(5, 5) = k33;

	return kg;
}

const Matrix &
RockingBC::getGlobalMatrixFromLocal(const Matrix &ml)
{
	static Matrix kg(6, 6);
	kg.Zero();
	this->compTransfMatrixLocalGlobal(Tlg);
	kg.addMatrixTripleProduct(0.0, Tlg, ml, 1.0);

	return kg;
}

int
RockingBC::state_determination(void)
{

	const Vector& ue = this->getLocalTrialDisp();
	const Vector& due = this->getLocalIncrDeltaDisp();
	const Vector& Due = this->getLocalIncrDisp();

	ueV = ue;
	dueV = due;
	DueV = Due;
	DW.Zero();

	if (usecomstiff == -1) {
		DW.Zero();
		W = Wcommit;
	}
	else if (usecomstiff) {
		DW = af * dW_due_com * DueV;
		W = Wcommit + DW;
	}
	else {
		DW = af * dW_due * dueV;
		W += DW;
	}

	Dt = getDt();
	isdynamic = is_analysis_dynamic();
	if (isdynamic && Dt > 0) {
		beta_Dt = betaK / Dt;
	}
	else if (isdynamic && Dt == 0) {
		beta_Dt = -1; // Sto commitState
	}
	else if (isdynamic && Dt <= 0 && dyncount >= 1) {
		std::cout << "Error in Dt in dynamic analysis, Dt= " << Dt << std::endl;
		return -1;
	}
	else {
		beta_Dt = 0;
	}

	int NLsolvesuccess = 0;

	Winit = W;
	slidmode = slidmode_com;
	NLsolvesuccess = NL_solve_dyn();
	if ((NLsolvesuccess == 0 && slidmode != newslidmode) || NLsolvesuccess != 0) {
		slidingmodes_try.clear();

		bool ex2 = false;
		for (int j = 0; j < slidingmodes.size(); j++) {
			if (slidingmodes[j] == 2) {
				ex2 = true;
			}
		}

		if (slidmode == 0) {
			if (!ex2) {
				slidingmodes_try.push_back(1);
				slidingmodes_try.push_back(2);
			}
			else {
				slidingmodes_try.push_back(2);
				slidingmodes_try.push_back(1);
			}
		}
		else if (slidmode == 1) {
			slidingmodes_try.push_back(0);
			slidingmodes_try.push_back(2);
		}
		else {
			slidingmodes_try.push_back(0);
			slidingmodes_try.push_back(1);
		}

		for (int j = 0; j < slidingmodes_try.size(); j++) {
			slidmode = slidingmodes_try[j];
			//std::cout << "Trying sliding mode " << slidmode << std::endl;
			NLsolvesuccess = NL_solve_dyn();
			if (NLsolvesuccess == 0 && slidmode == newslidmode) {
				slidmode_init = slidmode;
				break;
			}
		}
	}

	for (int i = 0; i != 6; i++) {
		Fe(i) = FeV(i);
		for (int j = 0; j != 6; j++) {
			ke(i, j) = DFe_Due(i, j);
		}
	}

	//if (NLsolvesuccess != 0) {
	//	writedbgfile();
	//}

	if (isdynamic && Fst > 0 && curtime > committedtime) {
		forceratioN = std::fabs(FnVec[0] - FnVec_com[0]) / Fst;
		forceratioT = std::fabs(Fe[0] - Fecommit[0]) / Fst;

		if (triesfromcommitstate == 1) {
			forceratioNmax = 0;
			forceratioTmax = 0;
		}
		if (forceratioN > forceratioNmax && Due.Norm() > 0) {
			forceratioNmax = forceratioN;
		}
		if (forceratioT > forceratioTmax && Due.Norm() > 0) {
			forceratioTmax = forceratioT;
		}

		Dtprev = Dt;
		if (Dtprev <= 0) {
			std::cout << "Error in Dtprev in dynamic analysis, Dtprev = "<< Dtprev << std::endl;
			return -1;
		}
	}
	if (!isdynamic) {
		Fst = std::fabs(Fe[3]);
		//std::cout << Fst << std::endl;
	}

	//std::cout << "Dynamic step, Forces: " << FnVec[0] << " " << FnVec_com[0] << " " << Fe[0] << " " << Fecommit[0] << std::endl;
	//std::cout << "Dynamic step: " << dyncount << " forceratios: " << forceratioN << " " << forceratioT << " " << dueV.norm() << " " << DueV.norm() << " " << Dt << std::endl;

	if (errorifNexceeds && Due.Norm() > 0 && Dt > 1.01 * Dtlim && ( forceratioN > NlimN || forceratioT > NlimT ) ) {
		//std::cout << "Large rate of change of axial force: forceratioN= " << forceratioN << ", forceratioT= " << forceratioT <<", Dtcur=" << Dt << std::endl;
		return -1; // Raise Error
	}
	
	return NLsolvesuccess;
}

const Matrix &
RockingBC::getTangentStiff(void)
{
  return this->getGlobalStiffMatrix(ke);
}

const Matrix &
RockingBC::getInitialStiff(void)
{
  double L = this->getInitialLength();

  static Matrix fb(3, 3);
  static Matrix kb(3, 3);

  fb.Zero();
  fb(0, 0) = L / EA;
  fb(1, 1) = fb(2, 2) = L / 3 / EI + alpha / GA / L;
  fb(1, 2) = fb(2, 1) = -L / 6 / EI + alpha / GA / L;
  kb = inverse3x3matrix(fb);
    
  return this->getInitialGlobalStiffMatrix(kb);
}

const Matrix &
RockingBC::getDamp(void)
{
	K.Zero();

	return K;
}

const Matrix &
RockingBC::getMass(void)
{ 
    K.Zero();
    
    return K;
}

void 
RockingBC::zeroLoad(void)
{
  return;
}

int 
RockingBC::addLoad(ElementalLoad *theLoad, double loadFactor)
{
	opserr << "RockingBC::addLoad()  -- load type unknown for element with tag: " << this->getTag() << endln;
	return 1;
}

int
RockingBC::addInertiaLoadToUnbalance(const Vector &accel)
{
    return 0;
}

const Vector &
RockingBC::getResistingForceIncInertia()
{	
  P = this->getResistingForce();
  
  //// add the damping forces if rayleigh damping
  //if (alphaM != 0.0 || betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0)
  //  P.addVector(1.0, this->getRayleighDampingForces(), 1.0);
  
  return P;
}


const Vector &
RockingBC::getResistingForce()
{
  
  P = this->getGlobalResistingForce(Fe);
  
  return P;
}

int
RockingBC::sendSelf(int cTag, Channel &theChannel)
{
	return 1;
}

int
RockingBC::recvSelf(int cTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
	return 1;
}

void
RockingBC::Print(OPS_Stream &s, int flag)
{
  // to update forces!
  this->getResistingForce();

  if (flag == -1) {
    int eleTag = this->getTag();
    s << "RockingBC\t" << eleTag << "\t";
    s << 0 << "\t" << 0 << "\t" << connectedExternalNodes(0) << "\t" << connectedExternalNodes(1) ;
    s << "0\t0.0000000\n";
  } else {
    this->getResistingForce();
    s << "\nRockingBC: " << this->getTag() << endln;
    s << "\tConnected Nodes: " << connectedExternalNodes ;
    //s << "\tCoordTransf: " << this->getTag() << endln;
    double P  = Fe(3);
    double M1 = Fe(5);
    double M2 = Fe(2);
    double L = this->getInitialLength();
    double V = (M1+M2)/L;
    s << "\tEnd 1 Forces (P V M): " << -P
      << " " << V << " " << M2 << endln;
    s << "\tEnd 2 Forces (P V M): " << P
      << " " << -V << " " << M1 << endln;
  }
}

Response*
RockingBC::setResponse(const char **argv, int argc, OPS_Stream &output)
{

	Response *theResponse = 0;

	output.tag("ElementOutput");
	output.attr("eleType", "RockingBC");
	output.attr("eleTag", this->getTag());
	output.attr("node1", connectedExternalNodes[0]);
	output.attr("node2", connectedExternalNodes[1]);

	// global forces
	if (strcmp(argv[0], "force") == 0 || strcmp(argv[0], "forces") == 0 ||
		strcmp(argv[0], "globalForce") == 0 || strcmp(argv[0], "globalForces") == 0) {

		output.tag("ResponseType", "Px_1");
		output.tag("ResponseType", "Py_1");
		output.tag("ResponseType", "Mz_1");
		output.tag("ResponseType", "Px_2");
		output.tag("ResponseType", "Py_2");
		output.tag("ResponseType", "Mz_2");

		theResponse = new ElementResponse(this, 2, P);

		// local forces
	}
	else if (strcmp(argv[0], "localForce") == 0 || strcmp(argv[0], "localForces") == 0) {

		output.tag("ResponseType", "N_1");
		output.tag("ResponseType", "V_1");
		output.tag("ResponseType", "M_1");
		output.tag("ResponseType", "N_2");
		output.tag("ResponseType", "V_2");
		output.tag("ResponseType", "M_2");

		theResponse = new ElementResponse(this, 3, P);

		// basic forces
	}
	else if (strcmp(argv[0], "basicForce") == 0 || strcmp(argv[0], "basicForces") == 0) {

		output.tag("ResponseType", "N");
		output.tag("ResponseType", "M_1");
		output.tag("ResponseType", "M_2");

		theResponse = new ElementResponse(this, 4, Vector(3));

	}

	//	// deformations
	//}
	//else if (strcmp(argv[0], "deformatons") == 0 ||
	//	strcmp(argv[0], "basicDeformations") == 0) {

	//	output.tag("ResponseType", "eps");
	//	output.tag("ResponseType", "theta1");
	//	output.tag("ResponseType", "theta2");
	//	theResponse = new ElementResponse(this, 5, Vector(3));

	//	// chord rotation -
	//}
	//else if (strcmp(argv[0], "chordRotation") == 0 || strcmp(argv[0], "chordDeformation") == 0
	//	|| strcmp(argv[0], "basicDeformation") == 0) {

	//	output.tag("ResponseType", "eps");
	//	output.tag("ResponseType", "theta1");
	//	output.tag("ResponseType", "theta2");

	//	theResponse = new ElementResponse(this, 5, Vector(3));
	//}
	else if (strcmp(argv[0], "localDisplacements") == 0) {

		output.tag("ResponseType", "u1");
		output.tag("ResponseType", "v1");
		output.tag("ResponseType", "theta1");
		output.tag("ResponseType", "u2");
		output.tag("ResponseType", "v2");
		output.tag("ResponseType", "theta2");
		theResponse = new ElementResponse(this, 5, Vector(6));
	}

	else if (strcmp(argv[0], "sL") == 0 || strcmp(argv[0], "slip") == 0) {

		output.tag("ResponseType", "sL_com");

		theResponse = new ElementResponse(this, 6, Vector(1));

	}

	else if (strcmp(argv[0], "forceratioN") == 0) {

		output.tag("ResponseType", "forceratioN");

		theResponse = new ElementResponse(this, 7, Vector(1));

	}

	else if (strcmp(argv[0], "forceratioT") == 0) {

		output.tag("ResponseType", "forceratioT");

		theResponse = new ElementResponse(this, 8, Vector(1));

	}

	else if (strcmp(argv[0], "Dtmax") == 0) {

		output.tag("ResponseType", "Dtmax");

		theResponse = new ElementResponse(this, 9, Vector(1));

		// basic forces
	}

	else if (strcmp(argv[0], "forceratioNmax") == 0) {

		output.tag("ResponseType", "forceratioNmax");

		theResponse = new ElementResponse(this, 10, Vector(1));

	}

	else if (strcmp(argv[0], "forceratioTmax") == 0) {

		output.tag("ResponseType", "forceratioTmax");

		theResponse = new ElementResponse(this, 11, Vector(1));

	}

	else{
		std::string fstr = argv[0];
		Yup_file = std::ofstream(fstr + "_Yup.txt");
		Up_file = std::ofstream(fstr + "_Up.txt");
		Ys_file = std::ofstream(fstr + "_Ys.txt");
		S_file = std::ofstream(fstr + "_S.txt");

		theResponse = new ElementResponse(this, 20, Vector(1));

	}

	output.endTag(); // ElementOutput

	return theResponse;
}

int
RockingBC::getResponse(int responseID, Information &eleInfo)
{
	Vector vec1(1);
	double L = this->getInitialLength();
	this->getResistingForce();
	
	Vector vtemp;

	switch (responseID) {
	case 1: // stiffness
		return eleInfo.setMatrix(this->getTangentStiff());

	case 2: // global forces
		return eleInfo.setVector(this->getGlobalResistingForce(Fecommit));

	case 3: // local forces
		return eleInfo.setVector(Fecommit);

	case 4: // basic forces
		return eleInfo.setVector(FnVec_com);

	case 5:
		return eleInfo.setVector(this->getLocalTrialDisp());

	case 6:
		vec1(0) = sLcommit*L;
		return eleInfo.setVector(vec1);

	case 7:
		vec1(0) = forceratioN;
		return eleInfo.setVector(vec1);

	case 8:
		vec1(0) = forceratioT;
		return eleInfo.setVector(vec1);

	case 9:
		if (NlimN == 0 || forceratioN < 1.0e-12) {
			DtmaxN = -1;
		}
		else {
			DtmaxN = (NlimN / forceratioN) * (Dtprev / Dtlim);
		}

		if (NlimT == 0 || forceratioT < 1.0e-12) {
			DtmaxT = -1;
		}
		else {
			DtmaxT = (NlimT / forceratioT) * (Dtprev / Dtlim);
		}

		if (DtmaxN < 0 && DtmaxT < 0) {
			Dtmax = -1;
		}
		else if (DtmaxN < 0 && DtmaxT >= 0) {
			Dtmax = DtmaxT;
		}
		else if (DtmaxN >= 0 && DtmaxT < 0) {
			Dtmax = DtmaxN;
		}
		else {
			Dtmax = std::fmin(DtmaxN, DtmaxT);
		}

		if (Dtmax < 0.00001) {
			Dtmax = 0;
		} else if (Dtmax > 1000.0) {
			Dtmax = -1;
		}

		//std::cout << "Dtmax inside element " << DtmaxN << " " << DtmaxT << std::endl;
		vec1(0) = Dtmax;
		return eleInfo.setVector(vec1);

	case 10:
		vec1(0) = forceratioNmax;
		return eleInfo.setVector(vec1);

	case 11:
		vec1(0) = forceratioTmax;
		return eleInfo.setVector(vec1);

	case 20:
		if (useUelNM) {
			Ys_com = interval_join(Ysi_com);
			S_com = interval_join(Si_com);
			Yup_com = interval_join(Yupi_com);
			Up_com = interval_join(Upi_com);
		}

		for (size_t i = 0; i != Yup_com.Size(); i++) {
			Yup_file << Yup_com(i) * b << " ";
		}
		Yup_file << std::endl;

		for (size_t i = 0; i != Up_com.Size(); i++) {
			Up_file << Up_com(i) * b << " ";
		}
		Up_file << std::endl;

		for (size_t i = 0; i != Ys_com.Size(); i++) {
			Ys_file << Ys_com(i) * b << " ";
		}
		Ys_file << std::endl;

		for (size_t i = 0; i != S_com.Size(); i++) {
			S_file << S_com(i) << " ";
		}
		S_file << std::endl;

		return eleInfo.setVector(0);

	default:
		return -1;
	}
}

const Matrix &
RockingBC::inverse3x3matrix(Matrix &A) const
{
	double determinant = +A(0, 0)*(A(1, 1)*A(2, 2) - A(2, 1)*A(1, 2))
		- A(0, 1)*(A(1, 0)*A(2, 2) - A(1, 2)*A(2, 0))
		+ A(0, 2)*(A(1, 0)*A(2, 1) - A(1, 1)*A(2, 0));
	if (std::fabs(determinant) < 1.0e-100) {
		opserr << "Rocking BC determinant close to zero\n";
	}
	double invdet = 1.0 / determinant;
	static Matrix result = Matrix(3, 3);
	result(0, 0) = (A(1, 1)*A(2, 2) - A(2, 1)*A(1, 2))*invdet;
	result(0, 1) = -(A(0, 1)*A(2, 2) - A(0, 2)*A(2, 1))*invdet;
	result(0, 2) = (A(0, 1)*A(1, 2) - A(0, 2)*A(1, 1))*invdet;
	result(1, 0) = -(A(1, 0)*A(2, 2) - A(1, 2)*A(2, 0))*invdet;
	result(1, 1) = (A(0, 0)*A(2, 2) - A(0, 2)*A(2, 0))*invdet;
	result(1, 2) = -(A(0, 0)*A(1, 2) - A(1, 0)*A(0, 2))*invdet;
	result(2, 0) = (A(1, 0)*A(2, 1) - A(2, 0)*A(1, 1))*invdet;
	result(2, 1) = -(A(0, 0)*A(2, 1) - A(2, 0)*A(0, 1))*invdet;
	result(2, 2) = (A(0, 0)*A(1, 1) - A(1, 0)*A(0, 1))*invdet;

	return result;

}


// RCSL functions

Vector RockingBC::find_in_dist(const Vector& X, const Vector& Y, const Vector& Xf) {
	static std::vector<double> Yf{}; Yf.clear();
	int ix = 0;
	for (size_t i = 0; i != Xf.Size(); i++) {
		while (Xf[i] != X[ix]) {
			ix += 1;
		}
		Yf.push_back(Y[ix]);
	}

	static Vector YfXd;
	YfXd = Vector(Yf.size());
	for (size_t i = 0; i != Yf.size(); i++) {
		YfXd[i] = Yf[i];
	}
	return YfXd;

}

void RockingBC::simplify_dist_up(const Vector& X, const Vector& Y, const Vector& XwXd, Vector& XnewXd, Vector& YnewXd)
{
	static std::vector<double> Xnew{}; Xnew.clear();
	static std::vector<double> Ynew{}; Ynew.clear();
	
	//(XwXd.data(), XwXd.data() + XwXd.rows() * XwXd.cols());
	std::vector<double> Xw(XwXd.Size());
	for (size_t i = 0; i != XwXd.Size(); i++) {
		Xw[i] = XwXd(i);
	}

	Xnew.push_back(X[0]);
	Ynew.push_back(Y[0]);
	for (size_t i = 1; i!=X.Size() - 1; i++) {
		if (std::find(Xw.begin(), Xw.end(), X[i]) != Xw.end() && std::find(Xnew.begin(), Xnew.end(), X[i]) == Xnew.end()) {
			Xnew.push_back(X[i]);
			Ynew.push_back(Y[i]);
		}
		else {
			double A = X[i - 1] * (Y[i + 1] - Y[i]) + X[i] * (Y[i - 1] - Y[i + 1]) + X[i + 1] * (Y[i] - Y[i - 1]);
			if (std::fabs(A) < 1.0e-16) {
				continue;
			}
			Xnew.push_back(X[i]);
			Ynew.push_back(Y[i]);
		}
	}

	Xnew.push_back(X[X.Size() - 1]);
	Ynew.push_back(Y[Y.Size() - 1]);
	XnewXd = Vector(Xnew.size());
	YnewXd = Vector(Ynew.size());
	for (size_t i = 0; i != Xnew.size(); i++) {
		XnewXd[i] = Xnew[i];
		YnewXd[i] = Ynew[i];
	}

	return;
}

void RockingBC::W_to_ua_upl()
{
	double DAMPC = 1.0;
	if (beta_Dt >= 0) {
		DAMPC = beta_Dt / (1.0 + beta_Dt);
	}

	for (size_t i = 0; i != Si_com.size(); i++) {
		Ec[i] = Si_com[i][0];
		El[i] = Si_com[i][0] * DAMPC;
		Upl[i] = Upi_com[i][0];
	}
	Ec[Nw-1] = Si_com[Nw-2][(Si_com[Nw - 2]).size()-1];
	El[Nw - 1] = Si_com[Nw - 2][(Si_com[Nw - 2]).size() - 1]*DAMPC;
	Upl[Nw-1] = Upi_com[Nw-2][(Upi_com[Nw - 2]).size() - 1];

	dUa_dW.Zero();
	for (size_t i = 0; i != W.Size(); i++)
	{
		if (W(i) > El[i])
		{
			Ua(i) = W(i)-El[i];
			dUa_dW(i, i) = 1.0;
		}
		else if (W(i) <= ey)
		{
			Ua(i) = W(i) - ey;
			dUa_dW(i, i) = 1.0;
		}
		else
		{
			Ua(i) = 0.0;
		}
	}
	return;
}

void RockingBC::W_to_ua_upl_K()
{
	double DAMPC = 1.0;
	if (beta_Dt >= 0) {
		DAMPC = beta_Dt / (1.0 + beta_Dt);
	}

	Ec = find_in_dist(Ys_com, S_com, Yw);
	Upl = find_in_dist(Yup_com, Up_com, Yw);
	El = Ec * DAMPC;

	dUa_dW.Zero();
	for (size_t i = 0; i != W.Size(); i++)
	{
		if (W(i) > El[i])
		{
			Ua(i) = W(i) - El[i];
			dUa_dW(i, i) = 1.0;
		}
		else if (W(i) <= ey)
		{
			Ua(i) = W(i) - ey;
			dUa_dW(i, i) = 1.0;
		}
		else
		{
			Ua(i) = 0.0;
		}
	}
	return;
}


void RockingBC::Youter_calc()
{

	bool Szeros = true;
	for (size_t i = 0; i != Ys_cats_dist.size(); i++) {
		if (Ys_cats_dist[i] > 0) {
			Szeros = false;
			break;
		}
	}

	int zl = 0;
	int zr = Ys.Size() - 1;

	if (!Szeros)
	{
		while (true)
		{
			if (Ys_cats_dist[zl]>0)
			{
				break;
			}
			else
			{
				zl += 1;
			}
		}

		while (true)
		{
			if (Ys_cats_dist[zr - 1]>0)
			{
				break;
			}
			else
			{
				zr -= 1;
			}
		}
	}

	Youter(0) = Ys(zl);
	Youter(1) = Ys(zr);
	for (size_t i = 0; i != W.Size(); i++) {
		dYouter_dW(0, i) = dYs_dW(zl, i);
		dYouter_dW(1, i) = dYs_dW(zr, i);
	}

	return;
}

void RockingBC::NM_calc()
{
	N = 0.0;
	M = 0.0;
	dN_dW.Zero();
	dM_dW.Zero();

	for (size_t i = 0; i != Nints.Size(); i++)
	{
		N += Nints[i];
		M += Mints[i];
		for (size_t j = 0; j != W.Size(); j++) {
			dN_dW(j) += dNints_dW(i,j);
			dM_dW(j) += dMints_dW(i,j);
		}
	}

	return;

}

void RockingBC::NM_calc_YS()
{
	N = 0.0;
	M = 0.0;
	dN_dW.Zero();
	dM_dW.Zero();

	double y1{};
	double y2{};
	double s1{};
	double s2{};

	for (size_t i = 0; i != Ys.Size() - 1; i++)
	{
		y1 = Ys(i);
		y2 = Ys(i + 1);
		s1 = S(i);
		s2 = S(i + 1);

		N += (y2 - y1) * (s1 + s2) / 2.;
		M += (y2 - y1) * (2 * s1 * y1 + s1 * y2 + s2 * y1 + 2 * s2 * y2) / 6.;

		for (size_t j = 0; j != W.Size(); j++) {
			dN_dW(j) += (-s1 / 2. - s2 / 2.) * dYs_dW(i,j) + (s1 / 2. + s2 / 2.) * dYs_dW(i + 1, j) + (y2 / 2. - y1 / 2.) * dS_dW(i, j) + (y2 / 2. - y1 / 2.) * dS_dW(i + 1, j);
			dM_dW(j) += (-(s1 * y1) / 3. - (s1 * y2) / 6. - (s2 * y1) / 6. - (s2 * y2) / 3. - ((2 * s1 + s2) * (y1 - y2)) / 6.) * dYs_dW(i, j) +
				((s1 * y1) / 3. + (s1 * y2) / 6. + (s2 * y1) / 6. + (s2 * y2) / 3. - ((s1 + 2 * s2) * (y1 - y2)) / 6.) * dYs_dW(i + 1, j) +
				(-((y1 - y2) * (2 * y1 + y2)) / 6.) * dS_dW(i, j) + (-((y1 - y2) * (y1 + 2 * y2)) / 6.) * dS_dW(i + 1, j);
		}

	}

	return;

}

void RockingBC::NM_calc_Fncom()
{
	
	double DAMPC = 1.0;
	if (beta_Dt >= 0) {
		DAMPC = beta_Dt / (1.0 + beta_Dt);
	}

	N = Fn_com[0] * DAMPC;
	M = Fn_com[1] * DAMPC;
	dN_dW.Zero();
	dM_dW.Zero();

	double y1{};
	double y2{};
	double s1{};
	double s2{};

	for (size_t i = 0; i != Ydks.Size() - 1; i++)
	{
		y1 = Ydks(i);
		y2 = Ydks(i + 1);
		s1 = DS(i);
		s2 = DS(i + 1);

		N += (y2 - y1) * (s1 + s2) / 2.;
		M += (y2 - y1) * (2 * s1 * y1 + s1 * y2 + s2 * y1 + 2 * s2 * y2) / 6.;

		for (size_t j = 0; j != W.Size(); j++) {
			dN_dW(j) += (-s1 / 2. - s2 / 2.) * dYdks_dW(i, j) + (s1 / 2. + s2 / 2.) * dYdks_dW(i + 1, j) + (y2 / 2. - y1 / 2.) * dDS_dW(i, j) + (y2 / 2. - y1 / 2.) * dDS_dW(i + 1, j);
			dM_dW(j) += (-(s1 * y1) / 3. - (s1 * y2) / 6. - (s2 * y1) / 6. - (s2 * y2) / 3. - ((2 * s1 + s2) * (y1 - y2)) / 6.) * dYdks_dW(i, j) +
				((s1 * y1) / 3. + (s1 * y2) / 6. + (s2 * y1) / 6. + (s2 * y2) / 3. - ((s1 + 2 * s2) * (y1 - y2)) / 6.) * dYdks_dW(i + 1, j) +
				(-((y1 - y2) * (2 * y1 + y2)) / 6.) * dDS_dW(i, j) + (-((y1 - y2) * (y1 + 2 * y2)) / 6.) * dDS_dW(i + 1, j);
		}
	}

	return;

}

void RockingBC::fr_calc()
{

	fr(0, 0) = L / 2.;
	fr(1, 1) = 3. / 4. * L / b;
	fr(1, 2) = L * L / 4. / b / b - alpha * (1. + nu);
	fr(2, 1) = -3. / 4. * L / b;
	fr(2, 2) = -L * L / 2. / b / b - alpha * (1. + nu);

	return;
}

void RockingBC::sL_Q_t_calc()
{
	N_com = Fn_com[0];
	Q_com = Fn_com[2];
	ND_com = FnD_com[0];
	QD_com = FnD_com[2];

	w1 = ueV(2);
	r = (ueV(4) - ueV(1)) / L;
	w2 = ueV(5);

	if (beta_Dt >= 0) {
		Ntot = (1.0 + beta_Dt) * N - beta_Dt * N_com;
		dNtot_dW = (1.0 + beta_Dt) * dN_dW;
		if (Ntot > 0) { Ntot = 0; dNtot_dW = 0 * dNtot_dW; }

		PA = (1. + beta_Dt) * (k1 * (w1 - r) - k2 * M) - beta_Dt * Q_com + r * Ntot;
		PB = Ntot - k1 * (1. + beta_Dt);
		dPA_dW = (1. + beta_Dt) * (-k2 * dM_dW) + r * dNtot_dW;
		dPB_dW = dNtot_dW;
		dPA_due = (1. + beta_Dt) * k1 * (dw1_due - dr_due) + dr_due * Ntot;
	}
	else {
		Ntot = N + ND_com;
		dNtot_dW = dN_dW;
		if (Ntot > 0) { Ntot = 0; dNtot_dW = 0 * dNtot_dW; }

		PA = k1 * (w1 - r) - k2 * M + QD_com + r * Ntot;
		PB = Ntot - k1;
		dPA_dW = (-k2 * dM_dW) + r * dNtot_dW;
		dPB_dW = dNtot_dW;
		dPA_due = k1 * (dw1_due - dr_due) + dr_due * Ntot;
	}

	bool notsliding{ false };

	if (mu==0) {
		notsliding = false;
		for (int i = 0; i != Ys_cats_dist.size(); i++) {
			if (Ys_cats_dist[i] > 0) {
				notsliding = true;
				break;
			}
		}
		lim1 = 0.;
		lim2 = 0.;
		dlim1_dW.Zero();
		dlim2_dW.Zero();
		dlim1_due.Zero();
		dlim2_due.Zero();
	}
	else {
		lim1 = Ntot*(w2 + mu) / (1 - mu*w2);
		lim2 = Ntot*(w2 - mu) / (1 + mu*w2);
		dlim1_dW = dNtot_dW*(w2 + mu) / (1 - mu*w2);
		dlim2_dW = dNtot_dW*(w2 - mu) / (1 + mu*w2);
		dlim1_due = Ntot*dw2_due / (1 - mu*w2) + Ntot*(w2 + mu) / (1 - mu*w2) / (1 - mu*w2) * mu*dw2_due;
		dlim2_due = Ntot*dw2_due / (1 + mu*w2) - Ntot*(w2 - mu) / (1 + mu*w2) / (1 + mu*w2) * mu*dw2_due;
	}

	cval = PA + PB*sLcommit;
	if (notsliding || (cval >= lim1 && cval <= lim2) || PB==0) {
		newslidmode = 0;
	}
	else if (cval <= lim1) {
		newslidmode = 1;
	}
	else {
		newslidmode = 2;
	}
	
	if ((usespecslidmode && slidmode==0) || (!usespecslidmode && newslidmode==0)) {
		sL = sLcommit;
		dsL_dW.Zero();
		dsL_due.Zero();
	}
	else if ((usespecslidmode && slidmode == 1) || (!usespecslidmode && newslidmode == 1)) {
		sL = (lim1 - PA) / PB;
		dsL_dW = (dlim1_dW - dPA_dW) / PB - (lim1 - PA) / PB / PB*dPB_dW;
		dsL_due = (dlim1_due-dPA_due) / PB;
	}
	else {
		sL = (lim2 - PA) / PB;
		dsL_dW = (dlim2_dW - dPA_dW) / PB - (lim2 - PA) / PB / PB*dPB_dW;
		dsL_due = (dlim2_due - dPA_due) / PB;
	}
	
	Q = k1*(w1 - r - sL) - k2*M;
	gQ = 2. / 3. * (Youter(1) - Youter(0));
	t = Q / gQ;
	dQ_dW = k1*(-1.0*dsL_dW) - k2*dM_dW;
	dQ_due = k1*(dw1_due - dr_due - dsL_due);
	for (size_t j = 0; j != W.Size(); j++) {
		dgQ_dW(j) = 2. / 3. * (dYouter_dW(1, j) - dYouter_dW(0, j));
	}
	dt_dW = dQ_dW / gQ - Q / gQ / gQ * dgQ_dW;
	dt_due = dQ_due / gQ;

	return;

}

void RockingBC::un_calc()
{
	ues = ueV;
	ues(3) -= sL*L*ueV(5);
	ues(4) += sL*L;

	dues_due.Zero();
	dues_due(0, 0) = dues_due(1, 1) = dues_due(2, 2) = dues_due(3, 3) = dues_due(4, 4) = dues_due(5, 5) = 1.0;
	dues_due(3, 5) -= sL*L;
	for (size_t j = 0; j != ue.Size(); j++) {
		dues_due(3,j) -= dsL_due(j) * L * ueV(5);
		dues_due(4,j) += dsL_due(j) * L;
	}

	dues_dW.Zero();
	for (size_t j = 0; j != W.Size(); j++) {
		dues_dW(3, j) -= dsL_dW(j) * L * ueV(5);
		dues_dW(4, j) += dsL_dW(j) * L;
	}

	Tn(0, 0) = -1.0;
	Tn(0, 1) = -0.5 * (ues(4) - ues(1)) / L;
	Tn(0, 3) = 1.0;
	Tn(0, 4) = +0.5 * (ues(4) - ues(1)) / L;
	Tn(1, 1) = +1. / L;
	Tn(1, 4) = -1. / L;
	Tn(1, 5) = 1.0;
	Tn(2, 1) = +1. / L;
	Tn(2, 2) = 1.0;
	Tn(2, 4) = -1. / L;

	un = Tn*ues;

	dun_dues(0, 0) = -1.0;
	dun_dues(0, 1) = -(ues(4) - ues(1)) / L;
	dun_dues(0, 3) = 1.0;
	dun_dues(0, 4) = +(ues(4) - ues(1)) / L;
	dun_dues(1, 1) = +1. / L;
	dun_dues(1, 4) = -1. / L;
	dun_dues(1, 5) = 1.0;
	dun_dues(2, 1) = +1. / L;
	dun_dues(2, 2) = 1.0;
	dun_dues(2, 4) = -1. / L;

	dun_due = dun_dues*dues_due;
	dun_dW = dun_dues*dues_dW;

	return;

}

void RockingBC::ut_calc()
{
	if (noshear) {
		ut.Zero();
		dut_dW.Zero();
		dut_due.Zero();
		return;
	}

	dutn_dW = dutn_dYouter*dYouter_dW;
	ut = t*utn;

	dut_dW = t*dutn_dW;
	for (int k = 0; k != dt_dW.Size(); k++)
	{
		for (int i = 0; i != utn.Size(); i++)
		{
			dut_dW(i, k) += dt_dW(k) * utn(i);
		}
	}

	dut_due.Zero();
	for (int k = 0; k != dt_due.Size(); k++)
	{
		for (int i = 0; i != utn.Size(); i++)
		{
			dut_due(i, k) = dt_due(k) * utn(i);
		}
	}

	return;

}

void RockingBC::urf_calc()
{
	th2 = un(2);
	for (size_t j = 0; j != W.Size(); j++) {
		dth2_dW(j) = dun_dW(2,j);
	}
	for (size_t j = 0; j != ue.Size(); j++) {
		dth2_due(j) = dun_due(2, j);
	}

	urth(0) = 0;
	urth(1) = fr(1, 2) / fr(2, 2) * th2;

	for (size_t j = 0; j != W.Size(); j++) {
		durth_dW(1, j) = fr(1, 2) / fr(2, 2) * dth2_dW(j);
	}
	for (size_t j = 0; j != ue.Size(); j++) {
		durth_due(1, j) = fr(1, 2) / fr(2, 2) * dth2_due(j);
	}

	Fn2(0) = N;
	Fn2(1) = M;
	for (size_t j = 0; j != W.Size(); j++) {
		dFn2_dW(0, j) = dN_dW(j);
		dFn2_dW(1, j) = dM_dW(j);
	}

	urf = frr*Fn2 + urth;
	durf_dW = frr*dFn2_dW + durth_dW;
	durf_due = durth_due;

	return;
}

void RockingBC::Uel_NM_calc()
{

	Uel = UN*Nints + UM*Mints;
	dUel_dW = UN*dNints_dW + UM*dMints_dW;

	return;
}

void RockingBC::Uel_K_calc()
{
	double DAMPC = 1.0;
	if (beta_Dt >= 0) {
		DAMPC = beta_Dt / (1.0 + beta_Dt);
	}

	rnotfound.clear();
	rfoundi.clear();
	ifound.clear();
	inotfound.clear();
	rnfi = 0;
	ifi = 0;

	for (int i = 0; i != Ydks.Size(); i++) {
		while (ifi < UB_R.size() && UB_R[ifi] < Ydks[i]) {
			ifi += 1;
		}
		if (ifi < UB_R.size() && UB_R[ifi] == Ydks[i]) {
			rfoundi.push_back(ifi);
			ifound.push_back(i);
		}
		else {
			rnotfound.push_back(Ydks[i]);
			inotfound.push_back(i);
			rnfi += 1;
		}
	}

	UBnew_R = Vec(Ydks.Size(), 0.0);
	UBnew = Matrix(Nw, Ydks.Size());
	dUBnew_dR = Matrix(Nw, Ydks.Size());

	if (rnotfound.size() > 0) {
		rnotfoundvec = Vector(rnotfound.size());
		for (size_t j = 0; j != rnotfound.size(); j++) {
			rnotfoundvec(j) = rnotfound[j];
		}
	}
	else {
		rnotfoundvec = Vector(0);
	}
	Unf = Matrix(Nw, rnotfound.size());
	dUnf_dR = Matrix(Nw, rnotfound.size());
	triangle_dispslope_disps_givenMat1(rnotfoundvec, Yw, Im1, Jm1, Unf, dUnf_dR);

	for (int i = 0; i != ifound.size(); i++) {
		UBnew_R[ifound[i]] = UB_R[rfoundi[i]];
		for (int k = 0; k != Nw; k++) {
			UBnew(k, ifound[i]) = UB(k, rfoundi[i]);
			dUBnew_dR(k, ifound[i]) = dUB_dR(k, rfoundi[i]);
		}
	}
	for (int i = 0; i != inotfound.size(); i++) {
		UBnew_R[inotfound[i]] = rnotfound[i];
		for (int k = 0; k != Nw; k++) {
			UBnew(k, inotfound[i]) = Unf(k, i);
			dUBnew_dR(k, inotfound[i]) = dUnf_dR(k, i);
		}
	}

	DDKs = Vector(Ydks.Size());
	dDDKs_dW = Matrix(Ydks.Size(), Nw);
	for (int i = 0; i != Dks.Size() - 1; i++) {
		DDKs[i + 1] = Dks[i + 1] - Dks[i];
		for (size_t j = 0; j != W.Size(); j++) {
			dDDKs_dW(i + 1,j) = dDks_dW(i + 1,j) - dDks_dW(i,j);
		}
	}

	Uel = Uel_com * DAMPC + UBnew * DDKs;
	dUel_dW = UBnew * dDDKs_dW;
	for (size_t i = 0; i != Nw; i++)
	{
		for (size_t k = 0; k != Nw; k++)
		{
			for (size_t j = 0; j != DDKs.Size(); j++)
			{
				dUel_dW(i, k) += dUBnew_dR(i, j) * dYdks_dW(j, k) * DDKs[j];
			}
		}
	}

return;
}

void RockingBC::disp_comb()
{
	utar(0) = un(0);
	utar(1) = un(1);
	for (size_t j = 0; j != W.Size(); j++) {
		dutar_dW(0,j) = dun_dW(0,j); dutar_dW(1,j) = dun_dW(1,j);
	}
	for (size_t j = 0; j != ue.Size(); j++) {
		dutar_due(0, j) = dun_due(0, j); dutar_due(1, j) = dun_due(1, j);
	}

	Ut = CC*ut;
	dUt_dW = CC*dut_dW;
	dUt_due = CC*dut_due;

	Urf = CB*urf;
	dUrf_dW = CB*durf_dW;
	dUrf_due = CB*durf_due;

	Utar = CB*utar;
	dUtar_dW = CB*dutar_dW;
	dUtar_due = CB*dutar_due;

	Ud = Uel + Ua + Upl + Ut + Urf - Utar;
	dUd_dW = dUel_dW + dUa_dW + dUt_dW + dUrf_dW - dUtar_dW;
	dUd_due = dUt_due + dUrf_due - dUtar_due;

	return;
}

void RockingBC::forces()
{
	Fn(0) = N; Fn(1) = M; Fn(2) = Q;
	for (size_t j = 0; j != W.Size(); j++) {
		dFn_dW(0, j) = dN_dW(j);
		dFn_dW(1, j) = dM_dW(j);
		dFn_dW(2, j) = dQ_dW(j);
	}
	for (size_t j = 0; j != ue.Size(); j++) {
		dFn_due(2, j) = dQ_due(j);
	}

	if (beta_Dt >= 0) {
		FnD = beta_Dt * (Fn - Fn_com);
		dFnD_dW = beta_Dt * dFn_dW;
		dFnD_due = beta_Dt * dFn_due;
	}
	else {
		FnD = FnD_com;
		dFnD_dW = 0 * dFn_dW;
		dFnD_due = 0 * dFn_due;
	}

	Fntot = Fn + FnD;
	dFntot_dW = dFn_dW + dFnD_dW;
	dFntot_due = dFn_due + dFnD_due;

	Fnntot = TF*Fntot;
	dFnntot_dW = TF*dFntot_dW;
	dFnntot_due = TF*dFntot_due;

	FnNN = TF * Fn;
	FnVec(0) = FnNN(0);
	FnVec(1) = FnNN(1);
	FnVec(2) = FnNN(2);

	Fes(0) = -Fnntot(0);
	Fes(1) = -(ues(4) - ues(1)) / L * Fnntot(0) + 1. / L * Fnntot(1) + 1. / L * Fnntot(2);
	Fes(2) = Fnntot(2);
	Fes(3) = Fnntot(0);
	Fes(4) = (ues(4) - ues(1)) / L * Fnntot(0) - 1. / L * Fnntot(1) - 1. / L * Fnntot(2);
	Fes(5) = Fnntot(1);

	for (size_t j = 0; j != W.Size(); j++) {
		dFes_dW(0,j) = -dFnntot_dW(0, j);
		dFes_dW(1, j) = -(ues(4) - ues(1)) / L * dFnntot_dW(0, j) + 1. / L * dFnntot_dW(1, j) + 1. / L * dFnntot_dW(2, j)
			- (dues_dW(4, j) - dues_dW(1, j)) / L * Fnntot(0);
		dFes_dW(2, j) = dFnntot_dW(2, j);
		dFes_dW(3, j) = dFnntot_dW(0, j);
		dFes_dW(4, j) = (ues(4) - ues(1)) / L * dFnntot_dW(0, j) - 1. / L * dFnntot_dW(1, j) - 1. / L * dFnntot_dW(2, j)
			+ (dues_dW(4, j) - dues_dW(1, j)) / L * Fnntot(0);
		dFes_dW(5, j) = dFnntot_dW(1, j);
	}

	for (size_t j = 0; j != ue.Size(); j++) {
		dFes_due(0,j) = -dFnntot_due(0, j);
		dFes_due(1, j) = -(ues(4) - ues(1)) / L * dFnntot_due(0, j) + 1. / L * dFnntot_due(1, j) + 1. / L * dFnntot_due(2, j)
			- (dues_due(4, j) - dues_due(1, j)) / L * Fnntot(0);
		dFes_due(2, j) = dFnntot_due(2, j);
		dFes_due(3, j) = dFnntot_due(0, j);
		dFes_due(4, j) = (ues(4) - ues(1)) / L * dFnntot_due(0, j) - 1. / L * dFnntot_due(1, j) - 1. / L * dFnntot_due(2, j)
			+ (dues_due(4, j) - dues_due(1, j)) / L * Fnntot(0);
		dFes_due(5, j) = dFnntot_due(1, j);
	}

	FeV = Fes;
	FeV(5) -= Fes(3) * sL*L;
	FeV(5) -= Fes(4) * sL*L*ueV(5);

	dFeV_dW = dFes_dW;
	for (size_t j = 0; j != W.Size(); j++) {
		dFeV_dW(5,j) -= (dFes_dW(3,j) * sL * L + Fes(3) * dsL_dW(j) * L);
		dFeV_dW(5,j) -= (dFes_dW(4,j) * sL * L * ueV(5) + Fes(4) * dsL_dW(j) * L * ueV(5));
	}
	dFeV_due = dFes_due;
	for (size_t j = 0; j != ue.Size(); j++) {
		dFeV_due(5,j) -= (dFes_due(3,j) * sL * L + Fes(3) * dsL_due(j) * L);
		dFeV_due(5,j) -= (dFes_due(4,j) * sL * L * ueV(5) + Fes(4) * dsL_due(j) * L * ueV(5) + Fes(4) * sL * L * due5_due(j));
	}

	//if (!RockingBC_dll) {
		int solvesuc = dUd_dW.Solve(-1.0 * dUd_due, dW_due);
	//}
	//else {
		
	//CHANGE FOR DLL CREATION
		
	//	Eigen::MatrixXd dUd_dW_2 = Eigen::MatrixXd(Ud.Size(), W.Size());
	//	for (size_t i = 0; i != Ud.Size(); i++) {
	//		for (size_t j = 0; j != W.Size(); j++) {
	//			dUd_dW_2(i, j) = dUd_dW(i, j);
	//		}
	//	}
	//	Eigen::MatrixXd dUd_due_2 = Eigen::MatrixXd(Ud.Size(), ue.Size());
	//	for (size_t i = 0; i != Ud.Size(); i++) {
	//		for (size_t j = 0; j != ue.Size(); j++) {
	//			dUd_due_2(i, j) = dUd_due(i, j);
	//		}
	//	}

	//	Eigen::MatrixXd dW_due_2 = -dUd_dW_2.colPivHouseholderQr().solve(dUd_due_2);

	//	for (size_t i = 0; i != W.Size(); i++) {
	//		for (size_t j = 0; j != ue.Size(); j++) {
	//			dW_due(i, j) = dW_due_2(i, j);
	//		}
	//	}
	//}

	
	DFe_Due = dFeV_due + dFeV_dW*dW_due;

	return;
}

void RockingBC::WZ_solve()
{

	if (useUelNM) {
		
		interval_dists(Yw, W, Yupi_com, Upi_com, Ysi_com, Si_com, ey, beta_Dt,
			Ysi, Si, Yupi, Upi, Ys_cats, Nints, Mints, dNints_dW, dMints_dW, Ua_pos, dYsi_dW, dSi_dW);
		
		Ys_cats_dist_calc(Ys_cats,Ys_cats_dist);
		W_to_ua_upl();
		if (!noshear) {
			Ys = interval_join(Ysi);
			dYs_dW = interval_join(dYsi_dW);
			Youter_calc();
		}
		NM_calc();
		sL_Q_t_calc();
		un_calc();
		
		if (!noshear) {
			se_shear_1der(Youter, utn, dutn_dYouter);
			ut_calc();
		}
		urf_calc();
		Uel_NM_calc();
	}
	else {
		interval_dists_K(Yw, W, Yup_com, Up_com, Kup_com, Ys_com, S_com, Ks_com, ey, beta_Dt,
			Ys, S, Ks, Yup, Up, Kup, dYs_dW, dS_dW, dKs_dW, Ys_cats_dist, Ydks, Dks, dYdks_dW, dDks_dW, DS, dDS_dW);
		W_to_ua_upl_K();
		if (!noshear) {
			Youter_calc();
		}
		NM_calc_Fncom();
		sL_Q_t_calc();
		un_calc();
		if (!noshear) {
			se_shear_1der(Youter, utn, dutn_dYouter);
			ut_calc();
		}
		urf_calc();
		Uel_K_calc();
	}
	disp_comb();
	forces();

	return;
}

int RockingBC::NL_solve_dyn()
{

	int tries = 0;
	int NLsolvesuccess = 0;
	double curconvlim = convlim;
	double curaf = af;
	slidingmodes.clear();

	W = Winit;

	int itry = 0;
	while (true) {

		//writedbgfile();
		WZ_solve();
		for (int i = 0; i < slidingmodes.size(); i++) {
			if (slidingmodes[i] != newslidmode) {
				slidingmodes.push_back(newslidmode);
			}
		}
		if (Ud.Norm() < curconvlim) {
			break;
		}

		tries += 1;

		//if (!RockingBC_dll) {
			int solvesuc = dUd_dW.Solve(-1.0 * Ud, DW);
		//}
		//else {
			
		//CHANGE FOR DLL CREATION
			
		//	Eigen::MatrixXd dUd_dW_2 = Eigen::MatrixXd(Ud.Size(), W.Size());
		//	for (size_t i = 0; i != Ud.Size(); i++) {
		//		for (size_t j = 0; j != W.Size(); j++) {
		//			dUd_dW_2(i, j) = dUd_dW(i, j);
		//		}
		//	}
		//	Eigen::VectorXd Ud_2 = Eigen::VectorXd(Ud.Size());
		//	for (size_t i = 0; i != Ud.Size(); i++) {
		//		Ud_2(i) = Ud(i);
		//	}

		//	Eigen::VectorXd DW_2 = -dUd_dW_2.colPivHouseholderQr().solve(Ud_2);

		//	for (size_t i = 0; i != W.Size(); i++) {
		//		DW(i) = DW_2(i);
		//	}
		//}

		if (tries < maxtries / curaf / curaf / curaf) {
			W += curaf*DW;
		}
		else {
			itry += 1;

			if (curaf <= aflim) {
				std::cout << "Maximum tries reached at NL_solve" << std::endl;

				NLsolvesuccess = -1;
				//std::exit(25);
				break;
			}

			if (itry == 1) {
				W = Winit;
				curaf = af * 0.5;
				curconvlim = convlimmult * convlim;
			}
			else if (itry == 2) {
				W.Zero();
				curaf = af * 0.5;
				curconvlim = convlimmult * convlim;
			}
			else {
				W.Zero();
				curaf = curaf * 0.5;
				curconvlim = convlimmult * curconvlim;
			}
			
		}

	}
	return NLsolvesuccess;

}

void
RockingBC::writedbgfile() {
	std::ofstream NLFfile = std::ofstream("NLsolvefailure.txt");

	if (useUelNM) {
		Ys_com = interval_join(Ysi_com);
		S_com = interval_join(Si_com);
		Yup_com = interval_join(Yupi_com);
		Up_com = interval_join(Upi_com);
	}
	NLFfile << "ue:" << &ue << std::endl;
	NLFfile << "W:" << &W << std::endl;
	NLFfile << "Yw:" << &Yw << std::endl;
	NLFfile << "Yw_len:" << Yw.Size() << std::endl;
	NLFfile << "E:" << std::setprecision(16) << E << std::endl;
	NLFfile << "nu:" << nu << std::endl;
	NLFfile << "ey:" << std::setprecision(16) << ey << std::endl;
	NLFfile << "L:" << std::setprecision(16) << L << std::endl;
	NLFfile << "b:" << std::setprecision(16) << b << std::endl;
	NLFfile << "w:" << std::setprecision(16) << w << std::endl;
	NLFfile << "Yup_com:" << &Yup_com << std::endl;
	NLFfile << "Up_com:" << &Up_com << std::endl;
	NLFfile << "Up_com_len:" << Up_com.Size() << std::endl;
	NLFfile << "Ys_com:" << &Ys_com << std::endl;
	NLFfile << "S_com:" << &S_com << std::endl;
	NLFfile << "S_com_len:" << S_com.Size() << std::endl;
	NLFfile << "Fn_com:" << &Fn_com << std::endl;
	NLFfile << "mu:" << std::setprecision(16) << mu << std::endl;
	NLFfile << "sL_com:" << std::setprecision(16) << sLcommit << std::endl;
	NLFfile << "beta_Dt:" << std::setprecision(16) << beta_Dt << std::endl;
	NLFfile << "useshear:" << useshear << std::endl;
	NLFfile << "blevery:" << blevery << std::endl;
	NLFfile << "slidmode:" << slidmode_com << std::endl;
}
