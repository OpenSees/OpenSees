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
//		response of inelastic deformable rocking bodies.” Earthquake Engineering and Structural Dynamics,
//      49(11), 1072-1094.

#include <RockingBC.h>
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
		Yup_file.open(fstr + "_Yup.txt");
		Up_file.open(fstr + "_Up.txt");
		Ys_file.open(fstr + "_Ys.txt");
		S_file.open(fstr + "_S.txt");

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
	std::ofstream NLFfile("NLsolvefailure.txt");

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

// SISfuncs

double RockingBC::YMXLOGYMX(double y, double p)
{
	double ymx{ p - y };
	if (std::fabs(ymx)<SISfunclim)
	{
		return 0.0;
	}
	else
	{
		return ymx*std::log(std::fabs(ymx));
	}
}

double RockingBC::OMXYLOGOMXYOXY(double yp)
{
	if (std::fabs(yp)<SISfunclim)
	{
		return -1.0;
	}
	else if (std::fabs(yp - 1) < SISfunclim)
	{
		return 0.0;
	}
	else
	{
		return (1.0 - yp)*std::log1p(-yp) / (yp);
	}
}

double RockingBC::J2(double yp)
{
	if (std::fabs(yp)<SISfunclim)
	{
		return 0.5;
	}
	else if (std::fabs(yp - 1) < SISfunclim)
	{
		return 1.0;
	}
	else
	{
		return (OMXYLOGOMXYOXY(yp) + 1.0) / (yp);
	}
}

double RockingBC::OMXATANYMOOXMO(double y, double p)
{
	if (std::fabs(y-1)<SISfunclim)
	{
		return 0.0;
	}
	else
	{ 
		return (1.0 - y)*std::atan((p - 1.0) / (y - 1.0));
	}
}

double RockingBC::OMYLOGSQ(double y, double p)
{
	if (std::fabs(p - 1)<SISfunclim)
	{
		return 0.0;
	}
	else
	{
		return (1.0 - p)*std::log((y - 1.0)*(y - 1.0) + (p - 1.0)*(p - 1.0));
	}
}

double RockingBC::I_FA(double y, double p)
{
	double FA1 = 2 * YMXLOGYMX(y, p);
	double FA2 = -OMXYLOGOMXYOXY(y*p)*p / 3. * (2 * y*y * p*p + 5 * y*p - 1);
	double FA3 = OMXYLOGOMXYOXY(-y*p)*p / 3. * (1 + y*p)*(2 * y*p - 1);
	double FA4 = 4. / 3. * y*p*p;

	return FA1 + FA2 + FA3 + FA4;
}

double RockingBC::J_FA(double y, double p)
{
	double FA1 = +YMXLOGYMX(y, p)*(p + y);
	double FA2 = -p*p / 6. * (OMXYLOGOMXYOXY(y*p) + YMXLOGYMX(y*p, 1.)*(3. * y*p + 7.) + J2(y*p));
	double FA3 = -p*p / 6. * (OMXYLOGOMXYOXY(-y*p) + YMXLOGYMX(-y*p, 1.)*(3. * y*p + 1.) + J2(-y*p));
	double FA4 = y*p*p*p + p*p / 3. - p*y;

	return FA1 + FA2 + FA3 + FA4;
}

double RockingBC::pImJ_FA(double y, double p)
{
	double FA1 = (p - y)*YMXLOGYMX(y, p);
	double FA2 = -OMXYLOGOMXYOXY(y*p)*p*p / 3. * (2. * y*y * p*p + 5. * y*p - 3. / 2.) + p*p / 6. * YMXLOGYMX(y*p, 1.)*(3 * y*p + 7) + p*p / 6. * J2(y*p);
	double FA3 = OMXYLOGOMXYOXY(-y*p)*p*p / 3. * ((1 + y*p)*(2 * y*p - 1) + 1. / 2.) + p*p / 6. * YMXLOGYMX(-y*p, 1.)*(3 * y*p + 1) + p*p / 6. * J2(-y*p);
	double FA4 = 1. / 3. * y*p*p*p - p*p / 3. + p*y;

	return FA1 + FA2 + FA3 + FA4;
}

double RockingBC::pImJ_FA_nochecks(double y, double p)
{
	double FA1 = (p - y)*(p - y) * log(abs(p - y));
	double FA2 = -(1. - y*p)*log1p(-y*p) / (y*p)*p*p / 3. * (2. * y*y * p*p + 5 * y*p - 3. / 2.) + p*p / 6. * (1 - y*p)*log(1 - y*p)*(3. * y*p + 7.) + p*p / 6. * ((1. - y*p)*log1p(-y*p) / (y*p) + 1.) / (y*p);
	double FA3 = -(1. + y*p)*log1p(y*p) / (y*p)*p*p / 3. * ((1 + y*p)*(2 * y*p - 1) + 1. / 2.) + p*p / 6. * (1 + y*p)*log(1 + y*p)*(3. * y*p + 1.) + p*p / 6. * ((1. + y*p)*log1p(y*p) / (y*p) - 1.) / (y*p);
	double FA4 = 1. / 3. * y*p*p*p - p*p / 3. + p*y;

	return FA1 + FA2 + FA3 + FA4;
}

double RockingBC::I_FB(double y, double p)
{
	static const double log2{ std::log(2) };
	static const double pi{ std::atan(1) * 4 };

	double FB1 = -4. * OMXATANYMOOXMO(y, p) - 2. * OMYLOGSQ(y, p);
	double FB2 = 4. * OMXATANYMOOXMO(-y, -p) + 2. * OMYLOGSQ(-y, -p);

	double FB3 = 3. / 2. * p*p * ((1 + y)*YMXLOGYMX(-y, 1.) - (1 - y)*YMXLOGYMX(y, 1.));
	double FB4 = y*((1 + p)*(1 + p) * YMXLOGYMX(-p, 1.) + (1 - p)*(1 - p) * YMXLOGYMX(p, 1.));

	double FB5 = std::log((y - 1)*(y - 1) + 4)*p / 4. * (3 * p*y*y - 6 * p*y + 3 * p - 8);
	double FB6 = -std::log((y + 1)*(y + 1) + 4)*p / 4. * (3 * p*y*y + 6 * p*y + 3 * p + 8);
	double FB7 = +std::log((p - 1)*(p - 1) + 4)*(2 - y / 2. - 2 * p + 3. / 2. * y*p - 3. / 2. * y*p*p + y*p*p*p / 2.);
	double FB8 = -std::log((p + 1)*(p + 1) + 4)*(2 + y / 2. + 2 * p + 3. / 2. * y*p + 3. / 2. * y*p*p + y*p*p*p / 2.);

	double FB9 = std::atan(y / 2. - 1. / 2.)*p*(3 * p + 2)*(y - 1);
	double FB10 = -std::atan(y / 2. + 1. / 2.)*p*(3 * p - 2)*(y + 1);
	double FB11 = -std::atan(p / 2. - 1. / 2.)*((3 * y + 1)*(-5 + 2 * p - p*p) + 8 * (1 + y));
	double FB12 = -std::atan(p / 2. + 1. / 2.)*((3 * y - 1)*(5 + 2 * p + p*p) + 8 * (1 - y));

	double FB13 = +2 * p*(6. * log2 - pi) + 3. * y*p*p * (2 * log2 + pi + 1);

	return FB1 + FB2 + FB3 + FB4 + FB5 + FB6 + FB7 + FB8 + FB9 + FB10 + FB11 + FB12 + FB13;
}

double RockingBC::J_FB(double y, double p)
{
	static const double log2{ std::log(2) };
	static const double pi{ std::atan(1) * 4 };

	double FB1 = -4. * OMXATANYMOOXMO(y, p) + (1 - y)*OMYLOGSQ(p, y) - (p + 1)*OMYLOGSQ(y, p);
	double FB2 = -4. * OMXATANYMOOXMO(-y, -p) + (1 + y)*OMYLOGSQ(-p, -y) + (p - 1)*OMYLOGSQ(-y, -p);

	double FB3 = p*p*p * ((1 + y)*YMXLOGYMX(-y, 1.) - (1 - y)*YMXLOGYMX(y, 1.));
	double FB4 = y / 4. * ((3. * p - 1)*(1 + p)*(1 + p) * YMXLOGYMX(-p, 1.) + (3 * p + 1)*(1 - p)*(1 - p) * YMXLOGYMX(p, 1.));

	double FB5 = std::log((y - 1)*(y - 1) + 4.)*p*p * (p*y*y - 2 * p*y + p - 2.) / 2.;
	double FB6 = -std::log((y + 1)*(y + 1) + 4.)*p*p * (p*y*y + 2 * p*y + p + 2.) / 2.;
	double FB7 = std::log((p - 1)*(p - 1) + 4.)*(-1. / 3. - p*p + 15. / 8. * y + 3. / 4. * y*p*p - y*p*p*p + 3. / 8. * y*pow(p,4));
	double FB8 = std::log((p + 1)*(p + 1) + 4.)*(-1. / 3. - p*p - 15. / 8. * y - 3. / 4. * y*p*p - y*p*p*p - 3. / 8. * y*pow(p,4));

	double FB9 = std::atan(y / 2. - 1. / 2.)*p*p * (2 * p + 1)*(y - 1);
	double FB10 = -std::atan(y / 2. + 1. / 2.)*p*p * (2 * p - 1)*(y + 1);
	double FB11 = -std::atan(p / 2. - 1. / 2.)*((y + 1. / 3.)*(-13. + 3 * p*p - 2 * p*p*p) + 8. * (1 + y));
	double FB12 = -std::atan(p / 2. + 1. / 2.)*((y - 1. / 3.)*(-13. + 3 * p*p + 2 * p*p*p) - 8. * (1 - y));

	double FB13 = 2. * y*(2. * log2 + pi + 1.)*p*p*p + (6. * log2 - pi + 2. / 3.)*p*p - (2. * y*p);

	return FB1 + FB2 + FB3 + FB4 + FB5 + FB6 + FB7 + FB8 + FB9 + FB10 + FB11 + FB12 + FB13;

}

double RockingBC::pImJ_FB(double y, double p)
{
	static const double log2{ std::log(2) };
	static const double pi{ std::atan(1) * 4 };

	double FB1 = 4. * (1 - p)*OMXATANYMOOXMO(y, p) - (1. - y)*OMYLOGSQ(p, y) + (1 - p)*OMYLOGSQ(y, p);
	double FB2 = 4. * (1 + p)*OMXATANYMOOXMO(-y, -p) - (1. + y)*OMYLOGSQ(-p, -y) + (1 + p)*OMYLOGSQ(-y, -p);

	double FB3 = p*p*p * ((1 + y)*YMXLOGYMX(-y, 1) - (1 - y)*YMXLOGYMX(y, 1)) / 2.;
	double FB4 = y*(pow((1 + p),3) * YMXLOGYMX(-p, 1) - pow((1 - p),3) * YMXLOGYMX(p, 1)) / 4.;

	double FB5 = log((y - 1)*(y - 1) + 4)*p*p / 4. * (p*y*y - 2 * p*y + p - 4);
	double FB6 = -log((y + 1)*(y + 1) + 4)*p*p / 4. * (p*y*y + 2 * p*y + p + 4);
	double FB7 = log((p - 1)*(p - 1) + 4)*(+1. / 3. - 15. / 8. * y + 2. * p - y*p / 2. - p*p + 3. / 4 * y*p*p - 1. / 2. * y*p*p*p + 1. / 8. * y*p*p*p*p);
	double FB8 = -log((p + 1)*(p + 1) + 4)*(-1. / 3. - 15. / 8. * y + 2. * p + y*p / 2. + p*p + 3. / 4 * y*p*p + 1. / 2. * y*p*p*p + 1. / 8. * y*p*p*p*p);

	double FB9 = atan(y / 2. - 1. / 2.)*p*p * (1 + p)*(y - 1);
	double FB10 = atan(y / 2. + 1. / 2.)*p*p * (1 - p)*(y + 1);
	double FB11 = atan(p / 2. - 1. / 2.)*(1 - p)*(2. * p - 15. * y + 6. * p*y - 3. * p*p * y - p*p + 11.) / 3.;
	double FB12 = -atan(p / 2. + 1. / 2.)*(1 + p)*(-2. * p + 15. * y + 6. * p*y + 3. * p*p * y - p*p + 11.) / 3.;

	double FB13 = p*p * (6. * log2 - pi - 2. / 3.) + y*p*p*p * (2. * log2 + pi + 1.) + (2 * y*p);

	return FB1 + FB2 + FB3 + FB4 + FB5 + FB6 + FB7 + FB8 + FB9 + FB10 + FB11 + FB12 + FB13;
}

double RockingBC::pImJ_FB_nochecks(double y, double p)
{
	static const double log2{ std::log(2) };
	static const double pi{ std::atan(1) * 4 };

	double FB1 = 4. * (1 - p)*(1 - y)*atan((p - 1) / (y - 1)) + ((1 - p)*(1 - p) - (1 - y)*(1 - y))*log((y - 1)*(y - 1) + (p - 1)*(p - 1));
	double FB2 = 4. * (1 + p)*(1 + y)*atan((p + 1) / (y + 1)) + ((1 + p)*(1 + p) - (1 + y)*(1 + y))*log((y + 1)*(y + 1) + (p + 1)*(p + 1));

	double FB3 = p*p*p * ((1 + y)*(1 + y) * log(1 + y) - (1 - y)*(1 - y) * log(1 - y)) / 2.;
	double FB4 = y*(pow((1 + p),4) * log(1 + p) - pow((1 - p),4) * log(1 - p)) / 4.;

	double FB5 = log((y - 1)*(y - 1) + 4)*p*p / 4. * (p*y*y - 2 * p*y + p - 4);
	double FB6 = -log((y + 1)*(y + 1) + 4)*p*p / 4. * (p*y*y + 2 * p*y + p + 4);
	double FB7 = log((p - 1)*(p - 1) + 4)*(+1. / 3. - 15. / 8. * y + 2. * p - y*p / 2. - p*p + 3. / 4 * y*p*p - 1. / 2. * y*p*p*p + 1. / 8. * y*p*p*p*p);
	double FB8 = -log((p + 1)*(p + 1) + 4)*(-1. / 3. - 15. / 8. * y + 2. * p + y*p / 2. + p*p + 3. / 4 * y*p*p + 1. / 2. * y*p*p*p + 1. / 8. * y*p*p*p*p);

	double FB9 = atan(y / 2. - 1. / 2.)*p*p * (1 + p)*(y - 1);
	double FB10 = atan(y / 2. + 1. / 2.)*p*p * (1 - p)*(y + 1);
	double FB11 = atan(p / 2. - 1. / 2.)*(1 - p)*(2. * p - 15. * y + 6. * p*y - 3. * p*p * y - p*p + 11.) / 3.;
	double FB12 = -atan(p / 2. + 1. / 2.)*(1 + p)*(-2. * p + 15. * y + 6. * p*y + 3. * p*p * y - p*p + 11.) / 3.;

	double FB13 = p*p * (6. * log2 - pi - 2. / 3.) + y*p*p*p * (2. * log2 + pi + 1.) + (2 * y*p);

	return FB1 + FB2 + FB3 + FB4 + FB5 + FB6 + FB7 + FB8 + FB9 + FB10 + FB11 + FB12 + FB13;
}

double RockingBC::I_FP(double y, double p)
{
	static const double q22{ -1.23991592 }, q42{ 1.08897876 }, q44{ -9.726553 }, q62{ -1.50465044 }, q64{ 18.273236 }, q66{ -38.99971412 }, q82{ 0.74180336 }, q84{ -9.64366612 }, q86{ 22.03387365 }, q88{ -13.05630027 };
	static const double q31{ 0.74952005 }, q51{ -0.08175407 }, q53{ 5.10578057 }, q71{ 0.04700608 }, q73{ -2.9709584 }, q75{ 9.15391675 };

	double FP1 = ((q88*pow(y,8)) / 9 + (q86*pow(y,6)) / 9 + (q84*pow(y,4)) / 9 + (q82*y*y) / 9 - q82 / 27 - q84 / 45 - q86 / 63 - q88 / 81)*pow(p,9) + ((q86*pow(y,8)) / 7 + (q66*pow(y,6)) / 7 + (q64*pow(y,4)) / 7 + (q62*y*y) / 7 - q62 / 21 - q64 / 35 - q66 / 49 - q86 / 63)*pow(p,7) + ((q84*pow(y,8)) / 5 + (q64*pow(y,6)) / 5 + (q44*pow(y,4)) / 5 + (q42*y*y) / 5 - q42 / 15 - q44 / 25 - q64 / 35 - q84 / 45)*pow(p,5) + ((q82*pow(y,8)) / 3 + (q62*pow(y,6)) / 3 + (q42*pow(y,4)) / 3 + (q22*y*y) / 3 - q22 / 9 - q42 / 15 - q62 / 21 - q82 / 27)*p*p*p + ((-q82 / 3 - q84 / 5 - q86 / 7 - q88 / 9)*pow(y,8) + (-q62 / 3 - q64 / 5 - q66 / 7 - q86 / 9)*pow(y,6) + (-q42 / 3 - q44 / 5 - q64 / 7 - q84 / 9)*pow(y,4) + (-q22 / 3 - q42 / 5 - q62 / 7 - q82 / 9)*y*y + q22 / 9 + (2 * q42) / 15 + q44 / 25 + (2 * q62) / 21 + (2 * q64) / 35 + q66 / 49 + (2 * q82) / 27 + (2 * q84) / 45 + (2 * q86) / 63 + q88 / 81)*p;
	double FP2 = ((q71*y) / 8 - (pow(y,7) * (3 * q71 + (9 * q73) / 5 + (9 * q75) / 7)) / 8 + (q73*y*y*y) / 8 + (q75*pow(y,5)) / 8)*pow(p,8) + ((q51*y) / 6 - (pow(y,5) * ((7 * q51) / 3 + (7 * q53) / 5 + (7 * q75) / 9)) / 6 + (q53*y*y*y) / 6 + (q75*pow(y,7)) / 6)*pow(p,6) + ((q31*y) / 4 - (y*y*y * ((5 * q31) / 3 + (5 * q53) / 7 + (5 * q73) / 9)) / 4 + (q53*pow(y,5)) / 4 + (q73*pow(y,7)) / 4)*pow(p,4) + ((q31*y*y*y) / 2 - (y*((3 * q31) / 5 + (3 * q51) / 7 + q71 / 3)) / 2 + (q51*pow(y,5)) / 2 + (q71*pow(y,7)) / 2)*p*p;

	return FP1 + FP2;
}

double RockingBC::J_FP(double y, double p)
{
	static const double q22{ -1.23991592 }, q42{ 1.08897876 }, q44{ -9.726553 }, q62{ -1.50465044 }, q64{ 18.273236 }, q66{ -38.99971412 }, q82{ 0.74180336 }, q84{ -9.64366612 }, q86{ 22.03387365 }, q88{ -13.05630027 };
	static const double q31{ 0.74952005 }, q51{ -0.08175407 }, q53{ 5.10578057 }, q71{ 0.04700608 }, q73{ -2.9709584 }, q75{ 9.15391675 };

	double FP1 = ((q88*pow(y,8)) / 10 + (q86*pow(y,6)) / 10 + (q84*pow(y,4)) / 10 + (q82*y*y) / 10 - q82 / 30 - q84 / 50 - q86 / 70 - q88 / 90)*pow(p,10) + ((q86*pow(y,8)) / 8 + (q66*pow(y,6)) / 8 + (q64*pow(y,4)) / 8 + (q62*y*y) / 8 - q62 / 24 - q64 / 40 - q66 / 56 - q86 / 72)*pow(p,8) + ((q84*pow(y,8)) / 6 + (q64*pow(y,6)) / 6 + (q44*pow(y,4)) / 6 + (q42*y*y) / 6 - q42 / 18 - q44 / 30 - q64 / 42 - q84 / 54)*pow(p,6) + ((q82*pow(y,8)) / 4 + (q62*pow(y,6)) / 4 + (q42*pow(y,4)) / 4 + (q22*y*y) / 4 - q22 / 12 - q42 / 20 - q62 / 28 - q82 / 36)*pow(p,4) + (q22 / 18 + q42 / 15 + q44 / 50 + q62 / 21 + q64 / 35 + q66 / 98 + q82 / 27 + q84 / 45 + q86 / 63 + q88 / 162 - (q22*y*y) / 6 - (q42*y*y) / 10 - (q42*pow(y,4)) / 6 - (q44*pow(y,4)) / 10 - (q62*y*y) / 14 - (q62*pow(y,6)) / 6 - (q64*pow(y,4)) / 14 - (q64*pow(y,6)) / 10 - (q66*pow(y,6)) / 14 - (q82*y*y) / 18 - (q84*pow(y,4)) / 18 - (q82*pow(y,8)) / 6 - (q84*pow(y,8)) / 10 - (q86*pow(y,6)) / 18 - (q86*pow(y,8)) / 14 - (q88*pow(y,8)) / 18)*p*p;
	double FP2 = (y*p*p*p * (p*p - 1)*(63 * q31 + 45 * q51 + 35 * q71 - 105 * q31*y*y - 105 * q51*pow(y,4) - 105 * q71*pow(y,6) + 45 * q51*p*p + 35 * q71*p*p + 35 * q71*pow(p,4) - 105 * q51*pow(y,4) * p*p + 45 * q53*y*y * p*p - 63 * q53*pow(y,4) * p*p + 35 * q73*y*y * p*p - 105 * q71*pow(y,6) * p*p + 35 * q73*y*y * pow(p,4) - 105 * q71*pow(y,6) * pow(p,4) - 63 * q73*pow(y,6) * p*p - 63 * q73*pow(y,6) * pow(p,4) + 35 * q75*pow(y,4) * pow(p,4) - 45 * q75*pow(y,6) * pow(p,4))) / 315;

	return FP1 + FP2;
}

double RockingBC::I_FP_alt(double y, double p)
{
	static const double q22{ -1.23991592 }, q42{ 1.08897876 }, q44{ -9.726553 }, q62{ -1.50465044 }, q64{ 18.273236 }, q66{ -38.99971412 }, q82{ 0.74180336 }, q84{ -9.64366612 }, q86{ 22.03387365 }, q88{ -13.05630027 };
	static const double q31{ 0.74952005 }, q51{ -0.08175407 }, q53{ 5.10578057 }, q71{ 0.04700608 }, q73{ -2.9709584 }, q75{ 9.15391675 };
	static const double a1 = -q82 / 27 - q84 / 45 - q86 / 63 - q88 / 81;
	static const double a2 = -q62 / 21 - q64 / 35 - q66 / 49 - q86 / 63;
	static const double a3 = -q42 / 15 - q44 / 25 - q64 / 35 - q84 / 45;
	static const double a4 = -q22 / 9 - q42 / 15 - q62 / 21 - q82 / 27;
	static const double a5 = -q82 / 3 - q84 / 5 - q86 / 7 - q88 / 9;
	static const double a6 = -q62 / 3 - q64 / 5 - q66 / 7 - q86 / 9;
	static const double a7 = -q42 / 3 - q44 / 5 - q64 / 7 - q84 / 9;
	static const double a8 = -q22 / 3 - q42 / 5 - q62 / 7 - q82 / 9;
	static const double a9 = +q22 / 9 + (2 * q42) / 15 + q44 / 25 + (2 * q62) / 21 + (2 * q64) / 35 + q66 / 49 + (2 * q82) / 27 + (2 * q84) / 45 + (2 * q86) / 63 + q88 / 81;
	static const double a10 = 3 * q71 + (9 * q73) / 5 + (9 * q75) / 7;
	static const double a11 = (7 * q51) / 3 + (7 * q53) / 5 + (7 * q75) / 9;
	static const double a12 = (5 * q31) / 3 + (5 * q53) / 7 + (5 * q73) / 9;
	static const double a13 = (3 * q31) / 5 + (3 * q51) / 7 + q71 / 3;
	
	double FP1 = ((q88*pow(y,8)) / 9 + (q86*pow(y,6)) / 9 + (q84*pow(y,4)) / 9 + (q82*y*y) / 9 + a1)*pow(p,9) + ((q86*pow(y,8)) / 7 + (q66*pow(y,6)) / 7 + (q64*pow(y,4)) / 7 + (q62*y*y) / 7 +a2)*pow(p,7) + ((q84*pow(y,8)) / 5 + (q64*pow(y,6)) / 5 + (q44*pow(y,4)) / 5 + (q42*y*y) / 5 +a3)*pow(p,5) + ((q82*pow(y,8)) / 3 + (q62*pow(y,6)) / 3 + (q42*pow(y,4)) / 3 + (q22*y*y) / 3 +a4)*p*p*p + ((a5)*pow(y,8) + (a6)*pow(y,6) + (a7)*pow(y,4) + (a8)*y*y +a9)*p;
	double FP2 = ((q71*y) / 8 - (pow(y,7) * (a10)) / 8 + (q73*y*y*y) / 8 + (q75*pow(y,5)) / 8)*pow(p,8) + ((q51*y) / 6 - (pow(y,5) * (a11)) / 6 + (q53*y*y*y) / 6 + (q75*pow(y,7)) / 6)*pow(p,6) + ((q31*y) / 4 - (y*y*y * (a12)) / 4 + (q53*pow(y,5)) / 4 + (q73*pow(y,7)) / 4)*pow(p,4) + ((q31*y*y*y) / 2 - (y*(a13)) / 2 + (q51*pow(y,5)) / 2 + (q71*pow(y,7)) / 2)*p*p;

	return FP1 + FP2;
}

double RockingBC::pImJ_FP(double y, double p)
{
	static const double q22{ -1.23991592 }, q42{ 1.08897876 }, q44{ -9.726553 }, q62{ -1.50465044 }, q64{ 18.273236 }, q66{ -38.99971412 }, q82{ 0.74180336 }, q84{ -9.64366612 }, q86{ 22.03387365 }, q88{ -13.05630027 };
	static const double q31{ 0.74952005 }, q51{ -0.08175407 }, q53{ 5.10578057 }, q71{ 0.04700608 }, q73{ -2.9709584 }, q75{ 9.15391675 };

	double FP1 = -p*p * ((p*p * q22) / 36 - q42 / 15 - q44 / 50 - q62 / 21 - q64 / 35 - q66 / 98 - q82 / 27 - q84 / 45 - q86 / 63 - q88 / 162 - q22 / 18 + (p*p * q42) / 60 + (pow(p,4) * q42) / 90 + (pow(p,4) * q44) / 150 + (p*p * q62) / 84 + (pow(p,4) * q64) / 210 + (pow(p,6) * q62) / 168 + (pow(p,6) * q64) / 280 + (pow(p,6) * q66) / 392 + (p*p * q82) / 108 + (pow(p,4) * q84) / 270 + (pow(p,8) * q82) / 270 + (pow(p,6) * q86) / 504 + (pow(p,8) * q84) / 450 + (pow(p,8) * q86) / 630 + (pow(p,8) * q88) / 810 + (q22*y*y) / 6 + (q42*y*y) / 10 + (q42*pow(y,4)) / 6 + (q44*pow(y,4)) / 10 + (q62*y*y) / 14 + (q62*pow(y,6)) / 6 + (q64*pow(y,4)) / 14 + (q64*pow(y,6)) / 10 + (q66*pow(y,6)) / 14 + (q82*y*y) / 18 + (q84*pow(y,4)) / 18 + (q82*pow(y,8)) / 6 + (q84*pow(y,8)) / 10 + (q86*pow(y,6)) / 18 + (q86*pow(y,8)) / 14 + (q88*pow(y,8)) / 18 - (p*p * q22*y*y) / 12 - (p*p * q42*pow(y,4)) / 12 - (pow(p,4) * q42*y*y) / 30 - (pow(p,4) * q44*pow(y,4)) / 30 - (p*p * q62*pow(y,6)) / 12 - (pow(p,6) * q62*y*y) / 56 - (pow(p,4) * q64*pow(y,6)) / 30 - (pow(p,6) * q64*pow(y,4)) / 56 - (pow(p,6) * q66*pow(y,6)) / 56 - (p*p * q82*pow(y,8)) / 12 - (pow(p,8) * q82*y*y) / 90 - (pow(p,4) * q84*pow(y,8)) / 30 - (pow(p,8) * q84*pow(y,4)) / 90 - (pow(p,6) * q86*pow(y,8)) / 56 - (pow(p,8) * q86*pow(y,6)) / 90 - (pow(p,8) * q88*pow(y,8)) / 90);
	double FP2 = -(p*p*p * y*(756 * q31 + 540 * q51 + 420 * q71 - 378 * p*p * q31 - 180 * pow(p,4) * q51 - 105 * pow(p,6) * q71 - 1260 * q31*y*y - 1260 * q51*pow(y,4) - 1260 * q71*pow(y,6) + 630 * p*p * q31*y*y + 270 * p*p * q53*y*y - 378 * p*p * q53*pow(y,4) + 420 * pow(p,4) * q51*pow(y,4) - 180 * pow(p,4) * q53*y*y + 252 * pow(p,4) * q53*pow(y,4) + 210 * p*p * q73*y*y - 378 * p*p * q73*pow(y,6) - 105 * pow(p,6) * q73*y*y + 140 * pow(p,4) * q75*pow(y,4) + 315 * pow(p,6) * q71*pow(y,6) - 180 * pow(p,4) * q75*pow(y,6) + 189 * pow(p,6) * q73*pow(y,6) - 105 * pow(p,6) * q75*pow(y,4) + 135 * pow(p,6) * q75*pow(y,6))) / 7560;

	return FP1 + FP2;
}

double RockingBC::I_calc(double y, double r)
{
	static const double pi{ std::atan(1) * 4 };
	static const double A = -1 / pi;
	static const double B = -0.19532775;

	return A*I_FA(y, r) + B*I_FB(y, r) + I_FP(y, r);
}

double RockingBC::J_calc(double y, double r)
{
	static const double pi{ std::atan(1.) * 4 };
	static const double A = -1 / pi;
	static const double B = -0.19532775;

	return A*J_FA(y, r) + B*J_FB(y, r) + J_FP(y, r);
}

double RockingBC::FAa(double y, double p)
{
	static const double pi{ std::atan(1) * 4 };
	static const double A = -1 / pi;
	return A*2.0*std::log(std::fabs(p - y));
}

double RockingBC::dFAa_dp(double y, double p)
{
	static const double pi{ std::atan(1) * 4 };
	static const double A = -1 / pi;
	return -A*2.0/(y-p);
}

double RockingBC::FA(double y, double p)
{
	double FA1 = 2. * std::log(fabs(p - y));
	double FA2 = 2. * std::log1p(-y*p)*(-1. + y*p + y*y * p*p);
	double FA3 = 2. * std::log1p(y*p)*(-y*p - y*y * p*p);
	double FA4 = 2. * (2 * y*p + 1);

	return FA1 + FA2 + FA3 + FA4;
}

double RockingBC::D_FA(double y, double p)
{
	double FA1 = 2. / (y - p);
	double FA2 = 2. * std::log1p(-y*p)*(2 * y*p*p + p) + (2 * p*(y*y * p*p + y*p - 1)) / (y*p - 1);
	double FA3 = -2. * std::log1p(y*p)*(2 * y*p*p + p) - (2 * p*(y*y * p*p + y*p)) / (y*p + 1);
	double FA4 = 4. * p;

	return FA1 + FA2 + FA3 + FA4;
}

double RockingBC::FB(double y, double p)
{
	static const double log2{ std::log(2) };
	static const double pi{ std::atan(1) * 4 };

	double FB1 = 2. * log((1. - y)*(1. - y) + (1. - p)*(1. - p));
	double FB2 = 2. * log((1. + y)*(1. + y) + (1. + p)*(1. + p));

	double FB3 = 3. * p*((1. + y)*YMXLOGYMX(-y, 1.) - (1. - y)*YMXLOGYMX(y, 1.));
	double FB4 = 3. * y*((1. + p)*YMXLOGYMX(-p, 1.) - (1. - p)*YMXLOGYMX(p, 1.));

	double FB5 = (3. / 2. * p - 3 * y*p + 3. / 2. * p*y*y - 2.)*log((1. - y)*(1. - y) + 4.);
	double FB6 = (-3. / 2. * p - 3 * y*p - 3. / 2. * p*y*y - 2.)*log((1. + y)*(1. + y) + 4.);
	double FB7 = (3. / 2. * y - 3 * p*y + 3. / 2. * y*p*p - 2.)*log((1. - p)*(1. - p) + 4.);
	double FB8 = (-3. / 2. * y - 3 * p*y - 3. / 2. * y*p*p - 2.)*log((1. + p)*(1. + p) + 4.);

	double FB9 = 2. * (1. - y + 3. * p - 3. * y*p)*atan((1. - y) / 2.);
	double FB10 = 2. * (1. + y - 3. * p - 3. * y*p)*atan((1. + y) / 2.);
	double FB11 = 2. * (1. - p + 3. * y - 3. * y*p)*atan((1. - p) / 2.);
	double FB12 = 2. * (1. + p - 3. * y - 3. * y*p)*atan((1. + p) / 2.);

	double FB13 = (6. * y*p*(pi + 2. * log2 + 1.) - 2. * (pi - 6. * log2 - 2.));

	return FB1 + FB2 + FB3 + FB4 + FB5 + FB6 + FB7 + FB8 + FB9 + FB10 + FB11 + FB12 + FB13;
}

double RockingBC::D_FB(double y, double p)
{
	static const double log2{ std::log(2) };
	static const double pi{ std::atan(1) * 4 };

	double FB1 = (2. * (2. * y - 2.)) / ((y - 1.)*(y - 1.) + (p - 1.)*(p - 1.));
	double FB2 = (2. * (2. * y + 2.)) / ((y + 1.)*(y + 1.) + (p + 1.)*(p + 1.));

	double FB3 = 6. * p*(YMXLOGYMX(-y, 1.) + YMXLOGYMX(y, 1.) + 1.);
	double FB4 = 3. * YMXLOGYMX(-p, 1.)*(p + 1.) - 3. * YMXLOGYMX(p, 1.)*(1. - p);

	double FB5 = ((2. * y - 2.)*((3. * p*y*y) / 2. - 3. * p*y + (3. * p) / 2. - 2.)) / ((y - 1)*(y - 1) + 4.) - log((y - 1)*(y - 1) + 4.)*(3. * p - 3. * y*p);
	double FB6 = -log((y + 1)*(y + 1) + 4.)*(3. * p + 3. * y*p) - ((2. * y + 2.)*((3. * p*y*y) / 2. + 3. * p*y + (3. * p) / 2. + 2.)) / ((y + 1)*(y + 1) + 4.);
	double FB7 = log((p - 1)*(p - 1) + 4.)*((3. * p*p) / 2. - 3. * p + 3. / 2.);
	double FB8 = -log((p + 1)*(p + 1) + 4.)*((3. * p*p) / 2. + 3. * p + 3. / 2.);

	double FB9 = atan(y / 2. - 1. / 2.)*(6. * p + 2.) + (4. * (3. * p + 1.)*(y - 1.)) / (y*y - 2. * y + 5.);
	double FB10 = -atan(y / 2. + 1. / 2.)*(6. * p - 2.) - (4. * (3. * p - 1.)*(y + 1.)) / (y*y + 2. * y + 5.);
	double FB11 = atan(p / 2. - 1. / 2.)*(6. * p - 6.);
	double FB12 = -atan(p / 2. + 1. / 2.)*(6. * p + 6.);

	double FB13 = 6. * p*(2. * log2 + pi + 1.);

	return FB1 + FB2 + FB3 + FB4 + FB5 + FB6 + FB7 + FB8 + FB9 + FB10 + FB11 + FB12 + FB13;
}

double RockingBC::FP(double y, double p)
{
	static const double q22{ -1.23991592 }, q42{ 1.08897876 }, q44{ -9.726553 }, q62{ -1.50465044 }, q64{ 18.273236 }, q66{ -38.99971412 }, q82{ 0.74180336 }, q84{ -9.64366612 }, q86{ 22.03387365 }, q88{ -13.05630027 };
	static const double q31{ 0.74952005 }, q51{ -0.08175407 }, q53{ 5.10578057 }, q71{ 0.04700608 }, q73{ -2.9709584 }, q75{ 9.15391675 };

	double FP1 = q22 / 9 + (2 * q42) / 15 + q44 / 25 + (2 * q62) / 21 + (2 * q64) / 35 + q66 / 49 + (2 * q82) / 27 + (2 * q84) / 45 + (2 * q86) / 63 + q88 / 81 - (y*y + p*p)*(q22 / 3 + q42 / 5 + q62 / 7 + q82 / 9) - (pow(y,4) + pow(p,4))*(q42 / 3 + q44 / 5 + q64 / 7 + q84 / 9) - (pow(y,6) + pow(p,6))*(q62 / 3 + q64 / 5 + q66 / 7 + q86 / 9) - (pow(y,8) + pow(p,8))*(q82 / 3 + q84 / 5 + q86 / 7 + q88 / 9) + q22*y*y * p*p + q44*pow(y,4) * pow(p,4) + q66*pow(y,6) * pow(p,6) + q88*pow(y,8) * pow(p,8) + q42*y*y * p*p * (y*y + p*p) + q62*y*y * p*p * (pow(y,4) + pow(p,4)) + q64*pow(y,4) * pow(p,4) * (y*y + p*p) + q82*y*y * p*p * (pow(y,6) + pow(p,6)) + q84*pow(y,4) * pow(p,4) * (pow(y,4) + pow(p,4)) + q86*pow(y,6) * pow(p,6) * (y*y + p*p);
	double FP2 = q53*y*y*y * p*p*p * (y*y + p*p) - y*y*y * p*p*p * ((5 * q31) / 3 + (5 * q53) / 7 + (5 * q73) / 9) - pow(y,5) * pow(p,5) * ((7 * q51) / 3 + (7 * q53) / 5 + (7 * q75) / 9) - pow(y,7) * pow(p,7) * (3 * q71 + (9 * q73) / 5 + (9 * q75) / 7) - y*p*((3 * q31) / 5 + (3 * q51) / 7 + q71 / 3) + q73*y*y*y * p*p*p * (pow(y,4) + pow(p,4)) + q75*pow(y,5) * pow(p,5) * (y*y + p*p) + q31*y*p*(y*y + p*p) + q51*y*p*(pow(y,4) + pow(p,4)) + q71*y*p*(pow(y,6) + pow(p,6));
	
	return FP1 + FP2;
}

double RockingBC::D_FP(double y, double p)
{
	static const double q22{ -1.23991592 }, q42{ 1.08897876 }, q44{ -9.726553 }, q62{ -1.50465044 }, q64{ 18.273236 }, q66{ -38.99971412 }, q82{ 0.74180336 }, q84{ -9.64366612 }, q86{ 22.03387365 }, q88{ -13.05630027 };
	static const double q31{ 0.74952005 }, q51{ -0.08175407 }, q53{ 5.10578057 }, q71{ 0.04700608 }, q73{ -2.9709584 }, q75{ 9.15391675 };

	double FP1 = 2 * q22*y*p*p - (2 * q42*y) / 5 - (2 * q62*y) / 7 - (2 * q82*y) / 9 - (4 * q42*y*y*y) / 3 - (4 * q44*y*y*y) / 5 - 2 * q62*pow(y,5) - (4 * q64*y*y*y) / 7 - (6 * q64*pow(y,5)) / 5 - (6 * q66*pow(y,5)) / 7 - (4 * q84*y*y*y) / 9 - (8 * q82*pow(y,7)) / 3 - (8 * q84*pow(y,7)) / 5 - (2 * q86*pow(y,5)) / 3 - (8 * q86*pow(y,7)) / 7 - (8 * q88*pow(y,7)) / 9 - (2 * q22*y) / 3 + 2 * q42*y*pow(p,4) + 2 * q62*y*pow(p,6) + 2 * q82*y*pow(p,8) + 4 * q42*y*y*y * p*p + 4 * q44*y*y*y * pow(p,4) + 6 * q62*pow(y,5) * p*p + 4 * q64*y*y*y * pow(p,6) + 6 * q64*pow(y,5) * pow(p,4) + 6 * q66*pow(y,5) * pow(p,6) + 8 * q82*pow(y,7) * p*p + 4 * q84*y*y*y * pow(p,8) + 8 * q84*pow(y,7) * pow(p,4) + 6 * q86*pow(y,5) * pow(p,8) + 8 * q86*pow(y,7) * pow(p,6) + 8 * q88*pow(y,7) * pow(p,8);
	double FP2 = 2 * q31*y*y * p - p*((3 * q31) / 5 + (3 * q51) / 7 + q71 / 3) + 4 * q51*pow(y,4) * p + 6 * q71*pow(y,6) * p + q31*p*(y*y + p*p) + q51*p*(pow(y,4) + pow(p,4)) + q71*p*(pow(y,6) + pow(p,6)) + 2 * q53*pow(y,4) * p*p*p + 4 * q73*pow(y,6) * p*p*p + 2 * q75*pow(y,6) * pow(p,5) - 3 * y*y * p*p*p * ((5 * q31) / 3 + (5 * q53) / 7 + (5 * q73) / 9) - 5 * pow(y,4) * pow(p,5) * ((7 * q51) / 3 + (7 * q53) / 5 + (7 * q75) / 9) - 7 * pow(y,6) * pow(p,7) * (3 * q71 + (9 * q73) / 5 + (9 * q75) / 7) + 3 * q53*y*y * p*p*p * (y*y + p*p) + 3 * q73*y*y * p*p*p * (pow(y,4) + pow(p,4)) + 5 * q75*pow(y,4) * pow(p,5) * (y*y + p*p);

	return FP1 + FP2;
}

double RockingBC::I_FAb(double y, double p)
{
	double FA2 = -OMXYLOGOMXYOXY(y*p)*p / 3. * (2 * y*y * p*p + 5 * y*p - 1);
	double FA3 = OMXYLOGOMXYOXY(-y*p)*p / 3. * (1 + y*p)*(2 * y*p - 1);
	double FA4 = 4. / 3. * y*p*p+2.*(p-y);

	return FA2 + FA3 + FA4;
}

double RockingBC::J_FAb(double y, double p)
{
	double FA2 = -p*p / 6. * (OMXYLOGOMXYOXY(y*p) + YMXLOGYMX(y*p, 1.)*(3. * y*p + 7.) + J2(y*p));
	double FA3 = -p*p / 6. * (OMXYLOGOMXYOXY(-y*p) + YMXLOGYMX(-y*p, 1.)*(3. * y*p + 1.) + J2(-y*p));
	double FA4 = y*p*p*p + 5.*p*p/6.;

	return FA2 + FA3 + FA4;
}

double RockingBC::Ib_calc(double y, double r)
{
	static const double pi{ std::atan(1) * 4 };
	static const double A = -1 / pi;
	static const double B = -0.19532775;

	return A*I_FAb(y, r) + B*I_FB(y, r) + I_FP(y, r);
}

double RockingBC::Jb_calc(double y, double r)
{
	static const double pi{ std::atan(1.) * 4 };
	static const double A = -1 / pi;
	static const double B = -0.19532775;

	return A*J_FAb(y, r) + B*J_FB(y, r) + J_FP(y, r);
}

double RockingBC::pImJ_calc(double y, double r)
{
	static const double pi{ std::atan(1.) * 4 };
	static const double A = -1 / pi;
	static const double B = -0.19532775;

	return A*pImJ_FA(y, r) + B*pImJ_FB(y, r) + pImJ_FP(y, r);
}

void RockingBC::Imat_calc(const Vector& Y, const Vector& R, Matrix& Imat)
{
	for (size_t i = 0; i != Y.Size(); i++)
	{
		for (size_t j = 0; j != R.Size(); j++)
		{
			Imat(i,j) = I_calc(Y[i], R[j]);
		}
	}
	return;
}

void RockingBC::Jmat_calc(const Vector& Y, const Vector& R, Matrix& Jmat)
{
	for (size_t i = 0; i != Y.Size(); i++)
	{
		for (size_t j = 0; j != R.Size(); j++)
		{
			Jmat(i,j) = J_calc(Y[i], R[j]);
		}
	}
	return;
}

void RockingBC::Imatb_calc(const Vector& Y, const Vector& R, Matrix& Imat)
{
	for (size_t i = 0; i != Y.Size(); i++)
	{
		for (size_t j = 0; j != R.Size(); j++)
		{
			Imat(i, j) = Ib_calc(Y[i], R[j]);
		}
	}
	return;
}

void RockingBC::Jmatb_calc(const Vector& Y, const Vector& R, Matrix& Jmat)
{
	for (size_t i = 0; i != Y.Size(); i++)
	{
		for (size_t j = 0; j != R.Size(); j++)
		{
			Jmat(i, j) = Jb_calc(Y[i], R[j]);
		}
	}
	return;
}

void RockingBC::pImJmat_calc(const Vector& Y, const Vector& R, Matrix& pImJmat)
{
	for (size_t i = 0; i != Y.Size(); i++)
	{
		for (size_t j = 0; j != R.Size(); j++)
		{
			pImJmat(i, j) = pImJ_calc(Y[i], R[j]);
		}
	}
	return;
}

void RockingBC::Im1_calc(const Vector& Y, Vector& Im1)
{
	for (size_t i = 0; i != Y.Size(); i++)
	{
		Im1(i) =I_calc(Y[i], -1.0);
	}
	return;
}

void RockingBC::Jm1_calc(const Vector& Y, Vector& Jm1)
{
	for (size_t i = 0; i != Y.Size(); i++)
	{
		Jm1(i) =J_calc(Y[i], -1.0);
	}
	return;
}

void RockingBC::Im1b_calc(const Vector& Y, Vector& Im1)
{
	for (size_t i = 0; i != Y.Size(); i++)
	{
		Im1(i) = Ib_calc(Y[i], -1.0);
	}
	return;
}

void RockingBC::Jm1b_calc(const Vector& Y, Vector& Jm1)
{
	for (size_t i = 0; i != Y.Size(); i++)
	{
		Jm1(i) = Jb_calc(Y[i], -1.0);
	}
	return;
}

void RockingBC::Usgm_trapz(const Vector& Yw, Matrix& Usgm)
{
	Matrix CC (Yw.Size(),Yw.Size());
	for (int i = 0; i != Yw.Size(); i++)
	{
		if (i > 0) {
			CC(i - 1,i) += -1.0 / (Yw(i - 1) - Yw(i));
			CC(i,i) += 1.0 / (Yw(i - 1) - Yw(i));
		}
		if (i < Yw.Size() - 1) {
			CC(i,i) += 1.0 / (Yw(i) - Yw(i + 1));
			CC(i + 1,i) += -1.0 / (Yw(i) - Yw(i + 1));
		}
	}

	Matrix Imat(Yw.Size(), Yw.Size());
	Matrix Jmat(Yw.Size(), Yw.Size());
	Vector Im1(Yw.Size());
	Vector Jm1(Yw.Size());
	Imat_calc(Yw, Yw, Imat);
	Jmat_calc(Yw, Yw, Jmat);
	Im1_calc(Yw, Im1);
	Jm1_calc(Yw, Jm1);

	Matrix Us(Yw.Size(), Yw.Size());
	for (size_t i = 0; i != Yw.Size(); i++)
	{
		for (size_t k = 0; k != Yw.Size(); k++)
		{
			Us(k,i) = (Yw(i) * Imat(k,i) - Jmat(k,i)) - Im1(k) * Yw(i) + Jm1(k);
		}
	}

	Usgm = Us* CC;

	return;
}

void RockingBC::triangle_dispslope_disps(const Vector& R, const Vector& Y, Matrix& U, Matrix& dU_dR)
{
	Matrix Imat(Y.Size(), R.Size());
	Matrix Jmat(Y.Size(), R.Size());
	Vector Im1(Y.Size());
	Vector Jm1(Y.Size());
	Imat_calc(Y, R, Imat);
	Jmat_calc(Y, R, Jmat);
	Im1_calc(Y, Im1);
	Jm1_calc(Y, Jm1);

	for (size_t i = 0; i != R.Size(); i++)
	{
		for (size_t k = 0; k != Y.Size(); k++)
		{
			U(k,i) = (R(i) * Imat(k,i) - Jmat(k,i)) - Im1(k) * R(i) + Jm1(k);
			dU_dR(k,i) = Imat(k,i) - Im1(k);
		}
	}

	return;
}

void RockingBC::triangle_dispslope_disps_givenMat1(const Vector& R, const Vector& Y, const Vector& Im1, const Vector& Jm1, Matrix& U, Matrix& dU_dR)
{
	Matrix Imat(Y.Size(), R.Size());
	Matrix Jmat(Y.Size(), R.Size());
	Imat_calc(Y, R, Imat);
	Jmat_calc(Y, R, Jmat);

	for (size_t i = 0; i != R.Size(); i++)
	{
		for (size_t k = 0; k != Y.Size(); k++)
		{
			U(k, i) = (R(i) * Imat(k, i) - Jmat(k, i)) - Im1(k) * R(i) + Jm1(k);
			dU_dR(k, i) = Imat(k, i) - Im1(k);
		}
	}

	return;
}

void RockingBC::triangle_dispslope_disps_2(const Vector& R, const Vector& Y, const Vector& Im1, const Vector& Jm1, Matrix& U, Matrix& dU_dR)
{
	Matrix pImJmat(Y.Size(), R.Size());
	Matrix Imat(Y.Size(), R.Size());
	pImJmat_calc(Y, R, pImJmat);
	Imat_calc(Y, R, Imat);

	for (size_t i = 0; i != R.Size(); i++)
	{
		for (size_t k = 0; k != Y.Size(); k++)
		{
			U(k, i) = pImJmat(k, i) - Im1(k) * R(i) + Jm1(k);
			dU_dR(k, i) = Imat(k, i) - Im1(k);
		}
	}

	return;
}

void RockingBC::UNM_trapz(const Vector& R2, const Vector& R1, const Vector& Y, Matrix& U)
{
	Matrix Imata(Y.Size(), R1.Size());
	Matrix Jmata(Y.Size(), R1.Size());
	Matrix Imatb(Y.Size(), R2.Size());
	Matrix Jmatb(Y.Size(), R2.Size());
	Vector Im1(Y.Size());
	Imat_calc(Y, R1, Imata);
	Jmat_calc(Y, R1, Jmata);
	Imat_calc(Y, R2, Imatb);
	Jmat_calc(Y, R2, Jmatb);
	Im1_calc(Y, Im1);

	U = Matrix(Y.Size(), R2.Size());
	for (size_t i = 0; i != R2.Size(); i++)
	{
		for (size_t k = 0; k != Y.Size(); k++)
		{
			U(k, i) = (R2(i) * Imatb(k, i) - Jmatb(k, i)) - (R1(i) * Imata(k, i) - Jmata(k, i)) - Im1(k) * (R2(i)-R1(i));
		}
	}

	return;
}

void RockingBC::UNM_rect(const Vector& R, const Vector& Y, Matrix& U)
{
	Matrix Imat(Y.Size(), R.Size());
	Vector Im1(Y.Size());
	Imat_calc(Y, R, Imat);
	Im1_calc(Y, Im1);

	U = Matrix(Y.Size(), R.Size());
	for (size_t i = 0; i != R.Size(); i++)
	{
		for (size_t k = 0; k != Y.Size(); k++)
		{
			U(k, i) = Imat(k, i) - Im1(k);
		}
	}

	return;
}

void RockingBC::UNM_calc(const Vector& Yw, Matrix& UN, Matrix& UM)
{
	Vector R1(Yw.Size() - 1);
	for (int i = 0; i < Yw.Size() - 1; i++) {
		R1(i) = Yw(i);
	}
	Vector R2(Yw.Size() - 1);
	for (int i = 0; i < Yw.Size() - 1; i++) {
		R2(i) = Yw(i+1);
	}

	Matrix Utr{};
	Matrix Ur{};
	UNM_trapz(R2, R1, Yw, Utr);
	UNM_rect(Yw, Yw, Ur);

	Matrix Ur1(Ur.noRows(), Ur.noCols() -1);
	for (int i = 0; i < Ur.noRows(); i++) {
		for (int j = 0; j < Ur.noCols() - 1; j++) {
			Ur1(i, j) = Ur(i, j);
		}
	}
	Matrix Ur2(Ur.noRows(), Ur.noCols() - 1);
	for (int i = 0; i < Ur.noRows(); i++) {
		for (int j = 0; j < Ur.noCols() - 1; j++) {
			Ur2(i, j) = Ur(i, j+1);
		}
	}

	UN = Matrix(Yw.Size(), Yw.Size() - 1);
	UM = Matrix(Yw.Size(), Yw.Size() - 1);
	for (size_t i = 0; i != Yw.Size()-1; i++)
	{
		for (size_t k = 0; k != Yw.Size(); k++)
		{
			UN(k, i) = 6. * (Yw[i + 1] + Yw[i]) / (Yw[i + 1] - Yw[i]) / (Yw[i + 1] - Yw[i]) / (Yw[i + 1] - Yw[i]) * Utr(k, i) - 2. * (2. * Yw[i + 1] + Yw[i]) / (Yw[i + 1] - Yw[i]) / (Yw[i + 1] - Yw[i]) * Ur1(k, i) - 2. * (Yw[i + 1] + 2. * Yw[i]) / (Yw[i + 1] - Yw[i]) / (Yw[i + 1] - Yw[i]) * Ur2(k, i);
			UM(k, i) = -12. / (Yw[i + 1] - Yw[i]) / (Yw[i + 1] - Yw[i]) / (Yw[i + 1] - Yw[i]) * Utr(k, i) + 6. / (Yw[i + 1] - Yw[i]) / (Yw[i + 1] - Yw[i]) * (Ur1(k, i) + Ur2(k, i));
		}
	}
	return;
}

void RockingBC::UNMb_trapz(const Vector& R2, const Vector& R1, const Vector& Y, Matrix& U)
{
	Matrix Imata(Y.Size(), R1.Size());
	Matrix Jmata(Y.Size(), R1.Size());
	Matrix Imatb(Y.Size(), R2.Size());
	Matrix Jmatb(Y.Size(), R2.Size());
	Vector Im1(Y.Size());
	Imatb_calc(Y, R1, Imata);
	Jmatb_calc(Y, R1, Jmata);
	Imatb_calc(Y, R2, Imatb);
	Jmatb_calc(Y, R2, Jmatb);
	Im1b_calc(Y, Im1);

	U = Matrix(Y.Size(), R2.Size());
	for (size_t i = 0; i != R2.Size(); i++)
	{
		for (size_t k = 0; k != Y.Size(); k++)
		{
			U(k, i) = (R2(i) * Imatb(k, i) - Jmatb(k, i)) - (R1(i) * Imata(k, i) - Jmata(k, i)) - Im1(k) * (R2(i) - R1(i));
		}
	}

	return;
}

void RockingBC::UNMb_rect(const Vector& R, const Vector& Y, Matrix& U)
{
	Matrix Imat(Y.Size(), R.Size());
	Vector Im1(Y.Size());
	Imatb_calc(Y, R, Imat);
	Im1b_calc(Y, Im1);

	U = Matrix(Y.Size(), R.Size());
	for (size_t i = 0; i != R.Size(); i++)
	{
		for (size_t k = 0; k != Y.Size(); k++)
		{
			U(k, i) = Imat(k, i) - Im1(k);
		}
	}

	return;
}

void RockingBC::UNMb_calc(const Vector& Yw, Matrix& UN, Matrix& UM)
{
	Vector R1(Yw.Size() - 1);
	for (int i = 0; i < Yw.Size() - 1; i++) {
		R1(i) = Yw(i);
	}
	Vector R2(Yw.Size() - 1);
	for (int i = 0; i < Yw.Size() - 1; i++) {
		R2(i) = Yw(i + 1);
	}

	Matrix Utr{};
	Matrix Ur{};
	UNMb_trapz(R2, R1, Yw, Utr);
	UNMb_rect(Yw, Yw, Ur);
	//Matrix Ur1 = Ur.leftCols(Yw.Size() - 1);
	//Matrix Ur2 = Ur.rightCols(Yw.Size() - 1);

	Matrix Ur1(Ur.noRows(), Ur.noCols() - 1);
	for (int i = 0; i < Ur.noRows(); i++) {
		for (int j = 0; j < Ur.noCols() - 1; j++) {
			Ur1(i, j) = Ur(i, j);
		}
	}
	Matrix Ur2(Ur.noRows(), Ur.noCols() - 1);
	for (int i = 0; i < Ur.noRows(); i++) {
		for (int j = 0; j < Ur.noCols() - 1; j++) {
			Ur2(i, j) = Ur(i, j + 1);
		}
	}

	UN = Matrix(Yw.Size(), Yw.Size() - 1);
	UM = Matrix(Yw.Size(), Yw.Size() - 1);
	for (size_t i = 0; i != Yw.Size() - 1; i++)
	{
		for (size_t k = 0; k != Yw.Size(); k++)
		{
			UN(k, i) = 6. * (Yw[i + 1] + Yw[i]) / (Yw[i + 1] - Yw[i]) / (Yw[i + 1] - Yw[i]) / (Yw[i + 1] - Yw[i]) * Utr(k, i) - 2. * (2. * Yw[i + 1] + Yw[i]) / (Yw[i + 1] - Yw[i]) / (Yw[i + 1] - Yw[i]) * Ur1(k, i) - 2. * (Yw[i + 1] + 2. * Yw[i]) / (Yw[i + 1] - Yw[i]) / (Yw[i + 1] - Yw[i]) * Ur2(k, i);
			UM(k, i) = -12. / (Yw[i + 1] - Yw[i]) / (Yw[i + 1] - Yw[i]) / (Yw[i + 1] - Yw[i]) * Utr(k, i) + 6. / (Yw[i + 1] - Yw[i]) / (Yw[i + 1] - Yw[i]) * (Ur1(k, i) + Ur2(k, i));
		}
	}
	return;
}

// SCfuncs

void RockingBC::Dt_calc(const Vector& P, double& d, Vector& dddP)
{
	double p{ P(0) };
	double q{ P(1) };

	double a = SC_A(p);
	double b = SC_B(p);
	double c = SC_C(p);

	double da_dp = SC_DA(p);
	double db_dp = SC_DB(p);
	double dc_dp = SC_DC(p);

	d = a*pow(1.0 - pow(q, b), c);

	double dddp{ 0.0 };
	double dddq{ 0.0 };
	if (q > 0 && q < 1)
	{
		dddp = a*(std::log(1 - pow(q,b))*pow(1 - pow(q, b),c)*dc_dp - log(q)*c*pow(1 - pow(q, b),c-1)*db_dp*pow(q, b)) + pow(1 - pow(q, b),c)*da_dp;
		dddq = -a*b*c*pow(q,b-1)*pow(1 - pow(q, b),c-1);
	}
	else if (q == 1)
	{
		dddp = a*(-std::log(q)*c*pow(1 - pow(q, b),c-1)*db_dp*pow(q, b)) + pow(1 - pow(q, b),c)*da_dp;
		dddq = -a*b*c*pow(q, b-1)*pow(1 - pow(q, b),c-1);
	}
	else
	{
		throw;
	}
	dddP(0) = dddp;
	dddP(1) = dddq;

	return;

}

void RockingBC::Rt_calc(const Vector& P, double& th, Vector& dthdP)
{
	double p{ P(0) };
	double q{ P(1) };

	double d = SC_D(p);
	double e = SC_E(p);
	double f = SC_F(p);

	double dd_dp = SC_DD(p);
	double de_dp = SC_DE(p);
	double df_dp = SC_DF(p);

	th = d*pow(1 - q, e) + f;

	double dthdp{ 0.0 };
	double dthdq{ 0.0 };
	if (q > 0 && q < 1)
	{
		dthdp = pow(1 - q,e)*dd_dp + df_dp + std::log(1 - q)*d*pow(1 - q,e)*de_dp;
		dthdq = -d*e*pow(1 - q, e - 1);
	}
	else if (q == 1)
	{
		dthdp = df_dp;
		dthdq = -d*e*pow(1 - q, e - 1);
	}
	else
	{
		throw;
	}
	dthdP(0) = dthdp;
	dthdP(1) = dthdq;

	return;

}

void RockingBC::se_shear_1der(const Vector& Youter, Vector& Ut, Matrix& dUt_dYouter)
{
	double d{ 0.0 };
	double th{ 0.0 };
	static Vector dddP = Vector(2);
	static Vector dthdP = Vector(2);
	static Vector dp_dYouter = Vector(2);
	static Vector dq_dYouter = Vector(2);
	static Vector P(2);
	static Matrix dP_dYouter(2,2);

	if ((Youter(0) + Youter(1)) / 2.0 <= 0) {
		double q = 1. + (Youter(0) + Youter(1)) / 2.;
		double p = (1. + Youter(0)) / (1. + Youter(1));
		P(0) = p;
		P(1) = q;

		dp_dYouter(0) = 1. / (1. + Youter(1));
		dp_dYouter(1) = -(1. + Youter(0)) / (1. + Youter(1)) / (1. + Youter(1));
		dq_dYouter(0) = 1. / 2.;
		dq_dYouter(1) = 1. / 2.;
		dP_dYouter(0,0) = dp_dYouter(0);
		dP_dYouter(0,1) = dp_dYouter(1);
		dP_dYouter(1,0) = dq_dYouter(0);
		dP_dYouter(1,1) = dq_dYouter(1);

		Dt_calc(P,d,dddP);
		Rt_calc(P,th,dthdP);
	}
	else {

		double q = 1. - (Youter(0) + Youter(1)) / 2.;
		double p = (1. - Youter(1)) / (1. - Youter(0));
		P(0) = p;
		P(1) = q;

		dq_dYouter(0) = -1. / 2.;
		dq_dYouter(1) = -1. / 2.;
		dp_dYouter(0) = (1. - Youter(1)) / (1. - Youter(0)) / (1. - Youter(0));
		dp_dYouter(1) = -1. / (1. - Youter(0));
		dP_dYouter(0, 0) = dp_dYouter(0);
		dP_dYouter(0, 1) = dp_dYouter(1);
		dP_dYouter(1, 0) = dq_dYouter(0);
		dP_dYouter(1, 1) = dq_dYouter(1);

		Dt_calc(P, d, dddP);
		Rt_calc(P, th, dthdP);
		d = -d;
		dddP(0) = -dddP(0);
		dddP(1) = -dddP(1);
	}

	Ut(0) = d;
	Ut(1) = th;
	static Matrix dUt_dP(2, 2);
	dUt_dP(0,0) = dddP(0);
	dUt_dP(0,1) = dddP(1);
	dUt_dP(1,0) = dthdP(0);
	dUt_dP(1,1) = dthdP(1);
	dUt_dYouter = dUt_dP*dP_dYouter;

	return;

}

// region funcs

void RockingBC::Up_interval_split(const Vector& Yup, const Vector& Up, const Vector& Yw,
	VecVec& Yup_ints, VecVec& Up_ints)
{
	static std::vector<int> Yind{};
	Yind.clear();
	int iy = 0;
	for (size_t iw = 0; iw != Yw.Size(); iw++) {
		while (true) {
			if (Yup[iy] == Yw[iw]) {
			//if (std::fabs(Yup[iy]-Yw[iw])<1.0e-12) {
				Yind.push_back(iy);
				iy += 1;
				break;
			}
			iy += 1;
		}
	}
	
	Yup_ints.clear();
	Up_ints.clear();
	for (size_t i = 0; i != Yind.size() - 1; i++) {
		Vec X1{};
		for (size_t j = Yind[i]; j != Yind[i + 1] +1; j++) {
			X1.push_back(Up[j]);
		}
		Up_ints.push_back(X1);
		Vec X2{};
		for (size_t j = Yind[i]; j != Yind[i + 1] +1; j++) {
			X2.push_back(Yup[j]);
		}
		Yup_ints.push_back(X2);
	}
	return;
}

Vector RockingBC::interval_join(const VecVec& X_ints) {
	static Vec X{};
	X.clear();

	for (size_t i = 0; i != X_ints.size(); i++) {
		for (size_t j = 0; j != X_ints.at(i).size() - 1; j++) {
			X.push_back(X_ints[i][j]);
		}
	}
	X.push_back(X_ints[X_ints.size()-1][X_ints.at(X_ints.size() - 1).size() - 1]);

	static Vector XX;
	XX = Vector(X.size());

	for (size_t i = 0; i != X.size(); i++) {
		XX[i] = X[i];
	}
	return XX;
}

Matrix RockingBC::interval_join(const VecMatOS& X_ints) {
	static std::vector<int> vecints;
	vecints.clear();
	vecints.push_back(0);

	for (size_t i = 0; i != X_ints.size(); i++) {
		vecints.push_back(vecints[vecints.size() - 1] + X_ints[i].noRows()-1);
	}
	
	static Matrix res;
	res = Matrix(vecints[vecints.size() - 1] + 1, X_ints.at(0).noCols());
	for (size_t i = 0; i != X_ints.size(); i++) {
		for (size_t k = 0; k != X_ints.at(i).noRows()-1; k++) {
			for (size_t l = 0; l != X_ints.at(i).noCols(); l++) {
				res(vecints[i] + k, l) = X_ints[i](k, l);
			}
		}
	}
	
	const Matrix& mm = X_ints[X_ints.size() - 1];
	for (size_t l = 0; l != mm.noCols(); l++) {
		res(res.noRows() - 1, l) = mm(mm.noRows() - 1, l);
	}
	return res;
}

Vector RockingBC::array_join(const VecVec& X_ints) {
	Vec X{};
	for (size_t i = 0; i != X_ints.size(); i++) {
		for (size_t j = 0; j != X_ints.at(i).size(); j++) {
			X.push_back(X_ints[i][j]);
		}
	}

	Vector XX(X.size());
	for (size_t i = 0; i != X.size(); i++) {
		XX[i] = X[i];
	}
	return XX;
}

Matrix RockingBC::array_join(const VecMatOS& X_ints) {
	std::vector<int> vecints{0};
	for (size_t i = 0; i != X_ints.size(); i++) {
		vecints.push_back(vecints[vecints.size() - 1] + X_ints[i].noRows());
	}
	Matrix res = Matrix(vecints[vecints.size() - 1], X_ints.at(0).noCols());
	for (size_t i = 0; i != X_ints.size(); i++) {
		for (size_t k = 0; k != X_ints.at(i).noRows(); k++) {
			for (size_t l = 0; l != X_ints.at(i).noCols(); l++) {
				res(vecints[i] + k, l) = X_ints[i](k, l);
			}
		}
	}

	return res;
}

void RockingBC::commony(const Vec& ya, const Vec& fa, const Vec& yb, const Vec& fb, Vec& Y, Vec& FA, Vec& FB)
{
	Y.clear();
	FA.clear();
	FB.clear();

	int ia = 0;
	int ib = 0;
	while ((ia < ya.size() - 1) || (ib < yb.size() - 1))
	{
		if (ya[ia] == yb[ib]) {
			Y.push_back(ya[ia]);
			FA.push_back(fa[ia]);
			FB.push_back(fb[ib]);
			ia += 1;
			ib += 1;
		} else if (ya[ia] < yb[ib]) {
			Y.push_back(ya[ia]);
			FA.push_back(fa[ia]);
			FB.push_back(fb[ib - 1] + (ya[ia] - yb[ib - 1]) / (yb[ib] - yb[ib - 1])*(fb[ib] - fb[ib - 1]));
			ia += 1;
		}
		else {
			Y.push_back(yb[ib]);
			FB.push_back(fb[ib]);
			FA.push_back(fa[ia - 1] + (yb[ib] - ya[ia - 1]) / (ya[ia] - ya[ia - 1])*(fa[ia] - fa[ia - 1]));
			ib += 1;
		}
	}
	Y.push_back(ya[ya.size() - 1]);
	FA.push_back(fa[fa.size() - 1]);
	FB.push_back(fb[fb.size() - 1]);

	return;
}

void RockingBC::interval_interior(double wl, double wr, double ey, double dy, const Vec& up_com, const Vec& yup_com,
	const Vec& ys_com, const Vec& s_com, double beta_Dt,
	Vec& ys_new, Vec& s_new, Vecint& ys_cats, Vec& yup_new, Vec& up_new, 
	Vec& dys_new_dwl, Vec& dys_new_dwr, Vec& ds_new_dwl, Vec& ds_new_dwr, Vec& ua_pos)
{
	
	static const double pi{ std::atan(1.) * 4 };
	double DU_DS = dy / pi;
	double DAMPC = 1.0;
	if (beta_Dt >= 0) {
		DAMPC = beta_Dt / (1.0 + beta_Dt);
	}
	double eyn = ey*DU_DS;

	static Vec Y{};
	static Vec Up{};
	static Vec S{};

	commony(yup_com,up_com,ys_com,s_com,Y,Up,S);

	// Plastic displacements differences
	static Vec Upd; Upd.clear();
	double yline{};
	double kyline = (Up[Up.size()-1]-Up[0])/(Y[Y.size()-1]-Y[0]);
	for (size_t iy = 0; iy != Y.size(); iy++)
	{
		yline = Up[0]+(Y[iy]-Y[0])*kyline;
		Upd.push_back(yline - Up[iy]);
	}
	
	// Limits
	static Vec Slim; Slim.clear();
	static Vec Slimn; Slimn.clear();
	for (size_t i = 0; i != Y.size(); i++)
	{
		Slim.push_back(S[i]*DAMPC);
		Slimn.push_back(S[i]*DAMPC*DU_DS);
	}

	// Edge stress normalization
	double wln{};
	double dwln_dwl{};
	double wrn{};
	double dwrn_dwr{};
	if (wl >= Slim[0]) {
		wln = Slimn[0]+(wl-Slim[0]);
		dwln_dwl = 1.0;
	} 
	else if (wl>=ey){
		wln = wl*DU_DS;
		dwln_dwl = DU_DS;
	}
	else {
		wln = eyn + (wl - ey);
		dwln_dwl = 1.0;
	}

	if (wr >= Slim[Slim.size()-1]) {
		wrn = Slimn[Slimn.size()-1]+(wr-Slim[Slimn.size()-1]);
		dwrn_dwr = 1.0;
	}
	else if (wr >= ey) {
		wrn = wr*DU_DS;
		dwrn_dwr = DU_DS;
	}
	else {
		wrn = eyn + (wr - ey);
		dwrn_dwr = 1.0;
	}

	double kwnline=(wrn-wln)/dy;
    double dkwnline_dwl=-dwln_dwl/dy;
    double dkwnline_dwr=dwrn_dwr/dy;

	// Plastic displacements into stresses insertion
	static Vec Wn; Wn.clear();
	static Vec dWn_dwl; dWn_dwl.clear();
	static Vec dWn_dwr; dWn_dwr.clear();
	double wline{};
	for (size_t iy = 0; iy != Y.size(); iy++)
	{
		Wn.push_back(wln + (Y[iy] - Y[0])*kwnline + Upd[iy]);
		dWn_dwl.push_back(dwln_dwl + (Y[iy] - Y[0])*dkwnline_dwl);
		dWn_dwr.push_back((Y[iy] - Y[0])*dkwnline_dwr);
	}

	// Crossings
	static Vec Yf{}; Yf.clear();
	static Vec Wnf{}; Wnf.clear();
	static Vec Upf{}; Upf.clear();
	static Vec Slimnf{}; Slimnf.clear();
	static Vec dYf_dwl{}; dYf_dwl.clear();
	static Vec dWnf_dwl{}; dWnf_dwl.clear();
	static Vec dSlimnf_dwl{}; dSlimnf_dwl.clear();
	static Vec dYf_dwr{}; dYf_dwr.clear();
	static Vec dWnf_dwr{}; dWnf_dwr.clear();
	static Vec dSlimnf_dwr{}; dSlimnf_dwr.clear();
	
	double wnf1{}, wnf2{};
	bool wnf1found = false;
	bool wnf2found = false;
	double yf1{}, upf1{}, slimnf1{}, dyf1_dwl{}, dyf1_dwr{}, dwnf1_dwl{}, dwnf1_dwr{}, dslimnf1_dwl{}, dslimnf1_dwr{};
	double yf2{}, upf2{}, slimnf2{}, dyf2_dwl{}, dyf2_dwr{}, dwnf2_dwl{}, dwnf2_dwr{}, dslimnf2_dwl{}, dslimnf2_dwr{};
	for (size_t i = 0; i != Wn.size()-1; i++)
	{
        Wnf.push_back(Wn[i]);
        Yf.push_back(Y[i]);
        Upf.push_back(Up[i]);
        Slimnf.push_back(Slimn[i]);
        dWnf_dwl.push_back(dWn_dwl[i]);
        dWnf_dwr.push_back(dWn_dwr[i]);
        dYf_dwl.push_back(0.);
        dYf_dwr.push_back(0.);   
        dSlimnf_dwl.push_back(0.);
        dSlimnf_dwr.push_back(0.);
		wnf1found = false;
		wnf2found = false;

		if ((Wn[i]<Slimn[i] && Wn[i+1]>Slimn[i+1]) || (Wn[i]>Slimn[i] && Wn[i+1]<Slimn[i+1])) {
			wnf1found = true;
            yf1=Y[i]-(Y[i+1]-Y[i])*(Wn[i]-Slimn[i])/(Wn[i+1]-Slimn[i+1]-Wn[i]+Slimn[i]);
            wnf1=Wn[i]-(Wn[i]-Slimn[i])/(Wn[i+1]-Slimn[i+1]-Wn[i]+Slimn[i])*(Wn[i+1]-Wn[i]);
            upf1=Up[i]+(yf1-Y[i])/(Y[i+1]-Y[i])*(Up[i+1]-Up[i]);
            slimnf1=Slimn[i]+(yf1-Y[i])/(Y[i+1]-Y[i])*(Slimn[i+1]-Slimn[i]);
            dyf1_dwl=-(Y[i+1]-Y[i])*(dWn_dwl[i])/(Wn[i+1]-Slimn[i+1]-Wn[i]+Slimn[i])+(Y[i+1]-Y[i])*(Wn[i]-Slimn[i])/pow((Wn[i+1]-Slimn[i+1]-Wn[i]+Slimn[i]),2)*(dWn_dwl[i+1]-dWn_dwl[i]);
            dyf1_dwr=-(Y[i+1]-Y[i])*(dWn_dwr[i])/(Wn[i+1]-Slimn[i+1]-Wn[i]+Slimn[i])+(Y[i+1]-Y[i])*(Wn[i]-Slimn[i])/pow((Wn[i+1]-Slimn[i+1]-Wn[i]+Slimn[i]),2)*(dWn_dwr[i+1]-dWn_dwr[i]);
            dwnf1_dwl=dWn_dwl[i]+(dyf1_dwl)/(Y[i+1]-Y[i])*(Wn[i+1]-Wn[i])+(yf1-Y[i])/(Y[i+1]-Y[i])*(dWn_dwl[i+1]-dWn_dwl[i]);
            dwnf1_dwr=dWn_dwr[i]+(dyf1_dwr)/(Y[i+1]-Y[i])*(Wn[i+1]-Wn[i])+(yf1-Y[i])/(Y[i+1]-Y[i])*(dWn_dwr[i+1]-dWn_dwr[i]);
            dslimnf1_dwl=(dyf1_dwl)/(Y[i+1]-Y[i])*(Slimn[i+1]-Slimn[i]);
            dslimnf1_dwr=(dyf1_dwr)/(Y[i+1]-Y[i])*(Slimn[i+1]-Slimn[i]);
		}
		if ((Wn[i]<eyn && Wn[i+1]>eyn) || (Wn[i]>eyn && Wn[i+1]<eyn)) {
			wnf2found = true;
            yf2=Y[i]-(Y[i+1]-Y[i])*(Wn[i]-eyn)/(Wn[i+1]-Wn[i]);
            wnf2=eyn;
            upf2=Up[i]+(yf2-Y[i])/(Y[i+1]-Y[i])*(Up[i+1]-Up[i]);
            slimnf2=Slimn[i]+(yf2-Y[i])/(Y[i+1]-Y[i])*(Slimn[i+1]-Slimn[i]);
            dyf2_dwl=-(Y[i+1]-Y[i])*(dWn_dwl[i])/(Wn[i+1]-Wn[i])+(Y[i+1]-Y[i])*(Wn[i]-eyn)/(Wn[i+1]-Wn[i])/(Wn[i+1]-Wn[i])*(dWn_dwl[i+1]-dWn_dwl[i]);
            dyf2_dwr=-(Y[i+1]-Y[i])*(dWn_dwr[i])/(Wn[i+1]-Wn[i])+(Y[i+1]-Y[i])*(Wn[i]-eyn)/(Wn[i+1]-Wn[i])/(Wn[i+1]-Wn[i])*(dWn_dwr[i+1]-dWn_dwr[i]);
            dwnf2_dwl=0.;
            dwnf2_dwr=0.;
            dslimnf2_dwl=(dyf2_dwl)/(Y[i+1]-Y[i])*(Slimn[i+1]-Slimn[i]);
            dslimnf2_dwr=(dyf2_dwr)/(Y[i+1]-Y[i])*(Slimn[i+1]-Slimn[i]);
		}
		
		if (wnf1found && !wnf2found) {
            Wnf.push_back(wnf1);
            Yf.push_back(yf1);
            Upf.push_back(upf1);
            Slimnf.push_back(slimnf1);
            dWnf_dwl.push_back(dwnf1_dwl);
            dYf_dwl.push_back(dyf1_dwl);
            dWnf_dwr.push_back(dwnf1_dwr);
            dYf_dwr.push_back(dyf1_dwr);
            dSlimnf_dwl.push_back(dslimnf1_dwl);
            dSlimnf_dwr.push_back(dslimnf1_dwr);
		}
		if (wnf2found && !wnf1found) {
            Wnf.push_back(wnf2);
            Yf.push_back(yf2);
            Upf.push_back(upf2);
            Slimnf.push_back(slimnf2);
            dWnf_dwl.push_back(dwnf2_dwl);
            dYf_dwl.push_back(dyf2_dwl);
            dWnf_dwr.push_back(dwnf2_dwr);
            dYf_dwr.push_back(dyf2_dwr);
            dSlimnf_dwl.push_back(dslimnf2_dwl);
            dSlimnf_dwr.push_back(dslimnf2_dwr);
		}
		if (wnf1found && wnf2found) {
			if (yf1 <= yf2) {
                Wnf.push_back(wnf1);
                Yf.push_back(yf1);
                Upf.push_back(upf1);
                Slimnf.push_back(slimnf1);
                Wnf.push_back(wnf2);
                Yf.push_back(yf2);
                Upf.push_back(upf2);
                Slimnf.push_back(slimnf2);
                dWnf_dwl.push_back(dwnf1_dwl);
                dYf_dwl.push_back(dyf1_dwl);
                dWnf_dwr.push_back(dwnf1_dwr);
                dYf_dwr.push_back(dyf1_dwr);
                dWnf_dwl.push_back(dwnf2_dwl);
                dYf_dwl.push_back(dyf2_dwl);
                dWnf_dwr.push_back(dwnf2_dwr);
                dYf_dwr.push_back(dyf2_dwr);
                dSlimnf_dwl.push_back(dslimnf1_dwl);
                dSlimnf_dwr.push_back(dslimnf1_dwr);
                dSlimnf_dwl.push_back(dslimnf2_dwl);
                dSlimnf_dwr.push_back(dslimnf2_dwr);
			}
			else {
                Wnf.push_back(wnf2);
                Yf.push_back(yf2);
                Upf.push_back(upf2);
                Slimnf.push_back(slimnf2);
                Wnf.push_back(wnf1);
                Yf.push_back(yf1);
                Upf.push_back(upf1);   
                Slimnf.push_back(slimnf1);
                dWnf_dwl.push_back(dwnf2_dwl);
                dYf_dwl.push_back(dyf2_dwl);
                dWnf_dwr.push_back(dwnf2_dwr);
                dYf_dwr.push_back(dyf2_dwr);
                dWnf_dwl.push_back(dwnf1_dwl);
                dYf_dwl.push_back(dyf1_dwl);
                dWnf_dwr.push_back(dwnf1_dwr);
                dYf_dwr.push_back(dyf1_dwr);
                dSlimnf_dwl.push_back(dslimnf2_dwl);
                dSlimnf_dwr.push_back(dslimnf2_dwr);
                dSlimnf_dwl.push_back(dslimnf1_dwl);
                dSlimnf_dwr.push_back(dslimnf1_dwr);
			}
		}
	}
	
	//std::cout << Eigen::Map<Vector>(&s_com[0], s_com.size()).transpose() << std::endl;
	
    Wnf.push_back(Wn[Wn.size()-1]);
    Yf.push_back(Y[Y.size()-1]);   
    Upf.push_back(Up[Up.size()-1]);
    Slimnf.push_back(Slimn[Slimn.size()-1]);
    dWnf_dwl.push_back(dWn_dwl[dWn_dwl.size()-1]);
    dWnf_dwr.push_back(dWn_dwr[dWn_dwr.size()-1]);
    dYf_dwl.push_back(0.);
    dYf_dwr.push_back(0.);
    dSlimnf_dwl.push_back(0.);
    dSlimnf_dwr.push_back(0.);
	
	//Categorization

	static std::vector<int> intcats{};
	intcats.clear();
	for (size_t i = 0; i != Wnf.size()-1; i++)
	{
		if (0.5*(Wnf[i] + Wnf[i + 1]) > 0.5*(Slimnf[i] + Slimnf[i + 1])) {
			intcats.push_back(0);
		}
		else if (0.5*(Wnf[i] + Wnf[i + 1]) > eyn) {
			intcats.push_back(1);
		}
		else {
			intcats.push_back(2);
		}
	}

	//Separation into stresses, plastic displacements
	
	static Vec Sf_new{}; Sf_new.clear();
	static Vec dSf_new_dwl{}; dSf_new_dwl.clear();
	static Vec dSf_new_dwr{}; dSf_new_dwr.clear();
	static Vec Upf_new{}; Upf_new.clear();
	static Vec Ua_pos{}; Ua_pos.clear();

	for (size_t i = 0; i != Wnf.size(); i++) {
		if (Wnf[i] > Slimnf[i]) {
            Sf_new.push_back(Slimnf[i]/DU_DS);
			Upf_new.push_back(Upf[i]);
			Ua_pos.push_back(Wnf[i] - Slimnf[i]);
            dSf_new_dwl.push_back(dSlimnf_dwl[i]/DU_DS);
            dSf_new_dwr.push_back(dSlimnf_dwr[i]/DU_DS);
		}
		else if (Wnf[i] > eyn) {
            Sf_new.push_back(Wnf[i]/DU_DS);
			Upf_new.push_back(Upf[i]);
			Ua_pos.push_back(0.);
            dSf_new_dwl.push_back(dWnf_dwl[i]/DU_DS);
            dSf_new_dwr.push_back(dWnf_dwr[i]/DU_DS);
		}
		else {
            Sf_new.push_back(ey);
			Upf_new.push_back(Upf[i]+Wnf[i]-eyn);
			Ua_pos.push_back(0.);
            dSf_new_dwl.push_back(0.);
            dSf_new_dwr.push_back(0.);
		}
	}
	
	// Simplification
	
	ys_new.clear();
	s_new.clear();
	ys_cats.clear();
	dys_new_dwl.clear();
	dys_new_dwr.clear();
	ds_new_dwl.clear();
	ds_new_dwr.clear();
	yup_new.clear();
	up_new.clear();
	ua_pos.clear();

	ys_new.push_back(Yf[0]);
	s_new.push_back(Sf_new[0]);	
	
	dys_new_dwl.push_back(dYf_dwl[0]);
	ds_new_dwl.push_back(dSf_new_dwl[0]);	
	dys_new_dwr.push_back(dYf_dwr[0]);
	ds_new_dwr.push_back(dSf_new_dwr[0]);	
	ys_cats.push_back(intcats[0]);
	ua_pos.push_back(Ua_pos[0]);
	
	for (size_t i = 1; i != Yf.size()-1; i++) {
		if (intcats[i-1]==2 && intcats[i]==2) {
			continue;
		}
		ys_new.push_back(Yf[i]);
		s_new.push_back(Sf_new[i]);
		ys_cats.push_back(intcats[i]);
		dys_new_dwl.push_back(dYf_dwl[i]);
		ds_new_dwl.push_back(dSf_new_dwl[i]);	
		dys_new_dwr.push_back(dYf_dwr[i]);
		ds_new_dwr.push_back(dSf_new_dwr[i]);	
		ua_pos.push_back(Ua_pos[i]);
	}
	ys_new.push_back(Yf[Yf.size()-1]);
	s_new.push_back(Sf_new[Sf_new.size()-1]);
	dys_new_dwl.push_back(dYf_dwl[dYf_dwl.size()-1]);
	ds_new_dwl.push_back(dSf_new_dwl[dSf_new_dwl.size()-1]);
	dys_new_dwr.push_back(dYf_dwr[dYf_dwr.size()-1]);
	ds_new_dwr.push_back(dSf_new_dwr[dSf_new_dwr.size()-1]);	
	ua_pos.push_back(Ua_pos[Ua_pos.size()-1]);
	
	yup_new.push_back(Yf[0]);
	up_new.push_back(Upf_new[0]);	
	for (size_t i = 1; i != Yf.size()-1; i++) {
		if (intcats[i - 1] == 2 && intcats[i] == 2) {
			continue;
		}
		yup_new.push_back(Yf[i]);
		up_new.push_back(Upf_new[i]);
	}
	yup_new.push_back(Yf[Yf.size()-1]);
	up_new.push_back(Upf_new[Upf_new.size()-1]);

	return;
}

void RockingBC::interval_dists(const Vector& Yw, const Vector& W, const VecVec& Yupi_com, const VecVec& Upi_com, const VecVec& Ysi_com, const VecVec& Si_com, double ey, double beta_Dt,
	VecVec& Ysi, VecVec& Si, VecVec& Yupi_new, VecVec& Upi_new,
	VecVecint& Ys_cats, Vector& Nints, Vector& Mints, Matrix& dNints_dW, Matrix& dMints_dW, VecVec& Ua_pos, VecMatOS& dYsi_dW, VecMatOS& dSi_dW)
{
	
	VecVec dys_dwl_list( W.Size() - 1, std::vector<double>{} );
	VecVec ds_dwl_list( W.Size() - 1, std::vector<double>{} );
	VecVec dys_dwr_list( W.Size() - 1, std::vector<double>{} );
	VecVec ds_dwr_list( W.Size() - 1, std::vector<double>{} );
	
	for (size_t i = 0; i != W.Size()-1; i++) {

		interval_interior(W[i], W[i + 1], ey, Yw[i + 1] - Yw[i], Upi_com[i], Yupi_com[i],Ysi_com[i], Si_com[i], beta_Dt,
			Ysi[i], Si[i], Ys_cats[i], Yupi_new[i], Upi_new[i],
			dys_dwl_list[i], dys_dwr_list[i], ds_dwl_list[i], ds_dwr_list[i],
			Ua_pos[i]);
	}

	static Vector dNdW{}, dMdW{};
	
	for (size_t i = 0; i != W.Size() - 1; i++) {
		
		Vec dwl_dW(W.Size()); dwl_dW[i] = 1.0;
		Vec dwr_dW(W.Size()); dwr_dW[i+1] = 1.0;
		Matrix dys_dW = Matrix(dys_dwl_list[i].size(), W.Size());
		Matrix ds_dW = Matrix(ds_dwl_list[i].size(), W.Size());
		for (size_t l = 0; l != W.Size(); l++) {
			for (size_t k = 0; k != dys_dwl_list[i].size(); k++) {
				dys_dW(k, l) += dys_dwl_list[i][k] * dwl_dW[l];
				dys_dW(k, l) += dys_dwr_list[i][k] * dwr_dW[l];
				ds_dW(k, l) += ds_dwl_list[i][k] * dwl_dW[l];
				ds_dW(k, l) += ds_dwr_list[i][k] * dwr_dW[l];
			}
		}
		dYsi_dW[i] = dys_dW;
		dSi_dW[i] = ds_dW;
		NM_calc_int(Ysi[i], dys_dW, Si[i], ds_dW, Nints[i], Mints[i], dNdW, dMdW);
		for (size_t j = 0; j != W.Size(); j++) {
			dNints_dW(i,j) = dNdW(j);
			dMints_dW(i,j) = dMdW(j);
		}
	}
	return;

}

void RockingBC::NM_calc_int(const Vec& Ys, const Matrix& dYs_dW, const Vec& S, const Matrix& dS_dW, double& N, double& M, Vector& dN_dW, Vector& dM_dW)
{
	N = 0;
	M = 0;
	dN_dW = Vector(dYs_dW.noCols());
	dM_dW = Vector(dS_dW.noCols());

	for (size_t i = 0; i != Ys.size() - 1; i++)
	{

		N += (Ys[i + 1] - Ys[i])*(S[i] + S[i + 1]) / 2.;
		M += (Ys[i + 1] - Ys[i])*(2 * S[i] * Ys[i] + S[i] * Ys[i + 1] + S[i + 1] * Ys[i] + 2 * S[i + 1] * Ys[i + 1]) / 6.;

		for (size_t j= 0 ; j != dN_dW.Size(); j++) {
			dN_dW(j) += (-S[i] / 2. - S[i + 1] / 2.) * dYs_dW(i,j) + (S[i] / 2. + S[i + 1] / 2.) * dYs_dW(i + 1,j) + (Ys[i + 1] / 2. - Ys[i] / 2.) * dS_dW(i,j) + (Ys[i + 1] / 2. - Ys[i] / 2.) * dS_dW(i + 1,j);
			dM_dW(j) += (-(S[i] * Ys[i]) / 3. - (S[i] * Ys[i + 1]) / 6. - (S[i + 1] * Ys[i]) / 6. - (S[i + 1] * Ys[i + 1]) / 3. - ((2. * S[i] + S[i + 1]) * (Ys[i] - Ys[i + 1])) / 6.) * dYs_dW(i,j) +
				((S[i] * Ys[i]) / 3. + (S[i] * Ys[i + 1]) / 6. + (S[i + 1] * Ys[i]) / 6. + (S[i + 1] * Ys[i + 1]) / 3. - ((S[i] + 2. * S[i + 1]) * (Ys[i] - Ys[i + 1])) / 6.) * dYs_dW(i + 1,j) +
				(-((Ys[i] - Ys[i + 1]) * (2. * Ys[i] + Ys[i + 1])) / 6.) * dS_dW(i,j) + (-((Ys[i] - Ys[i + 1]) * (Ys[i] + 2. * Ys[i + 1])) / 6.) * dS_dW(i + 1,j);
		}

	}
	return;
}

void RockingBC::critpoints(const Vec& y, const Vec& s, int rinit, int rend, Vecint& cp)
{
	cp.clear();

	for (size_t i = rinit + 1; i != rend; i++) {
		// Slope change
		if ((s[i] - s[i - 1])*(s[i + 1] - s[i]) <= 0 && (s[i] - s[i - 1] != 0 || s[i + 1] - s[i] != 0))
		{
			cp.push_back(i);
			continue;
		}
		// Curvature change
		//if (i != rend - 1) {
		//	if (((s[i + 1] - s[i])*(y[i] - y[i - 1])*(y[i + 2] - y[i + 1]) - (s[i] - s[i - 1])*(y[i + 1] - y[i])*(y[i + 2] - y[i + 1]))*(
		//		(s[i + 2] - s[i + 1])*(y[i + 1] - y[i])*(y[i] - y[i - 1]) - (s[i + 1] - s[i])*(y[i] - y[i - 1])*(y[i + 2] - y[i + 1])) < 0) {
		//		cp.push_back(i);
		//		cp.push_back(i + 1);
		//	}
		//}
	}
	return;
}

void RockingBC::int_bilin(const Vecint& ys_cats, const Vec& ys, const Vec& s, const Vec& yup, const Vec& up, const Vec& ua_pos, double ey,
	Vec& ys_new, Vec& s_new, Vec& yup_new, Vec& up_new)
{

	if ((ys.size() != yup.size()) || (ys_cats.size() != ys.size() - 1) || (ys.size() != ua_pos.size())) {
		ys_new = ys;
		s_new = s;
		yup_new = yup;
		up_new = up;
		return;
	}

	static const double pi{ std::atan(1.) * 4 };
	double DU_DS = (ys[ys.size()-1] - ys[0]) / pi;

	static Vec w{}; w.clear();
	static Vec sm{}; sm.clear();
	for (size_t i = 0; i != s.size(); i++) {
		w.push_back(s[i] * DU_DS + ua_pos[i]);
		sm.push_back(s[i] * DU_DS);
	}

	static Vecint regi{}; regi.clear();
	regi.push_back(0);
	for (size_t i = 0; i != ys_cats.size()-1; i++) {
		if (ys_cats[i + 1] != ys_cats[i]) {
			regi.push_back(i + 1);
		}
	}
	regi.push_back(ys.size() - 1);

	static VecVecint regs{}; regs.clear();
	static Vecint regs_cats{}; regs_cats.clear();
	for (size_t i = 0; i != regi.size()-1; i++) {
		if (ys_cats[regi[i]] != 2) {
			Vecint v = { regi[i], regi[i + 1] };
			regs.push_back(v);
			regs_cats.push_back(ys_cats[regi[i]]);
		}
	}

	static VecVecint regs2{}; regs2.clear();
	static Vecint regs2_cats{}; regs2_cats.clear();
	static Vecint m{}; m.clear();
	static Vecint cp{}; cp.clear();
	for (size_t ir = 0; ir != regs.size(); ir++) {
		Vecint r = regs[ir];
		m.clear();
		m.push_back(r[0]);
		critpoints(ys, s, r[0], r[1], cp);
		m.insert(m.end(), cp.begin(), cp.end());
		critpoints(yup, up, r[0], r[1], cp);
		m.insert(m.end(), cp.begin(), cp.end());
		m.push_back(r[1]);
		sort(m.begin(), m.end());
		m.erase(unique(m.begin(), m.end()), m.end());
		for (size_t j = 0; j != m.size() - 1; j++) {
			Vecint v = { m[j],m[j + 1] };
			regs2.push_back(v);
			regs2_cats.push_back(regs_cats[ir]);
		}
	}
	static VecVecint i_s_bl{}; i_s_bl.clear();
	static VecVecint i_up_bl{}; i_up_bl.clear();
	static VecVec s_bl{}; s_bl.clear();
	static VecVec ys_bl{}; ys_bl.clear();
	static VecVec up_bl{}; up_bl.clear();
	static VecVec yup_bl{}; yup_bl.clear();

	static Vec ysp{};
	static Vec smp{};
	static Vec sp{};
	static Vec wp{};
	static Vec ss{};
	static Vec ys_try{}, sm_try{}, yw_try{}, w_try{}, s_try{};
	static Vecint v{};

	for (size_t ir = 0; ir != regs2.size(); ir++) {
		Vecint r = regs2[ir];
		if (regs2_cats[ir] == 0) {
			ysp = Vec(ys.begin() + r[0], ys.begin() + r[1] + 1);
			smp = Vec(sm.begin() + r[0], sm.begin() + r[1] + 1);
			wp = Vec(w.begin() + r[0], w.begin() + r[1] + 1);
			ys_try.clear(); sm_try.clear(); yw_try.clear(); w_try.clear();
			bool suc = bilin_two(ysp, smp, ysp, wp, ys_try, sm_try, yw_try, w_try);
			if (suc) {

				if (yw_try.size() == 2) {
					v = { r[0],r[1] };
					i_s_bl.push_back(v);
					i_up_bl.push_back(v);
					ys_bl.push_back(ys_try);
					ss.clear();
					for (size_t i = 0; i != sm_try.size(); i++) {
						ss.push_back(sm_try[i] / DU_DS);
					}
					s_bl.push_back(ss);

					ss = { yup[r[0]],yup[r[r.size() - 1]] };
					yup_bl.push_back(ss);
					ss = { up[r[0]],up[r[r.size() - 1]] };
					up_bl.push_back(ss);
				}
				else {
					double yupn0 = yw_try[1];
					double upn0 = (up[r[0]] + (yupn0 - yup[r[0]]) / (yup[r[1]] - yup[r[0]])*(up[r[1]] - up[r[0]])) + (w_try[0] + (yupn0 - yw_try[0]) / (yw_try[yw_try.size() - 1] - yw_try[0])*(w_try[w_try.size() - 1] - w_try[0])) - (w_try[1]);
					
					if (sm_try[1] <= 0 && sm_try[1] / DU_DS >= ey && upn0 <= 0) {
						v = { r[0],r[1] };
						i_s_bl.push_back(v);
						i_up_bl.push_back(v);
						ys_bl.push_back(ys_try);
						ss.clear();
						for (size_t i = 0; i != sm_try.size(); i++) {
							ss.push_back(sm_try[i] / DU_DS);
						}
						s_bl.push_back(ss);

						ss = { yup[r[0]],yupn0,yup[r[1]] };
						yup_bl.push_back(ss);
						ss = { up[r[0]],upn0,up[r[1]] };
						up_bl.push_back(ss);
					}
				}
			}
		}
		else {
			ysp = Vec(ys.begin() + r[0], ys.begin() + r[1] + 1);
			sp = Vec(s.begin() + r[0], s.begin() + r[1] + 1);
			ys_try.clear(); s_try.clear();
			bool suc = bilin_one(ysp, sp, ys_try, s_try);
			if (suc) {
				if (ys_try.size() == 2) {

					v= { r[0],r[1] };
					i_s_bl.push_back(v);
					i_up_bl.push_back(v);
					ys_bl.push_back(ys_try);
					s_bl.push_back(s_try);

					ss = { yup[r[0]],yup[r[1]] };
					yup_bl.push_back(ss);
					ss = { up[r[0]],up[r[1]] };
					up_bl.push_back(ss);
				}
				else {
					double yupn0 = ys_try[1];
					double upn0 = (up[r[0]] + (yupn0 - yup[r[0]]) / (yup[r[r.size() - 1]] - yup[r[0]])*(up[r[r.size() - 1]] - up[r[0]])) - (s_try[1] - (s_try[0] + (yupn0 - ys_try[0]) / (ys_try[ys_try.size() - 1] - ys_try[0])*(s_try[s_try.size() - 1] - s_try[0])))*DU_DS;
					
					if (s_try[1] <= 0 && s_try[1] >= ey && upn0 <= 0) {
						v = { r[0],r[1] };
						i_s_bl.push_back(v);
						i_up_bl.push_back(v);
						ys_bl.push_back(ys_try);
						s_bl.push_back(s_try);

						ss = { yup[r[0]],yupn0,yup[r[1]] };
						yup_bl.push_back(ss);
						ss = { up[r[0]],upn0,up[r[1]] };
						up_bl.push_back(ss);
					}
				}
			}

		}
	}

	ys_new.clear();
	s_new.clear();
	int ir = 0;
	int i = 0;
	while (i < ys.size()) {
		if (ir < i_s_bl.size() && (i == i_s_bl[ir][0] || i == i_s_bl[ir][0] + 1)) {
			if (i == i_s_bl[ir][0]) {
				for (size_t k = 0; k != ys_bl[ir].size(); k++) {
					ys_new.push_back(ys_bl[ir][k]);
					s_new.push_back(s_bl[ir][k]);
				}
			}
			else {
				for (size_t k = 1; k != ys_bl[ir].size(); k++) {
					ys_new.push_back(ys_bl[ir][k]);
					s_new.push_back(s_bl[ir][k]);
				}
			}
			i = i_s_bl[ir][1] + 1;
			ir += 1;
		}
		else {
			ys_new.push_back(ys[i]);
			s_new.push_back(s[i]);
			i += 1;
		}
	}

	yup_new.clear();
	up_new.clear();
	ir = 0;
	i = 0;
	while (i < yup.size()) {
		if (ir < i_up_bl.size() && (i == i_up_bl[ir][0] || i == i_up_bl[ir][0] + 1)) {
			if (i == i_up_bl[ir][0]) {
				for (size_t k = 0; k != yup_bl[ir].size(); k++) {
					yup_new.push_back(yup_bl[ir][k]);
					up_new.push_back(up_bl[ir][k]);
				}
			}
			else {
				for (size_t k = 1; k != yup_bl[ir].size(); k++) {
					yup_new.push_back(yup_bl[ir][k]);
					up_new.push_back(up_bl[ir][k]);
				}
			}
			i = i_up_bl[ir][1] + 1;
			ir += 1;
		}
		else {
			yup_new.push_back(yup[i]);
			up_new.push_back(up[i]);
			i += 1;
		}
	}
}

void RockingBC::Up_interval_split_K(const Vector& Yup, const Vector& Up, const Vector& Kup, const Vector& Yw,
	VecVecOS& Yup_ints, VecVecOS& Up_ints, VecVecOS& Kup_ints)
{
	static std::vector<int> Yind{};
	Yind.clear();

	int iy = 0;
	for (size_t iw = 0; iw != Yw.Size(); iw++) {
		while (true) {
			if (Yup[iy] == Yw[iw]) {
				//if (std::fabs(Yup[iy]-Yw[iw])<1.0e-15) {
				Yind.push_back(iy);
				iy += 1;
				break;
			}
			iy += 1;
		}
	}

	Yup_ints.clear();
	Up_ints.clear();
	Kup_ints.clear();
	for (size_t i = 0; i != Yind.size() - 1; i++) {
		//Up_ints.push_back(Up.segment(Yind[i], Yind[i + 1] - Yind[i] + 1));
		//Yup_ints.push_back(Yup.segment(Yind[i], Yind[i + 1] - Yind[i] + 1));
		//Kup_ints.push_back(Kup.segment(Yind[i], Yind[i + 1] - Yind[i]));

		Vector upint(Yind[i + 1] - Yind[i] + 1);
		Vector yupint(Yind[i + 1] - Yind[i] + 1);
		for (size_t j = 0; j != Yind[i + 1] - Yind[i] + 1; j++) {
			upint(j) = Up(Yind[i] + j);
			yupint(j) = Yup(Yind[i] + j);
		}
		Vector kupint(Yind[i + 1] - Yind[i]);
		for (size_t j = 0; j != Yind[i + 1] - Yind[i]; j++) {
			kupint(j) = Kup(Yind[i] + j);
		}
		Up_ints.push_back(upint);
		Yup_ints.push_back(yupint);
		Kup_ints.push_back(kupint);
	}

	return;
}

void RockingBC::commony_K(const Vector& ya, const Vector& fa, const Vector& ka, const Vector& yb, const Vector& fb, const Vector& kb, Vec& Y, Vec& FA, Vec& FB, Vec& KA, Vec& KB)
{
	Y.clear();
	FA.clear();
	FB.clear();
	KA.clear();
	KB.clear();

	int ia = 0;
	int ib = 0;
	while ((ia < ya.Size() - 1) || (ib < yb.Size() - 1))
	{
		if (ya[ia] == yb[ib]) {
			Y.push_back(ya[ia]);
			FA.push_back(fa[ia]);
			FB.push_back(fb[ib]);
			KA.push_back(ka[ia]);
			KB.push_back(kb[ib]);
			ia += 1;
			ib += 1;
		}
		else if (ya[ia] < yb[ib]) {
			Y.push_back(ya[ia]);
			FA.push_back(fa[ia]);
			FB.push_back(fb[ib - 1] + (ya[ia] - yb[ib - 1]) / (yb[ib] - yb[ib - 1]) * (fb[ib] - fb[ib - 1]));
			KA.push_back(ka[ia]);
			KB.push_back(kb[ib - 1]);
			ia += 1;
		}
		else {
			Y.push_back(yb[ib]);
			FB.push_back(fb[ib]);
			FA.push_back(fa[ia - 1] + (yb[ib] - ya[ia - 1]) / (ya[ia] - ya[ia - 1]) * (fa[ia] - fa[ia - 1]));
			KB.push_back(kb[ib]);
			KA.push_back(ka[ia - 1]);
			ib += 1;
		}
	}
	Y.push_back(ya[ya.Size() - 1]);
	FA.push_back(fa[fa.Size() - 1]);
	FB.push_back(fb[fb.Size() - 1]);

	return;
}

void RockingBC::interval_interior_K(double wl, double wr, double ey, double dy, const Vector& up_com, const Vector& yup_com, const Vector& kup_com,
	const Vector& ys_com, const Vector& s_com, const Vector& ks_com, double beta_Dt,
	Vec& ys_new, Vec& s_new, Vec& ks_new, Vecint& ys_cats, Vec& yup_new, Vec& up_new, Vec& kup_new,
	Vec& dys_new_dwl, Vec& dys_new_dwr, Vec& ds_new_dwl, Vec& ds_new_dwr, Vec& dks_new_dwl, Vec& dks_new_dwr,
	Vec& ydks, Vec& dks, Vec& dydks_dwl, Vec& dydks_dwr, Vec& ddks_dwl, Vec& ddks_dwr,
	Vec& ds, Vec& dds_dwl, Vec& dds_dwr)
{
	static const double pi{ std::atan(1.) * 4 };
	double DU_DS = dy / pi;
	double DAMPC = beta_Dt / (1.0 + beta_Dt);
	double eyn = ey * DU_DS;

	static Vec Y{};
	static Vec Up{};
	static Vec S{};
	static Vec KUp{};
	static Vec KS{};

	commony_K(yup_com, up_com, kup_com, ys_com, s_com, ks_com, Y, Up, S, KUp, KS);

	// Plastic displacements differences
	static Vec Upd; Upd.clear();
	static Vec KUpd; KUpd.clear();
	double yline{};
	double kyline = (Up[Up.size() - 1] - Up[0]) / (Y[Y.size() - 1] - Y[0]);
	for (size_t iy = 0; iy != Y.size(); iy++)
	{
		yline = Up[0] + (Y[iy] - Y[0]) * kyline;
		Upd.push_back(yline - Up[iy]);
	}
	for (size_t iy = 0; iy != Y.size() - 1; iy++)
	{
		KUpd.push_back(kyline - KUp[iy]);
	}

	// Limits
	static Vec Slim; Slim.clear();
	static Vec Slimn; Slimn.clear();
	static Vec KSlim; KSlim.clear();
	for (size_t i = 0; i != Y.size(); i++)
	{
		Slim.push_back(S[i] * DAMPC);
		Slimn.push_back(S[i] * DAMPC * DU_DS);
	}
	for (size_t i = 0; i != Y.size() - 1; i++)
	{
		KSlim.push_back(KS[i] * DAMPC);
	}

	// Edge stress normalization
	double wln{};
	double dwln_dwl{};
	double wrn{};
	double dwrn_dwr{};
	if (wl >= Slim[0]) {
		wln = Slimn[0] + (wl - Slim[0]);
		dwln_dwl = 1.0;
	}
	else if (wl >= ey) {
		wln = wl * DU_DS;
		dwln_dwl = DU_DS;
	}
	else {
		wln = eyn + (wl - ey);
		dwln_dwl = 1.0;
	}

	if (wr >= Slim[Slim.size() - 1]) {
		wrn = Slimn[Slimn.size() - 1] + (wr - Slim[Slimn.size() - 1]);
		dwrn_dwr = 1.0;
	}
	else if (wr >= ey) {
		wrn = wr * DU_DS;
		dwrn_dwr = DU_DS;
	}
	else {
		wrn = eyn + (wr - ey);
		dwrn_dwr = 1.0;
	}

	double kwnline = (wrn - wln) / dy;
	double dkwnline_dwl = -dwln_dwl / dy;
	double dkwnline_dwr = dwrn_dwr / dy;

	// Plastic displacements into stresses insertion
	static Vec Wn; Wn.clear();
	static Vec dWn_dwl; dWn_dwl.clear();
	static Vec dWn_dwr; dWn_dwr.clear();
	static Vec KWn; KWn.clear();
	static Vec dKWn_dwl; dKWn_dwl.clear();
	static Vec dKWn_dwr; dKWn_dwr.clear();
	double wline{};
	for (size_t iy = 0; iy != Y.size(); iy++)
	{
		Wn.push_back(wln + (Y[iy] - Y[0]) * kwnline + Upd[iy]);
		dWn_dwl.push_back(dwln_dwl + (Y[iy] - Y[0]) * dkwnline_dwl);
		dWn_dwr.push_back((Y[iy] - Y[0]) * dkwnline_dwr);
	}
	for (size_t iy = 0; iy != Y.size() - 1; iy++)
	{
		KWn.push_back(kwnline + KUpd[iy]);
		dKWn_dwl.push_back(dkwnline_dwl);
		dKWn_dwr.push_back(dkwnline_dwr);
	}

	// Crossings
	static Vec Yf{}; Yf.clear();
	static Vec Wnf{}; Wnf.clear();
	static Vec Upf{}; Upf.clear();
	static Vec Slimnf{}; Slimnf.clear();
	static Vec dYf_dwl{}; dYf_dwl.clear();
	static Vec dWnf_dwl{}; dWnf_dwl.clear();
	static Vec dSlimnf_dwl{}; dSlimnf_dwl.clear();
	static Vec dYf_dwr{}; dYf_dwr.clear();
	static Vec dWnf_dwr{}; dWnf_dwr.clear();
	static Vec dSlimnf_dwr{}; dSlimnf_dwr.clear();
	static Vec KWnf{}; KWnf.clear();
	static Vec KUpf{}; KUpf.clear();
	static Vec KSlimf{}; KSlimf.clear();
	static Vec dKWnf_dwl{}; dKWnf_dwl.clear();
	static Vec dKWnf_dwr{}; dKWnf_dwr.clear();

	double wnf1{}, wnf2{};
	bool wnf1found = false;
	bool wnf2found = false;
	double yf1{}, upf1{}, slimnf1{}, dyf1_dwl{}, dyf1_dwr{}, dwnf1_dwl{}, dwnf1_dwr{}, dslimnf1_dwl{}, dslimnf1_dwr{};
	double yf2{}, upf2{}, slimnf2{}, dyf2_dwl{}, dyf2_dwr{}, dwnf2_dwl{}, dwnf2_dwr{}, dslimnf2_dwl{}, dslimnf2_dwr{};
	for (size_t i = 0; i != Wn.size() - 1; i++)
	{
		Wnf.push_back(Wn[i]);
		Yf.push_back(Y[i]);
		Upf.push_back(Up[i]);
		Slimnf.push_back(Slimn[i]);
		dWnf_dwl.push_back(dWn_dwl[i]);
		dWnf_dwr.push_back(dWn_dwr[i]);
		dYf_dwl.push_back(0.);
		dYf_dwr.push_back(0.);
		dSlimnf_dwl.push_back(0.);
		dSlimnf_dwr.push_back(0.);
		KWnf.push_back(KWn[i]);
		dKWnf_dwl.push_back(dKWn_dwl[i]);
		dKWnf_dwr.push_back(dKWn_dwr[i]);
		KUpf.push_back(KUp[i]);
		KSlimf.push_back(KSlim[i]);
		wnf1found = false;
		wnf2found = false;

		if ((Wn[i]<Slimn[i] && Wn[i + 1]>Slimn[i + 1]) || (Wn[i] > Slimn[i] && Wn[i + 1] < Slimn[i + 1])) {
			wnf1found = true;
			yf1 = Y[i] - (Y[i + 1] - Y[i]) * (Wn[i] - Slimn[i]) / (Wn[i + 1] - Slimn[i + 1] - Wn[i] + Slimn[i]);
			wnf1 = Wn[i] - (Wn[i] - Slimn[i]) / (Wn[i + 1] - Slimn[i + 1] - Wn[i] + Slimn[i]) * (Wn[i + 1] - Wn[i]);
			upf1 = Up[i] + (yf1 - Y[i]) / (Y[i + 1] - Y[i]) * (Up[i + 1] - Up[i]);
			slimnf1 = Slimn[i] + (yf1 - Y[i]) / (Y[i + 1] - Y[i]) * (Slimn[i + 1] - Slimn[i]);
			dyf1_dwl = -(Y[i + 1] - Y[i]) * (dWn_dwl[i]) / (Wn[i + 1] - Slimn[i + 1] - Wn[i] + Slimn[i]) + (Y[i + 1] - Y[i]) * (Wn[i] - Slimn[i]) / pow((Wn[i + 1] - Slimn[i + 1] - Wn[i] + Slimn[i]), 2) * (dWn_dwl[i + 1] - dWn_dwl[i]);
			dyf1_dwr = -(Y[i + 1] - Y[i]) * (dWn_dwr[i]) / (Wn[i + 1] - Slimn[i + 1] - Wn[i] + Slimn[i]) + (Y[i + 1] - Y[i]) * (Wn[i] - Slimn[i]) / pow((Wn[i + 1] - Slimn[i + 1] - Wn[i] + Slimn[i]), 2) * (dWn_dwr[i + 1] - dWn_dwr[i]);
			dwnf1_dwl = dWn_dwl[i] + (dyf1_dwl) / (Y[i + 1] - Y[i]) * (Wn[i + 1] - Wn[i]) + (yf1 - Y[i]) / (Y[i + 1] - Y[i]) * (dWn_dwl[i + 1] - dWn_dwl[i]);
			dwnf1_dwr = dWn_dwr[i] + (dyf1_dwr) / (Y[i + 1] - Y[i]) * (Wn[i + 1] - Wn[i]) + (yf1 - Y[i]) / (Y[i + 1] - Y[i]) * (dWn_dwr[i + 1] - dWn_dwr[i]);
			dslimnf1_dwl = (dyf1_dwl) / (Y[i + 1] - Y[i]) * (Slimn[i + 1] - Slimn[i]);
			dslimnf1_dwr = (dyf1_dwr) / (Y[i + 1] - Y[i]) * (Slimn[i + 1] - Slimn[i]);
		}
		if ((Wn[i]<eyn && Wn[i + 1]>eyn) || (Wn[i] > eyn&& Wn[i + 1] < eyn)) {
			wnf2found = true;
			yf2 = Y[i] - (Y[i + 1] - Y[i]) * (Wn[i] - eyn) / (Wn[i + 1] - Wn[i]);
			wnf2 = eyn;
			upf2 = Up[i] + (yf2 - Y[i]) / (Y[i + 1] - Y[i]) * (Up[i + 1] - Up[i]);
			slimnf2 = Slimn[i] + (yf2 - Y[i]) / (Y[i + 1] - Y[i]) * (Slimn[i + 1] - Slimn[i]);
			dyf2_dwl = -(Y[i + 1] - Y[i]) * (dWn_dwl[i]) / (Wn[i + 1] - Wn[i]) + (Y[i + 1] - Y[i]) * (Wn[i] - eyn) / (Wn[i + 1] - Wn[i]) / (Wn[i + 1] - Wn[i]) * (dWn_dwl[i + 1] - dWn_dwl[i]);
			dyf2_dwr = -(Y[i + 1] - Y[i]) * (dWn_dwr[i]) / (Wn[i + 1] - Wn[i]) + (Y[i + 1] - Y[i]) * (Wn[i] - eyn) / (Wn[i + 1] - Wn[i]) / (Wn[i + 1] - Wn[i]) * (dWn_dwr[i + 1] - dWn_dwr[i]);
			dwnf2_dwl = 0.;
			dwnf2_dwr = 0.;
			dslimnf2_dwl = (dyf2_dwl) / (Y[i + 1] - Y[i]) * (Slimn[i + 1] - Slimn[i]);
			dslimnf2_dwr = (dyf2_dwr) / (Y[i + 1] - Y[i]) * (Slimn[i + 1] - Slimn[i]);
		}

		if (wnf1found && !wnf2found) {
			Wnf.push_back(wnf1);
			Yf.push_back(yf1);
			Upf.push_back(upf1);
			Slimnf.push_back(slimnf1);
			dWnf_dwl.push_back(dwnf1_dwl);
			dYf_dwl.push_back(dyf1_dwl);
			dWnf_dwr.push_back(dwnf1_dwr);
			dYf_dwr.push_back(dyf1_dwr);
			dSlimnf_dwl.push_back(dslimnf1_dwl);
			dSlimnf_dwr.push_back(dslimnf1_dwr);
			KWnf.push_back(KWn[i]);
			dKWnf_dwl.push_back(dKWn_dwl[i]);
			dKWnf_dwr.push_back(dKWn_dwr[i]);
			KUpf.push_back(KUp[i]);
			KSlimf.push_back(KSlim[i]);
		}
		if (wnf2found && !wnf1found) {
			Wnf.push_back(wnf2);
			Yf.push_back(yf2);
			Upf.push_back(upf2);
			Slimnf.push_back(slimnf2);
			dWnf_dwl.push_back(dwnf2_dwl);
			dYf_dwl.push_back(dyf2_dwl);
			dWnf_dwr.push_back(dwnf2_dwr);
			dYf_dwr.push_back(dyf2_dwr);
			dSlimnf_dwl.push_back(dslimnf2_dwl);
			dSlimnf_dwr.push_back(dslimnf2_dwr);
			KWnf.push_back(KWn[i]);
			dKWnf_dwl.push_back(dKWn_dwl[i]);
			dKWnf_dwr.push_back(dKWn_dwr[i]);
			KUpf.push_back(KUp[i]);
			KSlimf.push_back(KSlim[i]);
		}
		if (wnf1found && wnf2found) {
			if (yf1 <= yf2) {
				Wnf.push_back(wnf1);
				Yf.push_back(yf1);
				Upf.push_back(upf1);
				Slimnf.push_back(slimnf1);
				Wnf.push_back(wnf2);
				Yf.push_back(yf2);
				Upf.push_back(upf2);
				Slimnf.push_back(slimnf2);
				dWnf_dwl.push_back(dwnf1_dwl);
				dYf_dwl.push_back(dyf1_dwl);
				dWnf_dwr.push_back(dwnf1_dwr);
				dYf_dwr.push_back(dyf1_dwr);
				dWnf_dwl.push_back(dwnf2_dwl);
				dYf_dwl.push_back(dyf2_dwl);
				dWnf_dwr.push_back(dwnf2_dwr);
				dYf_dwr.push_back(dyf2_dwr);
				dSlimnf_dwl.push_back(dslimnf1_dwl);
				dSlimnf_dwr.push_back(dslimnf1_dwr);
				dSlimnf_dwl.push_back(dslimnf2_dwl);
				dSlimnf_dwr.push_back(dslimnf2_dwr);
			}
			else {
				Wnf.push_back(wnf2);
				Yf.push_back(yf2);
				Upf.push_back(upf2);
				Slimnf.push_back(slimnf2);
				Wnf.push_back(wnf1);
				Yf.push_back(yf1);
				Upf.push_back(upf1);
				Slimnf.push_back(slimnf1);
				dWnf_dwl.push_back(dwnf2_dwl);
				dYf_dwl.push_back(dyf2_dwl);
				dWnf_dwr.push_back(dwnf2_dwr);
				dYf_dwr.push_back(dyf2_dwr);
				dWnf_dwl.push_back(dwnf1_dwl);
				dYf_dwl.push_back(dyf1_dwl);
				dWnf_dwr.push_back(dwnf1_dwr);
				dYf_dwr.push_back(dyf1_dwr);
				dSlimnf_dwl.push_back(dslimnf2_dwl);
				dSlimnf_dwr.push_back(dslimnf2_dwr);
				dSlimnf_dwl.push_back(dslimnf1_dwl);
				dSlimnf_dwr.push_back(dslimnf1_dwr);
			}
			KWnf.push_back(KWn[i]);
			dKWnf_dwl.push_back(dKWn_dwl[i]);
			dKWnf_dwr.push_back(dKWn_dwr[i]);
			KUpf.push_back(KUp[i]);
			KSlimf.push_back(KSlim[i]);
			KWnf.push_back(KWn[i]);
			dKWnf_dwl.push_back(dKWn_dwl[i]);
			dKWnf_dwr.push_back(dKWn_dwr[i]);
			KUpf.push_back(KUp[i]);
			KSlimf.push_back(KSlim[i]);
		}
	}

	//std::cout << Eigen::Map<Vector>(&s_com[0], s_com.size()).transpose() << std::endl;

	Wnf.push_back(Wn[Wn.size() - 1]);
	Yf.push_back(Y[Y.size() - 1]);
	Upf.push_back(Up[Up.size() - 1]);
	Slimnf.push_back(Slimn[Slimn.size() - 1]);
	dWnf_dwl.push_back(dWn_dwl[dWn_dwl.size() - 1]);
	dWnf_dwr.push_back(dWn_dwr[dWn_dwr.size() - 1]);
	dYf_dwl.push_back(0.);
	dYf_dwr.push_back(0.);
	dSlimnf_dwl.push_back(0.);
	dSlimnf_dwr.push_back(0.);

	//Categorization

	static std::vector<int> intcats{};
	intcats.clear();
	for (size_t i = 0; i != Wnf.size() - 1; i++)
	{
		if (0.5 * (Wnf[i] + Wnf[i + 1]) > 0.5 * (Slimnf[i] + Slimnf[i + 1])) {
			intcats.push_back(0);
		}
		else if (0.5 * (Wnf[i] + Wnf[i + 1]) > eyn) {
			intcats.push_back(1);
		}
		else {
			intcats.push_back(2);
		}
	}

	//Separation into stresses, plastic displacements

	static Vec Sf_new{}; Sf_new.clear();
	static Vec dSf_new_dwl{}; dSf_new_dwl.clear();
	static Vec dSf_new_dwr{}; dSf_new_dwr.clear();
	static Vec Upf_new{}; Upf_new.clear();

	for (size_t i = 0; i != Wnf.size(); i++) {
		if (Wnf[i] > Slimnf[i]) {
			Sf_new.push_back(Slimnf[i] / DU_DS);
			Upf_new.push_back(Upf[i]);
			dSf_new_dwl.push_back(dSlimnf_dwl[i] / DU_DS);
			dSf_new_dwr.push_back(dSlimnf_dwr[i] / DU_DS);
		}
		else if (Wnf[i] > eyn) {
			Sf_new.push_back(Wnf[i] / DU_DS);
			Upf_new.push_back(Upf[i]);
			dSf_new_dwl.push_back(dWnf_dwl[i] / DU_DS);
			dSf_new_dwr.push_back(dWnf_dwr[i] / DU_DS);
		}
		else {
			Sf_new.push_back(ey);
			Upf_new.push_back(Upf[i] + Wnf[i] - eyn);
			dSf_new_dwl.push_back(0.);
			dSf_new_dwr.push_back(0.);
		}
	}

	static Vec DSf{}; DSf.clear();
	static Vec dDSf_dwl{}; dDSf_dwl.clear();
	static Vec dDSf_dwr{}; dDSf_dwr.clear();

	for (size_t i = 0; i != Sf_new.size(); i++) {
		DSf.push_back(Sf_new[i] - Slimnf[i] / DU_DS);
		dDSf_dwl.push_back(dSf_new_dwl[i] - dSlimnf_dwl[i] / DU_DS);
		dDSf_dwr.push_back(dSf_new_dwr[i] - dSlimnf_dwr[i] / DU_DS);
	}

	//Slopes

	static Vec KSf_new{}; KSf_new.clear();
	static Vec dKSf_new_dwl{}; dKSf_new_dwl.clear();
	static Vec dKSf_new_dwr{}; dKSf_new_dwr.clear();
	static Vec KUpf_new{}; KUpf_new.clear();
	static Vec DKSf{}; DKSf.clear();
	static Vec dDKSf_dwl{}; dDKSf_dwl.clear();
	static Vec dDKSf_dwr{}; dDKSf_dwr.clear();

	for (size_t i = 0; i != intcats.size(); i++) {
		if (intcats[i] == 0) {
			KSf_new.push_back(KSlimf[i]);
			dKSf_new_dwl.push_back(0.);
			dKSf_new_dwr.push_back(0.);
			KUpf_new.push_back(KUpf[i]);
			DKSf.push_back(0.);
			dDKSf_dwl.push_back(0.);
			dDKSf_dwr.push_back(0.);
		}
		else if (intcats[i] == 1) {
			KSf_new.push_back(KWnf[i] / DU_DS);
			dKSf_new_dwl.push_back(dKWnf_dwl[i] / DU_DS);
			dKSf_new_dwr.push_back(dKWnf_dwr[i] / DU_DS);
			KUpf_new.push_back(KUpf[i]);
			DKSf.push_back(KSf_new[i] - KSlimf[i]);
			dDKSf_dwl.push_back(dKSf_new_dwl[i]);
			dDKSf_dwr.push_back(dKSf_new_dwr[i]);
		}
		else {
			KSf_new.push_back(0.);
			dKSf_new_dwl.push_back(0.);
			dKSf_new_dwr.push_back(0.);
			KUpf_new.push_back(KUpf[i] + KWnf[i]);
			DKSf.push_back(KSf_new[i] - KSlimf[i]);
			dDKSf_dwl.push_back(dKSf_new_dwl[i]);
			dDKSf_dwr.push_back(dKSf_new_dwr[i]);
		}
	}

	// Simplification

	static double REMLIM = 1.0e-16;

	ys_new.clear();
	s_new.clear();
	ks_new.clear();
	ys_cats.clear();
	dys_new_dwl.clear();
	dys_new_dwr.clear();
	ds_new_dwl.clear();
	ds_new_dwr.clear();
	dks_new_dwl.clear();
	dks_new_dwr.clear();
	yup_new.clear();
	up_new.clear();
	kup_new.clear();
	ydks.clear();
	dydks_dwl.clear();
	dydks_dwr.clear();
	dks.clear();
	ddks_dwl.clear();
	ddks_dwr.clear();
	ds.clear();
	dds_dwl.clear();
	dds_dwr.clear();

	ys_new.push_back(Yf[0]);
	s_new.push_back(Sf_new[0]);
	ks_new.push_back(KSf_new[0]);
	dys_new_dwl.push_back(dYf_dwl[0]);
	ds_new_dwl.push_back(dSf_new_dwl[0]);
	dks_new_dwl.push_back(dKSf_new_dwl[0]);
	dys_new_dwr.push_back(dYf_dwr[0]);
	ds_new_dwr.push_back(dSf_new_dwr[0]);
	dks_new_dwr.push_back(dKSf_new_dwr[0]);
	ys_cats.push_back(intcats[0]);
	for (size_t i = 1; i != Yf.size() - 1; i++) {
		if (std::fabs(KSf_new[i - 1] - KSf_new[i]) <= REMLIM) {
			continue;
		}
		ys_new.push_back(Yf[i]);
		s_new.push_back(Sf_new[i]);
		ks_new.push_back(KSf_new[i]);
		ys_cats.push_back(intcats[i]);
		dys_new_dwl.push_back(dYf_dwl[i]);
		ds_new_dwl.push_back(dSf_new_dwl[i]);
		dks_new_dwl.push_back(dKSf_new_dwl[i]);
		dys_new_dwr.push_back(dYf_dwr[i]);
		ds_new_dwr.push_back(dSf_new_dwr[i]);
		dks_new_dwr.push_back(dKSf_new_dwr[i]);
	}
	ys_new.push_back(Yf[Yf.size() - 1]);
	s_new.push_back(Sf_new[Sf_new.size() - 1]);
	dys_new_dwl.push_back(dYf_dwl[dYf_dwl.size() - 1]);
	ds_new_dwl.push_back(dSf_new_dwl[dSf_new_dwl.size() - 1]);
	dys_new_dwr.push_back(dYf_dwr[dYf_dwr.size() - 1]);
	ds_new_dwr.push_back(dSf_new_dwr[dSf_new_dwr.size() - 1]);

	yup_new.push_back(Yf[0]);
	up_new.push_back(Upf_new[0]);
	kup_new.push_back(KUpf_new[0]);
	for (size_t i = 1; i != Yf.size() - 1; i++) {
		if (std::fabs(KUpf_new[i - 1] - KUpf_new[i]) <= REMLIM) {
			continue;
		}
		yup_new.push_back(Yf[i]);
		up_new.push_back(Upf_new[i]);
		kup_new.push_back(KUpf_new[i]);
	}
	yup_new.push_back(Yf[Yf.size() - 1]);
	up_new.push_back(Upf_new[Upf_new.size() - 1]);

	ydks.push_back(Yf[0]);
	dks.push_back(DKSf[0]);
	ds.push_back(DSf[0]);
	dydks_dwl.push_back(dYf_dwl[0]);
	ddks_dwl.push_back(dDKSf_dwl[0]);
	dydks_dwr.push_back(dYf_dwr[0]);
	ddks_dwr.push_back(dDKSf_dwr[0]);
	dds_dwl.push_back(dDSf_dwl[0]);
	dds_dwr.push_back(dDSf_dwr[0]);
	for (size_t i = 1; i != Yf.size() - 1; i++) {
		if (std::fabs(DKSf[i - 1] - DKSf[i]) <= REMLIM) {
			continue;
		}
		ydks.push_back(Yf[i]);
		ds.push_back(DSf[i]);
		dks.push_back(DKSf[i]);
		dydks_dwl.push_back(dYf_dwl[i]);
		ddks_dwl.push_back(dDKSf_dwl[i]);
		dydks_dwr.push_back(dYf_dwr[i]);
		ddks_dwr.push_back(dDKSf_dwr[i]);
		dds_dwl.push_back(dDSf_dwl[i]);
		dds_dwr.push_back(dDSf_dwr[i]);
	}
	ydks.push_back(Yf[Yf.size() - 1]);
	ds.push_back(DSf[DSf.size() - 1]);
	dydks_dwl.push_back(dYf_dwl[dYf_dwl.size() - 1]);
	dydks_dwr.push_back(dYf_dwr[dYf_dwr.size() - 1]);
	dds_dwl.push_back(dDSf_dwl[dDSf_dwl.size() - 1]);
	dds_dwr.push_back(dDSf_dwr[dDSf_dwr.size() - 1]);

	return;
}

void RockingBC::interval_dists_K(const Vector& Yw, const Vector& W, const Vector& Yup_com, const Vector& Up_com, const Vector& Kup_com, const Vector& Ys_com, const Vector& S_com, const Vector& Ks_com, double ey, double beta_Dt,
	Vector& Ys, Vector& S, Vector& Ks, Vector& Yup_new, Vector& Up_new, Vector& Kup_new,
	Matrix& dYs_dW, Matrix& dS_dW, Matrix& dKs_dW, Vecint& Ys_cats, Vector& Ydks, Vector& Dks, Matrix& dYdks_dW, Matrix& dDks_dW, Vector& DS, Matrix& dDS_dW)
{
	static VecVecOS Yup_ints{}; Yup_ints.clear();
	static VecVecOS Up_ints{}; Up_ints.clear();
	static VecVecOS Kup_ints{}; Kup_ints.clear();
	static VecVecOS Ys_ints{}; Ys_ints.clear();
	static VecVecOS S_ints{}; S_ints.clear();
	static VecVecOS Ks_ints{}; Ks_ints.clear();

	Up_interval_split_K(Yup_com, Up_com, Kup_com, Yw, Yup_ints, Up_ints, Kup_ints);
	Up_interval_split_K(Ys_com, S_com, Ks_com, Yw, Ys_ints, S_ints, Ks_ints);

	//std::cout << Ys_ints[3].transpose() << std::endl;
	//std::cout << S_ints[3].transpose() << std::endl;
	//std::cout << Ks_ints[3].transpose() << std::endl;

	VecVec ys_list(W.Size() - 1, std::vector<double>{});
	VecVec s_list(W.Size() - 1, std::vector<double>{});
	VecVec ks_list(W.Size() - 1, std::vector<double>{});
	VecVec yup_new_list(W.Size() - 1, std::vector<double>{});
	VecVec up_new_list(W.Size() - 1, std::vector<double>{});
	VecVec kup_new_list(W.Size() - 1, std::vector<double>{});

	VecVec dys_dwl_list(W.Size() - 1, std::vector<double>{});
	VecVec ds_dwl_list(W.Size() - 1, std::vector<double>{});
	VecVec dks_dwl_list(W.Size() - 1, std::vector<double>{});
	VecVec dys_dwr_list(W.Size() - 1, std::vector<double>{});
	VecVec ds_dwr_list(W.Size() - 1, std::vector<double>{});
	VecVec dks_dwr_list(W.Size() - 1, std::vector<double>{});
	VecVecint ys_cats_list(W.Size() - 1, std::vector<int>{});

	VecVec ydks_list(W.Size() - 1, std::vector<double>{});
	VecVec dks_list(W.Size() - 1, std::vector<double>{});
	VecVec ds_list(W.Size() - 1, std::vector<double>{});
	VecVec dydks_dwl_list(W.Size() - 1, std::vector<double>{});
	VecVec dydks_dwr_list(W.Size() - 1, std::vector<double>{});
	VecVec ddks_dwl_list(W.Size() - 1, std::vector<double>{});
	VecVec ddks_dwr_list(W.Size() - 1, std::vector<double>{});
	VecVec dds_dwl_list(W.Size() - 1, std::vector<double>{});
	VecVec dds_dwr_list(W.Size() - 1, std::vector<double>{});

	for (size_t i = 0; i != W.Size() - 1; i++) {

		interval_interior_K(W[i], W[i + 1], ey, Yw[i + 1] - Yw[i], Up_ints[i], Yup_ints[i], Kup_ints[i],
			Ys_ints[i], S_ints[i], Ks_ints[i], beta_Dt,
			ys_list[i], s_list[i], ks_list[i], ys_cats_list[i], yup_new_list[i], up_new_list[i], kup_new_list[i],
			dys_dwl_list[i], dys_dwr_list[i], ds_dwl_list[i], ds_dwr_list[i], dks_dwl_list[i], dks_dwr_list[i],
			ydks_list[i], dks_list[i], dydks_dwl_list[i], dydks_dwr_list[i], ddks_dwl_list[i], ddks_dwr_list[i],
			ds_list[i], dds_dwl_list[i], dds_dwr_list[i]);
	}

	Ys = interval_join(ys_list);
	S = interval_join(s_list);
	Ks = array_join(ks_list);
	Yup_new = interval_join(yup_new_list);
	Up_new = interval_join(up_new_list);
	Kup_new = array_join(kup_new_list);

	Ydks = interval_join(ydks_list);
	Dks = array_join(dks_list);
	DS = interval_join(ds_list);

	static VecMatOS dYs_dW_list{}; dYs_dW_list.clear();
	static VecMatOS dS_dW_list{}; dS_dW_list.clear();
	static VecMatOS dKs_dW_list{}; dKs_dW_list.clear();
	static VecMatOS dYdks_dW_list{}; dYdks_dW_list.clear();
	static VecMatOS dDks_dW_list{}; dDks_dW_list.clear();
	static VecMatOS dDS_dW_list{}; dDS_dW_list.clear();

	for (size_t i = 0; i != W.Size() - 1; i++) {

		Vec dwl_dW(W.Size()); dwl_dW[i] = 1.0;
		Vec dwr_dW(W.Size()); dwr_dW[i + 1] = 1.0;
		Matrix dys_dW = Matrix(dys_dwl_list[i].size(), W.Size());
		Matrix ds_dW = Matrix(ds_dwl_list[i].size(), W.Size());
		Matrix dks_dW = Matrix(dks_dwl_list[i].size(), W.Size());
		Matrix dydks_dW = Matrix(dydks_dwl_list[i].size(), W.Size());
		Matrix ddks_dW = Matrix(ddks_dwl_list[i].size(), W.Size());
		Matrix dds_dW = Matrix(dds_dwl_list[i].size(), W.Size());
		for (size_t l = 0; l != W.Size(); l++) {
			for (size_t k = 0; k != dys_dwl_list[i].size(); k++) {
				dys_dW(k, l) += dys_dwl_list[i][k] * dwl_dW[l];
				dys_dW(k, l) += dys_dwr_list[i][k] * dwr_dW[l];
				ds_dW(k, l) += ds_dwl_list[i][k] * dwl_dW[l];
				ds_dW(k, l) += ds_dwr_list[i][k] * dwr_dW[l];
			}
			for (size_t k = 0; k != dks_dwl_list[i].size(); k++) {
				dks_dW(k, l) += dks_dwl_list[i][k] * dwl_dW[l];
				dks_dW(k, l) += dks_dwr_list[i][k] * dwr_dW[l];
			}
			for (size_t k = 0; k != dydks_dwl_list[i].size(); k++) {
				dydks_dW(k, l) += dydks_dwl_list[i][k] * dwl_dW[l];
				dydks_dW(k, l) += dydks_dwr_list[i][k] * dwr_dW[l];
			}
			for (size_t k = 0; k != ddks_dwl_list[i].size(); k++) {
				ddks_dW(k, l) += ddks_dwl_list[i][k] * dwl_dW[l];
				ddks_dW(k, l) += ddks_dwr_list[i][k] * dwr_dW[l];
			}
			for (size_t k = 0; k != dds_dwl_list[i].size(); k++) {
				dds_dW(k, l) += dds_dwl_list[i][k] * dwl_dW[l];
				dds_dW(k, l) += dds_dwr_list[i][k] * dwr_dW[l];
			}
		}

		dYs_dW_list.push_back(dys_dW);
		dS_dW_list.push_back(ds_dW);
		dKs_dW_list.push_back(dks_dW);
		dYdks_dW_list.push_back(dydks_dW);
		dDks_dW_list.push_back(ddks_dW);
		dDS_dW_list.push_back(dds_dW);
	}
	dYs_dW = interval_join(dYs_dW_list);
	dS_dW = interval_join(dS_dW_list);
	dKs_dW = array_join(dKs_dW_list);
	dYdks_dW = interval_join(dYdks_dW_list);
	dDks_dW = array_join(dDks_dW_list);
	dDS_dW = interval_join(dDS_dW_list);

	Ys_cats.clear();
	for (size_t i = 0; i != ys_cats_list.size(); i++) {
		for (size_t j = 0; j != ys_cats_list[i].size(); j++) {
			Ys_cats.push_back(ys_cats_list[i][j]);
		}
	}

	return;

}

void RockingBC::Ys_cats_dist_calc(const VecVecint& Ys_cats, Vecint& Ys_cats_dist) {
	Ys_cats_dist.clear();
	for (size_t i = 0; i != Ys_cats.size(); i++) {
		for (size_t j = 0; j != Ys_cats[i].size(); j++) {
			Ys_cats_dist.push_back(Ys_cats[i][j]);
		}
	}
}

// bilin funcs 

void RockingBC::commony_BL(const Vec& ya, const Vec& fa, const Vec& yb, const Vec& fb, Vec& Y, Vec& FA, Vec& FB)
{
	Y.clear();
	FA.clear();
	FB.clear();

	int ia = 0;
	int ib = 0;
	while ((ia < ya.size() - 1) || (ib < yb.size() - 1))
	{
		if (ya[ia] == yb[ib]) {
			Y.push_back(ya[ia]);
			FA.push_back(fa[ia]);
			FB.push_back(fb[ib]);
			ia += 1;
			ib += 1;
		}
		else if (ya[ia] < yb[ib]) {
			Y.push_back(ya[ia]);
			FA.push_back(fa[ia]);
			FB.push_back(fb[ib - 1] + (ya[ia] - yb[ib - 1]) / (yb[ib] - yb[ib - 1]) * (fb[ib] - fb[ib - 1]));
			ia += 1;
		}
		else {
			Y.push_back(yb[ib]);
			FB.push_back(fb[ib]);
			FA.push_back(fa[ia - 1] + (yb[ib] - ya[ia - 1]) / (ya[ia] - ya[ia - 1]) * (fa[ia] - fa[ia - 1]));
			ib += 1;
		}
	}
	Y.push_back(ya[ya.size() - 1]);
	FA.push_back(fa[fa.size() - 1]);
	FB.push_back(fb[fb.size() - 1]);

	return;
}

bool RockingBC::distintersec(const Vec& YP, const Vec& P, const Vec& YQ, const Vec& Q)
{
	static Vec Y{};
	static Vec PT{};
	static Vec QT{};
	commony_BL(YP, P, YQ, Q, Y, PT, QT);

	int sgn = 0;
	int s = 0;
	for (size_t i = 0; i != PT.size(); i++) {
		if (PT[i] < QT[i]) {
			s = -1;
		}
		else if (PT[i] == QT[i]) {
			s = 0;
		}
		else {
			s = 1;
		}
		if (s != 0 && sgn == 0) {
			sgn = s;
		}
		else if (s != 0 && s != sgn) {
			return true;
		}
	}
	return false;
}

bool RockingBC::twobilinintersec(double y1, double y2, double p1, double p2, double q1, double q2, double yp, double p0, double yq, double q0)
{
	double pe{};
	double qe{};
	if (yp <= yq)
	{
		pe = p0 + (yq - yp) / (y2 - yp) * (p2 - p0);
		qe = q1 + (yp - y1) / (yq - y1) * (q0 - q1);
	}
	else {
		qe = q0 + (yp - yq) / (y2 - yq) * (q2 - q0);
		pe = p1 + (yq - y1) / (yp - y1) * (p0 - p1);
	}
	if (p1 <= q1 && p2 <= q2 && p0 <= qe && pe <= q0) {
		return false;
	}
	else if (p1 >= q1 && p2 >= q2 && p0 >= qe && pe >= q0) {
		return false;
	}
	else {
		return true;
	}
}

void RockingBC::NM_BL(const Vec& Y, const Vec& S, double& N, double& M, double& Nd, double& Md)
{
	N = 0.;
	M = 0.;

	for (size_t i = 0; i != Y.size() - 1; i++) {
		N += (Y[i + 1] - Y[i]) * (S[i + 1] + S[i]) / 2.;
		M += (Y[i + 1] - Y[i]) * (2 * S[i] * Y[i] + S[i] * Y[i + 1] + S[i + 1] * Y[i] + 2 * S[i + 1] * Y[i + 1]) / 6.;
	}

	double Nlin = (Y[Y.size() - 1] - Y[0]) * (S[S.size() - 1] + S[0]) / 2.;
	double Mlin = (Y[Y.size() - 1] - Y[0]) * (2 * S[0] * Y[0] + S[0] * Y[Y.size() - 1] + S[S.size() - 1] * Y[0] + 2 * S[S.size() - 1] * Y[Y.size() - 1]) / 6.;

	Nd = N - Nlin;
	Md = M - Mlin;
}

bool RockingBC::bilinable(double Nd, double Md, double y1, double y2, double BILINLIM) {
	if (fabs(Nd) < BILINLIM && fabs(Md) > BILINLIM) {
		return false;
	}
	else if ((fabs(Nd) < BILINLIM && fabs(Md) < BILINLIM) || (2. * y1 + y2 < 3. * Md / Nd && 3. * Md / Nd < y1 + 2. * y2)) {
		return true;
	}
	else {
		return false;
	}
}

void RockingBC::bilindist(const Vec& Y, const Vec& S, double Nd, double Md, Vec& Ybl, Vec& Sbl, double BILINLIM)
{
	Ybl.clear();
	Sbl.clear();
	if (fabs(Nd) < BILINLIM && fabs(Md) < BILINLIM) {
		Ybl = { Y[0], Y[Y.size() - 1] };
		Sbl = { S[0], S[S.size() - 1] };
		return;
	}

	double s = 2. * Nd / (Y[Y.size() - 1] - Y[0]);
	double y0 = 3. * Md / Nd - Y[0] - Y[Y.size() - 1];
	double k = (S[S.size() - 1] - S[0]) / (Y[Y.size() - 1] - Y[0]);
	double snew = S[0] + (y0 - Y[0]) * k + s;

	Ybl = { Y[0], y0, Y[Y.size() - 1] };
	Sbl = { S[0], snew, S[S.size() - 1] };
	return;

}

bool RockingBC::bilin_two(const Vec& YP, const Vec& P, const Vec& YQ, const Vec& Q, Vec& YPn, Vec& Pn, Vec& YQn, Vec& Qn)
{
	double NP{}, MP{}, NPd{}, MPd{}, NQ{}, MQ{}, NQd{}, MQd{};
	NM_BL(YP, P, NP, MP, NPd, MPd);
	NM_BL(YQ, Q, NQ, MQ, NQd, MQd);

	if (!bilinable(NPd, MPd, YP[0], YP[YP.size() - 1]) || !bilinable(NQd, MQd, YQ[0], YQ[YQ.size() - 1])) {
		return false;
	}

	bilindist(YP, P, NPd, MPd, YPn, Pn);
	bilindist(YQ, Q, NQd, MQd, YQn, Qn);

	double yp{}, p0{}, yq{}, q0{};
	if (YPn.size() == 3) {
		yp = YPn[1];
		p0 = Pn[1];
	}
	else {
		yp = 0.5 * (YPn[0] + YPn[1]);
		p0 = 0.5 * (Pn[0] + Pn[1]);
	}
	if (YQn.size() == 3) {
		yq = YQn[1];
		q0 = Qn[1];
	}
	else {
		yq = 0.5 * (YQn[0] + YQn[1]);
		q0 = 0.5 * (Qn[0] + Qn[1]);
	}

	bool intersec = twobilinintersec(YPn[0], YPn[YPn.size() - 1], Pn[0], Pn[Pn.size() - 1], Qn[0], Qn[Qn.size() - 1], yp, p0, yq, q0);
	if (!intersec) {
		return true;
	}
	else {
		return false;
	}

}

bool RockingBC::bilin_one(const Vec& YP, const Vec& P, Vec& YPn, Vec& Pn) {
	double NP{}, MP{}, NPd{}, MPd{};
	NM_BL(YP, P, NP, MP, NPd, MPd);

	if (!bilinable(NPd, MPd, YP[0], YP[YP.size() - 1])) {
		return false;
	}

	bilindist(YP, P, NPd, MPd, YPn, Pn);
	return true;

}

int 
RockingBC::displaySelf(Renderer& theViewer, int displayMode, float fact, const char** modes, int numMode)
{
	static Vector v1(3);
	static Vector v2(3);

	theNodes[0]->getDisplayCrds(v1, fact, displayMode);
	theNodes[1]->getDisplayCrds(v2, fact, displayMode);

	return theViewer.drawLine(v1, v2, 1.0, 1.0, this->getTag());
}
