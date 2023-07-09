/* ****************************************************************** **
**    Openers - Open System for Earthquake Engineering Simulation    **
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
                                                                        
#include "BeamColumnwLHNMYS.h"
#include <ElementalLoad.h>
#include <algorithm>
#include "GPYS2d.h"

#include <Domain.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>

#include <CrdTransf.h>
#include <SectionForceDeformation.h>
#include <Information.h>
#include <Parameter.h>
#include <ElementResponse.h>
#include <Renderer.h>

#include <math.h>
#include <stdlib.h>
#include <elementAPI.h>
#include <string>
#include <ElementIter.h>
#include <iostream>

Matrix BeamColumnwLHNMYS::K(6,6);
Vector BeamColumnwLHNMYS::P(6);
Matrix BeamColumnwLHNMYS::k(3,3);

void* OPS_BeamColumnwLHNMYS(void)
{
    if(OPS_GetNumRemainingInputArgs() < 11) {
	opserr<<"insufficient arguments:eleTag,iNode,jNode,A,E,Iz,NpI NpJ,MpI MpJ,transfTag\n";
	return 0;
    }

    int ndm = OPS_GetNDM();
    int ndf = OPS_GetNDF();
    if(ndm != 2 || ndf != 3) {
	opserr<<"ndm must be 2 and ndf must be 3\n";
	return 0;
    }

    // inputs: 
    int iData[3];
    int numData = 3;
    if(OPS_GetIntInput(&numData,&iData[0]) < 0) {
	opserr<<"WARNING failed to read integers\n";
	return 0;
    }

    double data[7];
    // Read A, E, I, NpI, NpJ, MpI, MpJ
    numData = 7;
    if(OPS_GetDoubleInput(&numData,&data[0]) < 0) {
	opserr<<"WARNING failed to read doubles\n";
	return 0;
      }

    numData = 1;
    int transfTag;
    if(OPS_GetIntInput(&numData,&transfTag) < 0) {
	opserr<<"WARNING transfTag is not integer\n";
	return 0;
    }
    
    // options
    double HirI = 0.0, HirJ = 0.0, HkrA = 0.0, HkrI = 0.0, HkrJ =0.0;
    double mass = 0.0, Wtol = 1e-16, yftol = 1e-8;
    int cMass = 0, MaxIter = 15, nrow = 4;
    double coefData[3*4];
    for (int i=0; i<3*nrow; i++){
        coefData[i] = 0.0;
    }
    coefData[0] = 1.0;
    coefData[1] = 1.0;
    coefData[2] = 3.5;
    coefData[3] = -1.0;
    coefData[4] = 2.0;
    coefData[6] = 2.0;
    coefData[9] = 2.0;
    coefData[10] = 2.0;
    
    
    while(OPS_GetNumRemainingInputArgs() > 0) {
	std::string type = OPS_GetString();
	if(type == "-Wtol") {
	    if(OPS_GetNumRemainingInputArgs() > 0) {
		if(OPS_GetDoubleInput(&numData,&Wtol) < 0) return 0;
	    }
    } else if(type == "-yftol") {
        if(OPS_GetNumRemainingInputArgs() > 0) {
            if(OPS_GetDoubleInput(&numData,&yftol) < 0) return 0;
        }
    } else if(type == "-MaxIter") {
        if(OPS_GetNumRemainingInputArgs() > 0) {
            if(OPS_GetIntInput(&numData,&MaxIter) < 0) return 0;
        }
	} else if(type == "-HirI") {
	    if(OPS_GetNumRemainingInputArgs() > 0) {
		if(OPS_GetDoubleInput(&numData,&HirI) < 0) return 0;
	    }
    } else if(type == "-HirJ") {
        if(OPS_GetNumRemainingInputArgs() > 0) {
            if(OPS_GetDoubleInput(&numData,&HirJ) < 0) return 0;
        }
    } else if(type == "-HkrA") {
        if(OPS_GetNumRemainingInputArgs() > 0) {
            if(OPS_GetDoubleInput(&numData,&HkrA) < 0) return 0;
        }
    } else if(type == "-HkrI") {
        if(OPS_GetNumRemainingInputArgs() > 0) {
            if(OPS_GetDoubleInput(&numData,&HkrI) < 0) return 0;
        }
    } else if(type == "-HkrJ") {
        if(OPS_GetNumRemainingInputArgs() > 0) {
            if(OPS_GetDoubleInput(&numData,&HkrJ) < 0) return 0;
        }
	} else if(type == "-mass") {
	    if(OPS_GetNumRemainingInputArgs() > 0) {
            if(OPS_GetDoubleInput(&numData,&mass) < 0) return 0;
	    }
	} else if(type == "-cMass") {
	    cMass = 1;
    } else if(type == "-ySurf") {
        if(OPS_GetNumRemainingInputArgs() > 0) {
            if(OPS_GetIntInput(&numData,&nrow) < 0) return 0;
            numData = 3*nrow;
            if(OPS_GetDoubleInput(&numData,&coefData[0]) < 0) return 0;
            numData = 1;
        }
    }
    }
    
    // Create vector for isotropic hardening ratio
    Vector Hir(2);
    Hir(0) = HirI;
    Hir(1) = HirJ;
    
    // Create vector for kinematic hardening ratio
    Vector Hkr(3);
    Hkr(0) = HkrA;
    Hkr(1) = HkrI;
    Hkr(2) = HkrJ;
        
    // Create matrix with coefficients of yield surface
    Matrix GPYSC(coefData,nrow,3);

    // check transf
    CrdTransf* theTransf = OPS_getCrdTransf(transfTag);
    if(theTransf == 0) {
	opserr<<"coord transformation not found\n";
	return 0;
    }
    
    return new BeamColumnwLHNMYS(iData[0],iData[1],iData[2],data[0],data[1],data[2],data[3],data[4], data[5],data[6],*theTransf,yftol,Wtol,MaxIter,Hir,Hkr,nrow,GPYSC,mass,cMass);
    
}


BeamColumnwLHNMYS::BeamColumnwLHNMYS()
  :Element(0,ELE_TAG_ElasticBeam2d),
  A(0.0), E(0.0), I(0.0), rho(0.0), cMass(0),
  Q(6), q(3), connectedExternalNodes(2), theCoordTransf(0)
{
  // does nothing
  q0[0] = 0.0;
  q0[1] = 0.0;
  q0[2] = 0.0;

  p0[0] = 0.0;
  p0[1] = 0.0;
  p0[2] = 0.0;

  // set node pointers to NULL
  for (int i=0; i<2; i++)
    theNodes[i] = 0;      
}

BeamColumnwLHNMYS::BeamColumnwLHNMYS(int tag, int Nd1, int Nd2, double a, double e,
                                     double i, double npi, double npj, double mpi, double mpj,
                                     CrdTransf &coordTransf, double ytol,
                                     double wtol, double maxiter, Vector &hir, Vector &hkr,
                                     double nr, Matrix &gpysc, double r, int cm)
  :Element(tag,ELE_TAG_ElasticBeam2d), Wtol(wtol),MaxIter(maxiter), yftol(ytol),
A(a), E(e), I(i), NpI(npi),NpJ(npj), MpI(mpi), MpJ(mpj), Hir(2), Hkr(3), nrow(nr), GPYSC(nr,3),
rho(r), cMass(cm), vppast(3), vppres(3),alphapast(2),alphapres(2),qbpast(3),qbpres(3),
Q(6), q(3), connectedExternalNodes(2), theCoordTransf(0)
{
    connectedExternalNodes(0) = Nd1;
    connectedExternalNodes(1) = Nd2;

    Hir = hir;
    Hkr = hkr;
    GPYSC = gpysc;
    
    theCoordTransf = coordTransf.getCopy2d();
    
    if (!theCoordTransf) {
        opserr << "BeamColumnwLHNMYS::BeamColumnwLHNMYS -- failed to get copy of coordinate transformation\n";
        exit(01);
    }

    q0[0] = 0.0;
    q0[1] = 0.0;
    q0[2] = 0.0;

    p0[0] = 0.0;
    p0[1] = 0.0;
    p0[2] = 0.0;

    // set node pointers to NULL
    theNodes[0] = 0;
    theNodes[1] = 0;
    
    // Initialize history variables
    vppast.Zero();
    vppres.Zero();
    alphapast.Zero();
    alphapres.Zero();
    qbpast.Zero();
    qbpres.Zero();
    
}

BeamColumnwLHNMYS::~BeamColumnwLHNMYS()
{
    if (theCoordTransf)
	delete theCoordTransf;
}

int
BeamColumnwLHNMYS::getNumExternalNodes(void) const
{
    return 2;
}

const ID &
BeamColumnwLHNMYS::getExternalNodes(void)
{
    return connectedExternalNodes;
}

Node **
BeamColumnwLHNMYS::getNodePtrs(void)
{
  return theNodes;
}

int
BeamColumnwLHNMYS::getNumDOF(void)
{
    return 6;
}

void
BeamColumnwLHNMYS::setDomain(Domain *theDomain)
{
  if (theDomain == 0) {
    opserr << "BeamColumnwLHNMYS::setDomain -- Domain is null\n";
    exit(-1);
  }
    
    theNodes[0] = theDomain->getNode(connectedExternalNodes(0));
    theNodes[1] = theDomain->getNode(connectedExternalNodes(1));    
    
    if (theNodes[0] == 0) {
      opserr << "BeamColumnwLHNMYS::setDomain -- Node 1: " << connectedExternalNodes(0) << " does not exist\n";
      exit(-1);
    }
			      
    if (theNodes[1] == 0) {
      opserr << "BeamColumnwLHNMYS::setDomain -- Node 2: " << connectedExternalNodes(1) << " does not exist\n";
      exit(-1);
    }

    int dofNd1 = theNodes[0]->getNumberDOF();
    int dofNd2 = theNodes[1]->getNumberDOF();    
    
    if (dofNd1 != 3) {
      opserr << "BeamColumnwLHNMYS::setDomain -- Node 1: " << connectedExternalNodes(0)
	     << " has incorrect number of DOF\n";
      exit(-1);
    }
    
    if (dofNd2 != 3) {
      opserr << "BeamColumnwLHNMYS::setDomain -- Node 2: " << connectedExternalNodes(1)
	     << " has incorrect number of DOF\n";
      exit(-1);
    }
	
    this->DomainComponent::setDomain(theDomain);
    
    if (theCoordTransf->initialize(theNodes[0], theNodes[1]) != 0) {
	opserr << "BeamColumnwLHNMYS::setDomain -- Error initializing coordinate transformation\n";
	exit(-1);
    }
    
    double L = theCoordTransf->getInitialLength();

    if (L == 0.0) {
      opserr << "BeamColumnwLHNMYS::setDomain -- Element has zero length\n";
      exit(-1);
    }
}

int
BeamColumnwLHNMYS::commitState()
{
  vppast = vppres;
  alphapast = alphapres;
  qbpast = qbpres;

  int retVal = 0;
  // call element commitState to do any base class stuff
  if ((retVal = this->Element::commitState()) != 0) {
    opserr << "BeamColumnwLHNMYS::commitState () - failed in base class";
  }    
  retVal += theCoordTransf->commitState();
  return retVal;
}

int
BeamColumnwLHNMYS::revertToLastCommit()
{
    return theCoordTransf->revertToLastCommit();
}

int
BeamColumnwLHNMYS::revertToStart()
{
    return theCoordTransf->revertToStart();
}

int
BeamColumnwLHNMYS::update(void)
{
    int ok = theCoordTransf->update();
    const Vector &v = theCoordTransf->getBasicTrialDisp();
    
    double L = theCoordTransf->getInitialLength();
    double EA = E*A;
    double EI = E*I;
    double EAoverL = EA/L;
    double EIoverL = EI/L;
    double EIoverL2 = 2.0*EIoverL;        // 2EI/L
    double EIoverL4 = 2.0*EIoverL2;        // 4EI/L
    double LoverEA = L/EA;
    double LoverEI = L/EI;
    double LoverEI3 = LoverEI/3.0;      //  L/(3EI)
    double mLoverEI6 = -LoverEI3/2.0;   // -L/(6EI)
    
    // elastic stiffness matrix
    static Matrix ke(3,3);
    ke.Zero();
    ke(0,0) = EAoverL;
    ke(1,1) = ke(2,2) = EIoverL4;
    ke(2,1) = ke(1,2) = EIoverL2;
    
    // elastic flexibility matrix
    static Matrix fe(3,3);
    fe.Zero();
    fe(0,0) = LoverEA;
    fe(1,1) = fe(2,2) = LoverEI3;
    fe(2,1) = fe(1,2) = mLoverEI6;
    
    // kinematic hardening matrix
    static Matrix Hk(3,3);
    Hk.Zero();
    Hk(0,0) = EAoverL*Hkr(0);
    Hk(1,1) = 6.0*EIoverL*Hkr(1);
    Hk(2,2) = 6.0*EIoverL*Hkr(2);
    
    // isotropic hardening matrix
    static Matrix Hi(2,2);
    Hi.Zero();
    Hi(0,0) = Hir(0);
    Hi(1,1) = Hir(1);
    
    // Recover information from past
    Vector vp(vppast);
    Vector alpha(alphapast);
    Vector qb(qbpast);

    // trial elastic step
    static Vector vmvp(3);
    vmvp = v - vp;
    static Vector qtr(3);
    qtr = ke*vmvp;
    Vector xyrefI(2), xyrefJ(2);
    Vector ScVecI(2), ScVecJ(2);
    xyrefI(0) = qtr(0)-qb(0);
    xyrefI(1) = qtr(1)-qb(1);
    xyrefJ(0) = qtr(0)-qb(0);
    xyrefJ(1) = qtr(2)-qb(2);
    ScVecI(0) = NpI;
    ScVecI(1) = MpI;
    ScVecJ(0) = NpJ;
    ScVecJ(1) = MpJ;
    double yfnI, yfnJ;
    Vector gI(2), gJ(2);
    Matrix HI(2,2), HJ(2,2);
    GPYS2d(GPYSC,xyrefI,ScVecI,yfnI,gI,HI);
    GPYS2d(GPYSC,xyrefJ,ScVecJ,yfnJ,gJ,HJ);
    Vector yftr(2);
    yftr(0) = yfnI - (Hir(0)*alpha(0));
    yftr(1) = yfnJ - (Hir(1)*alpha(1));
    Matrix g(3,2);
    g.Zero();
    g(0,0) = gI(0);
    g(1,0) = gI(1);
    g(0,1) = gJ(0);
    g(2,1) = gJ(1);
    Matrix Hmat1(3,3);
    Hmat1.Zero();
    Hmat1(0,0) = HI(0,0);
    Hmat1(1,0) = HI(1,0);
    Hmat1(0,1) = HI(0,1);
    Hmat1(1,1) = HI(1,1);
    Matrix Hmat2(3,3);
    Hmat2.Zero();
    Hmat2(0,0) = HJ(0,0);
    Hmat2(2,0) = HJ(1,0);
    Hmat2(0,2) = HJ(0,1);
    Hmat2(2,2) = HJ(1,1);
    // Check location of trial force relative to yield surface
    if (yftr(0) <= yftol && yftr(1) <= yftol) {
        // trial state is admissible (falls inside the yield surface)
        k = ke;
        q = qtr;
    }
    else {
        // else, iterative solution
        int iter = 0;
        int le;
        double DWr = Wtol + 1.0;
        double DW0,DW;
        ID id3(3),Jact(2),Jinact(1);
        id3(0) = 0;
        id3(1) = 1;
        id3(2) = 2;
        Vector Dbeta(2),Rv(3),Ra(2),Rb(3),Dvp(3),R(10),x(10);
        Vector DDq(3),DDvp(3),Dalpha(2),Dqb(3),DDbeta(2),Dbeta1(2);
        Dbeta.Zero();
        Rv.Zero();
        Ra.Zero();
        Rb.Zero();
        Dvp.Zero();
        Matrix g1(3,2),Qk(3,3), feQk(3,3), A(3,3), invA(3,3), HQ(3,3);
        Matrix eye3(3,3), z3(3,3),IHQ(3,3),B(3,3),invB(3,3);
        Matrix J(10,10),Frc(3,7),Fcr(7,3),Fcc(7,7),FccFcr(7,3),Frr(3,3);
        Matrix zle, eyele, zle3, z3le, Hkg;
        z3.Zero();
        eye3.Zero();
        eye3(0,0) = eye3(1,1) = eye3(2,2) = 1.0;
        // Find active end node
        if ((yftr(0) > yftol) && (yftr(1) > yftol)){
            le = 2;
            Jact.resize(le);
            Jact(0) = 0;
            Jact(1) = 1;
        }
        else {
            le = 1;
            Jact.resize(le);
            if (yftr(0) > yftol) {
                Jact(0) = 0;
            }
            else {
                Jact(0) = 1;
            }
        }
        while (DWr > Wtol) {
            iter += 1;
            if (iter > MaxIter) {
                opserr << "Element reached maximum no. of iterations \n";
                return -1;
            }
            // Update Hessian matrix
            Qk = (Dbeta(0)*Hmat1) + (Dbeta(1)*Hmat2);
            feQk = fe + Qk;
            feQk.Invert(A);
            A.Invert(invA);
            HQ = Hk*Qk;
            // Determine plastic flow increment
            DDbeta.Zero();
            zle.resize(le,le);
            eyele.resize(le,le);
            zle3.resize(le,3);
            z3le.resize(3,le);
            zle.Zero();
            zle3.Zero();
            z3le.Zero();
            eyele.Zero();
            for (int i=0; i<le; i++){
                eyele(i,i) = 1.0;
            }
            R.resize(6+le+le);
            R.Zero();
            R.Assemble(Rv,0);
            R.Assemble(yftr(Jact),3);
            R.Assemble(Ra(Jact),3+le);
            R.Assemble(Rb,3+le+le);
            g1.resize(3,le);
            g1 = g(id3,Jact);
            Hkg = Hk*g1;
            IHQ = eye3 + HQ;
            // Jacobian
            J.resize(6+le+le,6+le+le);
            J.Zero();
            J.Assemble(g1,0,0,-1.0);
            J.Assemble(invA,0,le,-1.0);
            J.Assemble(z3le,0,3+le);
            J.Assemble(Qk,0,3+le+le);
            J.Assemble(zle,3,0);
            J.AssembleTranspose(g1,3,le,-1.0);
            J.Assemble(Hi(Jact,Jact),3,3+le);
            J.AssembleTranspose(g1,3,3+le+le);
            J.Assemble(eyele,3+le,0,-1.0);
            J.Assemble(zle3,3+le,le);
            J.Assemble(eyele,3+le,3+le);
            J.Assemble(zle3,3+le,3+le+le);
            J.Assemble(Hkg,3+le+le,0,-1.0);
            J.Assemble(HQ,3+le+le,le,-1.0);
            J.Assemble(z3le,3+le+le,3+le);
            J.Assemble(IHQ,3+le+le,3+le+le);
            // Compute increments
            x.resize(6+le+le);
            J.Solve(R,x);
            for (int i=0; i<le; i++){
                DDbeta(Jact(i)) = x(i);
            }
            
            // Test for unloading of active node
            Dbeta1.resize(le);
            Dbeta1 = Dbeta + DDbeta;
            if ((le == 2) && ((Dbeta1(0) < 0) || (Dbeta1(1) < 0))){
                le = 1;
                Jact.resize(le);
                Jinact.resize(le);
                if ((Dbeta1(0) > 0)){
                    Jact(0) = 0;
                    Jinact(0) = 1;
                } else {
                    Jact(0) = 1;
                    Jinact(0) = 0;
                }
                //  Define residual vector
                R.resize(6+le+le);
                R.Zero();
                R.Assemble(Rv,0);
                R.Assemble(yftr(Jact),3);
                R.Assemble(Ra(Jact),3+le);
                R.Assemble(Rb,3+le+le);
                g1.resize(3,le);
                g1 = g(id3,Jact);
                Hkg = Hk*g1;
                IHQ = eye3 + HQ;
                // Jacobian
                J.resize(6+le+le,6+le+le);
                J.Zero();
                J.Assemble(g1,0,0,-1.0);
                J.Assemble(invA,0,le,-1.0);
                J.Assemble(z3le,0,3+le);
                J.Assemble(Qk,0,3+le+le);
                J.Assemble(zle,3,0);
                J.AssembleTranspose(g1,3,le,-1.0);
                J.Assemble(Hi(Jact,Jact),3,3+le);
                J.AssembleTranspose(g1,3,3+le+le);
                J.Assemble(eyele,3+le,0,-1.0);
                J.Assemble(zle3,3+le,le);
                J.Assemble(eyele,3+le,3+le);
                J.Assemble(zle3,3+le,3+le+le);
                J.Assemble(Hkg,3+le+le,0,-1.0);
                J.Assemble(HQ,3+le+le,le,-1.0);
                J.Assemble(z3le,3+le+le,3+le);
                J.Assemble(IHQ,3+le+le,3+le+le);
                // Compute increments
                x.resize(6+le+le);
                J.Solve(R,x);
                for (int i=0; i<le; i++){
                    DDbeta(Jact(i)) = x(i);
                    // Reset plastic deformation increment of non-active end
                    DDbeta(Jinact(i)) = 0.0;
                }
            }
            
            // Determine element force increment
            Vector DDq(3);
            DDq.Extract(x,le);
            // Determine plastic deformation increment
            Vector DDvp(3);
            DDvp = -1.0*(fe*DDq);
            // Update plastic deformation and plastic flow increment
            Dvp += DDvp;
            Dbeta += DDbeta;
            // Update hardening variables
            Vector Dalpha(2);
            Dalpha.Zero();
            for (int i=0; i<le; i++){
                Dalpha(Jact(i)) = x(i+le+3);
            }
            Vector Dqb(3);
            Dqb.Extract(x,3+le+le);
            alpha += Dalpha;
            qb += Dqb;
            // Correct end force vector
            q = qtr - ke*Dvp;
            // Update yield function value, gradient and Hessian
            xyrefI(0) = q(0)-qb(0);
            xyrefI(1) = q(1)-qb(1);
            xyrefJ(0) = q(0)-qb(0);
            xyrefJ(1) = q(2)-qb(2);
            GPYS2d(GPYSC,xyrefI,ScVecI,yfnI,gI,HI);
            GPYS2d(GPYSC,xyrefJ,ScVecJ,yfnJ,gJ,HJ);
            yftr(0) = yfnI - (Hir(0)*alpha(0));
            yftr(1) = yfnJ - (Hir(1)*alpha(1));
            g(0,0) = gI(0);
            g(1,0) = gI(1);
            g(0,1) = gJ(0);
            g(2,1) = gJ(1);
            Hmat1(0,0) = HI(0,0);
            Hmat1(1,0) = HI(1,0);
            Hmat1(0,1) = HI(0,1);
            Hmat1(1,1) = HI(1,1);
            Hmat2(0,0) = HJ(0,0);
            Hmat2(2,0) = HJ(1,0);
            Hmat2(0,2) = HJ(0,1);
            Hmat2(2,2) = HJ(1,1);
            // Plastic deformation residual
            Rv = g*Dbeta - Dvp;
            Ra = alphapast - alpha + Dbeta;
            Rb = qbpast - qb + (Hk*g*Dbeta);
            DW = fabs((DDq^Rv) + (DDbeta^yftr) + (Dalpha^Ra) + (Dqb^Rb));
            if (iter==1) DW0 = std::max(DW,1e-6);
            DWr = DW/DW0;
        }
        // Update plastic deformation
        vp += Dvp;
        // algorithmic tangent element stiffness
        g1 = g(id3,Jact);
        Frc.resize(3,3+le+le);
        Fcr.resize(3+le+le,3);
        Fcc.resize(3+le+le,3+le+le);
        Frc.Zero();
        Fcr.Zero();
        Fcc.Zero();
        Frc.Assemble(g1,0,0,-1.0);
        Frc.Assemble(z3le,0,le);
        Frc.Assemble(Qk,0,le+le);
        Fcr.AssembleTranspose(g1,0,0,-1.0);
        Fcr.Assemble(zle3,le,0);
        Fcr.Assemble(HQ,le+le,0,-1.0);
        Fcc.Assemble(zle,0,0);
        Fcc.Assemble(Hi(Jact,Jact),0,le);
        Fcc.AssembleTranspose(g1,0,le+le);
        Fcc.Assemble(eyele,le,0,-1.0);
        Fcc.Assemble(eyele,le,le);
        Fcc.Assemble(zle3,le,le+le);
        Fcc.Assemble(Hkg,le+le,0,-1.0);
        Fcc.Assemble(z3le,le+le,le);
        Fcc.Assemble(IHQ,le+le,le+le);
        FccFcr.resize(3+le+le,3);
        Frr = -1.0*invA;
        Fcc.Solve(Fcr,FccFcr);
        B = Frr - Frc*FccFcr;
        B.Invert(invB);
        k = -1.0*invB;
    }
    
    // save element history before returning
    vppres = vp;
    alphapres = alpha;
    qbpres = qb;
    
 return ok;
    
}

const Matrix &
BeamColumnwLHNMYS::getTangentStiff(void)
{
  return theCoordTransf->getGlobalStiffMatrix(k, q);
}

const Matrix &
BeamColumnwLHNMYS::getInitialStiff(void)
{
    double L = theCoordTransf->getInitialLength();

    double EA = E*A;
    double EI = E*I;
    double EAoverL = EA/L;
    double EIoverL = EI/L;
    double EIoverL2 = 2.0*EIoverL;		// 2EI/L
    double EIoverL4 = 2.0*EIoverL2;		// 4EI/L
    
    // elastic stiffness matrix
    static Matrix ke(3,3);
    ke.Zero();
    ke(0,0) = EAoverL;
    ke(1,1) = ke(2,2) = EIoverL4;
    ke(2,1) = ke(1,2) = EIoverL2;
  
    return theCoordTransf->getInitialGlobalStiffMatrix(ke);
}

const Matrix &
BeamColumnwLHNMYS::getMass(void)
{ 
    K.Zero();
    
    if (rho > 0.0)  {
        // get initial element length
        double L = theCoordTransf->getInitialLength();
        if (cMass == 0)  {
            // lumped mass matrix
            double m = 0.5*rho*L;
            K(0,0) = K(1,1) = K(3,3) = K(4,4) = m;
        } else  {
            // consistent mass matrix
            static Matrix ml(6,6);
            double m = rho*L/420.0;
            ml(0,0) = ml(3,3) = m*140.0;
            ml(0,3) = ml(3,0) = m*70.0;

            ml(1,1) = ml(4,4) = m*156.0;
            ml(1,4) = ml(4,1) = m*54.0;
            ml(2,2) = ml(5,5) = m*4.0*L*L;
            ml(2,5) = ml(5,2) = -m*3.0*L*L;
            ml(1,2) = ml(2,1) = m*22.0*L;
            ml(4,5) = ml(5,4) = -ml(1,2);
            ml(1,5) = ml(5,1) = -m*13.0*L;
            ml(2,4) = ml(4,2) = -ml(1,5);
            
            // transform local mass matrix to global system
            K = theCoordTransf->getGlobalMatrixFromLocal(ml);
        }
    }
    
    return K;
}

void 
BeamColumnwLHNMYS::zeroLoad(void)
{
  Q.Zero();

  q0[0] = 0.0;
  q0[1] = 0.0;
  q0[2] = 0.0;

  p0[0] = 0.0;
  p0[1] = 0.0;
  p0[2] = 0.0;

  return;
}

int 
BeamColumnwLHNMYS::addLoad(ElementalLoad *theLoad, double loadFactor)
{
  int type;
  const Vector &data = theLoad->getData(type, loadFactor);
  double L = theCoordTransf->getInitialLength();

  if (type == LOAD_TAG_Beam2dUniformLoad) {
    double wt = data(0)*loadFactor;  // Transverse (+ve upward)
    double wa = data(1)*loadFactor;  // Axial (+ve from node I to J)

    double V = 0.5*wt*L;
    double M = V*L/6.0; // wt*L*L/12
    double P = wa*L;

    // Reactions in basic system
    p0[0] -= P;
    p0[1] -= V;
    p0[2] -= V;

    // Fixed end forces in basic system
    q0[0] -= 0.5*P;
    q0[1] -= M;
    q0[2] += M;
  }

  else if (type == LOAD_TAG_Beam2dPointLoad) {
    double P = data(0)*loadFactor;
    double N = data(1)*loadFactor;
    double aOverL = data(2);

    if (aOverL < 0.0 || aOverL > 1.0)
      return 0;

    double a = aOverL*L;
    double b = L-a;

    // Reactions in basic system
    p0[0] -= N;
    double V1 = P*(1.0-aOverL);
    double V2 = P*aOverL;
    p0[1] -= V1;
    p0[2] -= V2;

    double L2 = 1.0/(L*L);
    double a2 = a*a;
    double b2 = b*b;

    // Fixed end forces in basic system
    q0[0] -= N*aOverL;
    double M1 = -a * b2 * P * L2;
    double M2 = a2 * b * P * L2;
    q0[1] += M1;
    q0[2] += M2;
  }
  
  else if (type == LOAD_TAG_Beam2dTempLoad) {
    double Ttop1 = data(0)* loadFactor;
    double Tbot1 = data(1)* loadFactor;
    double Ttop2 = data(2)* loadFactor;
    double Tbot2 = data(3)* loadFactor;
        
    // fixed end forces due to a linear thermal load
    double dT1 = Ttop1-Tbot1;
    double dT = (Ttop2-Tbot2)-(Ttop1-Tbot1);
    ///////
    double alpha = 0;
    double d = 1;
    ///////
    double a = alpha/d;  // constant based on temp difference at top and bottom, 
    // coefficient of thermal expansion and beam depth
    double M1 = a*E*I*(-dT1+(4.0/3.0)*dT); //Fixed End Moment end 1
    double M2 = a*E*I*(dT1+(5.0/3.0)*dT); //Fixed End Moment end 2
    double F = alpha*(((Ttop2+Ttop1)/2+(Tbot2+Tbot1)/2)/2)*E*A; // Fixed End Axial Force
    double M1M2divL =(M1+M2)/L; // Fixed End Shear
    
    // Reactions in basic system
    p0[0] += 0;
    p0[1] += M1M2divL;
    p0[2] -= M1M2divL;

    // Fixed end forces in basic system
    q0[0] -= F;
    q0[1] += M1;
    q0[2] += M2;
  }

  else {
    opserr << "BeamColumnwLHNMYS::addLoad()  -- load type unknown for element with tag: " << this->getTag() << endln;
    return -1;
  }

  return 0;
}

int
BeamColumnwLHNMYS::addInertiaLoadToUnbalance(const Vector &accel)
{
  if (rho == 0.0)
    return 0;

  // get R * accel from the nodes
  const Vector &Raccel1 = theNodes[0]->getRV(accel);
  const Vector &Raccel2 = theNodes[1]->getRV(accel);
	
  if (3 != Raccel1.Size() || 3 != Raccel2.Size()) {
    opserr << "BeamColumnwLHNMYS::addInertiaLoadToUnbalance matrix and vector sizes are incompatible\n";
    return -1;
  }
    
  // want to add ( - fact * M R * accel ) to unbalance
  if (cMass == 0)  {
    // take advantage of lumped mass matrix
    double L = theCoordTransf->getInitialLength();
    double m = 0.5*rho*L;

    Q(0) -= m * Raccel1(0);
    Q(1) -= m * Raccel1(1);

    Q(3) -= m * Raccel2(0);
    Q(4) -= m * Raccel2(1);
  } else  {
    // use matrix vector multip. for consistent mass matrix
    static Vector Raccel(6);
    for (int i=0; i<3; i++)  {
      Raccel(i)   = Raccel1(i);
      Raccel(i+3) = Raccel2(i);
    }
    Q.addMatrixVector(1.0, this->getMass(), Raccel, -1.0);
  }
  
  return 0;
}

const Vector &
BeamColumnwLHNMYS::getResistingForceIncInertia()
{	
  P = this->getResistingForce();
  
  // subtract external load P = P - Q
  P.addVector(1.0, Q, -1.0);
  
  // add the damping forces if rayleigh damping
  if (alphaM != 0.0 || betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0)
    P.addVector(1.0, this->getRayleighDampingForces(), 1.0);
    
  if (rho == 0.0)
    return P;

  // add inertia forces from element mass
  const Vector &accel1 = theNodes[0]->getTrialAccel();
  const Vector &accel2 = theNodes[1]->getTrialAccel();    
  
  if (cMass == 0)  {
    // take advantage of lumped mass matrix
    double L = theCoordTransf->getInitialLength();
    double m = 0.5*rho*L;

    P(0) += m * accel1(0);
    P(1) += m * accel1(1);

    P(3) += m * accel2(0);
    P(4) += m * accel2(1);
  } else  {
    // use matrix vector multip. for consistent mass matrix
    static Vector accel(6);
    for (int i=0; i<3; i++)  {
      accel(i)   = accel1(i);
      accel(i+3) = accel2(i);
    }
    P.addMatrixVector(1.0, this->getMass(), accel, 1.0);
  }
  
  return P;
}


const Vector &
BeamColumnwLHNMYS::getResistingForce()
{
  theCoordTransf->update();
  
  q(0) += q0[0];
  q(1) += q0[1];
  q(2) += q0[2];
  
  // Vector for reactions in basic system
  Vector p0Vec(p0, 3);
  
  P = theCoordTransf->getGlobalResistingForce(q, p0Vec);

  return P;
}

int
BeamColumnwLHNMYS::sendSelf(int cTag, Channel &theChannel)
{
    return -1;
}

int
BeamColumnwLHNMYS::recvSelf(int cTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
    return -1;
}

void
BeamColumnwLHNMYS::Print(OPS_Stream &s, int flag)
{
  // to update forces!
  this->getResistingForce();

  if (flag == -1) {
    int eleTag = this->getTag();
    s << "EL_BEAM\t" << eleTag << "\t";
    s << 0 << "\t" << 0 << "\t" << connectedExternalNodes(0) << "\t" << connectedExternalNodes(1) ;
    s << "0\t0.0000000\n";
  }

  if (flag == OPS_PRINT_CURRENTSTATE) {
    this->getResistingForce();
    s << "\nBeamColumnwLHNMYS: " << this->getTag() << endln;
    s << "\tConnected Nodes: " << connectedExternalNodes ;
    s << "\tCoordTransf: " << theCoordTransf->getTag() << endln;
    s << "\tmass density:  " << rho << ", cMass: " << cMass << endln;
    s << "\tAxial capacity End 1: " << NpI << ", End 2: " << NpJ << endln;
    s << "\tMoment capacity End 1: " << MpI << ", End 2: " << MpJ << endln;
    double P  = q(0);
    double M1 = q(1);
    double M2 = q(2);
    double L = theCoordTransf->getInitialLength();
    double V = (M1+M2)/L;
    s << "\tEnd 1 Forces (P V M): " << -P+p0[0]
      << " " << V+p0[1] << " " << M1 << endln;
    s << "\tEnd 2 Forces (P V M): " << P
      << " " << -V+p0[2] << " " << M2 << endln;
  }

  if (flag == OPS_PRINT_PRINTMODEL_JSON) {
    s << "\t\t\t{";
	s << "\"name\": " << this->getTag() << ", ";
	s << "\"type\": \"BeamColumnwLHNMYS\", ";
    s << "\"nodes\": [" << connectedExternalNodes(0) << ", " << connectedExternalNodes(1) << "], ";
	s << "\"E\": " << E << ", ";
	s << "\"A\": "<< A << ", ";
    s << "\"Iz\": "<< I << ", ";
    s << "\"massperlength\": "<< rho << ", ";
    s << "\"crdTransformation\": \"" << theCoordTransf->getTag() << "\"}";
  }
}

int
BeamColumnwLHNMYS::displaySelf(Renderer &theViewer, int displayMode, float fact, const char **modes, int numMode)
{
  static Vector v1(3);
  static Vector v2(3);
  static Vector vp(3);

  theNodes[0]->getDisplayCrds(v1, fact);
  theNodes[1]->getDisplayCrds(v2, fact);

  float d1 = 0.0;
  float d2 = 0.0;
  float d3 = 0.0;

  int res = 0;

  if (displayMode > 0 && numMode == 0) {

    res += theViewer.drawLine(v1, v2, d1, d1, this->getTag(), 0);
    
  } else if (displayMode < 0) {
    
    theNodes[0]->getDisplayCrds(v1, 0.);
    theNodes[1]->getDisplayCrds(v2, 0.);
    
    // add eigenvector values
    int mode = displayMode  *  -1;

    const Matrix &eigen1 = theNodes[0]->getEigenvectors();
    const Matrix &eigen2 = theNodes[1]->getEigenvectors();
    if (eigen1.noCols() >= mode) {
      for (int i = 0; i < 2; i++) {
	v1(i) += eigen1(i,mode-1)*fact;
	v2(i) += eigen2(i,mode-1)*fact;    
      }    
    }

    res = theViewer.drawLine (v1, v2, 0.0, 0.0, this->getTag(), 0);
  }

  if (numMode > 0) {
    // calculate q for potential need below
    this->getResistingForce();
    vp = theCoordTransf->getBasicTrialDisp();
  }
  
  for (int i=0; i<numMode; i++) {

    const char *theMode = modes[i];
    if (strcmp(theMode, "axialForce") == 0) {
      d1 = q(0); 
      res +=theViewer.drawLine(v1, v2, d1, d1, this->getTag(), i);
      
    } else if (strcmp(theMode, "endMoments") == 0) {

      d1 = q(1);
      d2 = q(2);
      static Vector delta(3); delta = v2-v1; delta/=20.;
      res += theViewer.drawPoint(v1+delta, d1, this->getTag(), i);
      res += theViewer.drawPoint(v2-delta, d2, this->getTag(), i);

    } else if (strcmp(theMode, "localForces") == 0) {
      d1 = q(0);
      d2 = q(1);
      d3 = q(2);
      static Vector delta(3); delta = v2-v1; delta/=20;
      res += theViewer.drawPoint(v1+delta, d2, this->getTag(), i);
      res += theViewer.drawPoint(v2-delta, d3, this->getTag(), i);
      res +=theViewer.drawLine(v1, v2, d1, d1, this->getTag(), i);

    } else if (strcmp(theMode, "axialDeformation") == 0) {
      d1 = vp(0); 
      res +=theViewer.drawLine(v1, v2, d1, d1, this->getTag(), i);
      
    } else if (strcmp(theMode, "endRotations") == 0) {

      d1 = vp(1);
      d2 = vp(2);
      static Vector delta(3); delta = v2-v1; delta/=20.;
      res += theViewer.drawPoint(v1+delta, d1, this->getTag(), i);
      res += theViewer.drawPoint(v2-delta, d2, this->getTag(), i);

    } else if (strcmp(theMode, "localDeformations") == 0) {
      d1 = vp(0);
      d2 = vp(1);
      d3 = vp(2);
      static Vector delta(3); delta = v2-v1; delta/=20;
      res += theViewer.drawPoint(v1+delta, d2, this->getTag(), i);
      res += theViewer.drawPoint(v2-delta, d3, this->getTag(), i);
      res +=theViewer.drawLine(v1, v2, d1, d1, this->getTag(), i);

    } else if (strcmp(theMode, "plasticDeformations") == 0) {
      d1 = 0.;
      d2 = 0.;
      d3 = 0.;
      static Vector delta(3); delta = v2-v1; delta/=20;
      res += theViewer.drawPoint(v1+delta, d2, this->getTag(), i);
      res += theViewer.drawPoint(v2-delta, d3, this->getTag(), i);
      res +=theViewer.drawLine(v1, v2, d1, d1, this->getTag(), i);
    }

  }    

  return res;
}

Response*
BeamColumnwLHNMYS::setResponse(const char **argv, int argc, OPS_Stream &output)
{

  Response *theResponse = 0;

  output.tag("ElementOutput");
  output.attr("eleType","BeamColumnwLHNMYS");
  output.attr("eleTag",this->getTag());
  output.attr("node1",connectedExternalNodes[0]);
  output.attr("node2",connectedExternalNodes[1]);

    // global forces
  if (strcmp(argv[0],"force") == 0 || strcmp(argv[0],"forces") == 0 ||
      strcmp(argv[0],"globalForce") == 0 || strcmp(argv[0],"globalForces") == 0) {

    output.tag("ResponseType","Px_1");
    output.tag("ResponseType","Py_1");
    output.tag("ResponseType","Mz_1");
    output.tag("ResponseType","Px_2");
    output.tag("ResponseType","Py_2");
    output.tag("ResponseType","Mz_2");

    theResponse =  new ElementResponse(this, 2, P);
  
  // local forces
  }    else if (strcmp(argv[0],"localForce") == 0 || strcmp(argv[0],"localForces") == 0) {

    output.tag("ResponseType","N_1");
    output.tag("ResponseType","V_1");
    output.tag("ResponseType","M_1");
    output.tag("ResponseType","N_2");
    output.tag("ResponseType","V_2");
    output.tag("ResponseType","M_2");
    
    theResponse = new ElementResponse(this, 3, P);

  // basic forces
  }    else if (strcmp(argv[0],"basicForce") == 0 || strcmp(argv[0],"basicForces") == 0) {

    output.tag("ResponseType","N");
    output.tag("ResponseType","M_1");
    output.tag("ResponseType","M_2");
    
    theResponse = new ElementResponse(this, 4, Vector(3));

    // deformations
  }  else if (strcmp(argv[0],"deformations") == 0 ||
	      strcmp(argv[0],"basicDeformations") == 0) {
    
    output.tag("ResponseType","eps");
    output.tag("ResponseType","theta1");
    output.tag("ResponseType","theta2");
    theResponse = new ElementResponse(this, 5, Vector(3));
  
  // chord rotation -
  } else if (strcmp(argv[0],"chordRotation") == 0 || strcmp(argv[0],"chordDeformation") == 0 
	     || strcmp(argv[0],"basicDeformation") == 0) {

    output.tag("ResponseType","eps");
    output.tag("ResponseType","theta1");
    output.tag("ResponseType","theta2");

    theResponse =  new ElementResponse(this, 5, Vector(3));
  
    // plastic rotation -
  } else if (strcmp(argv[0],"plasticRotation") == 0 || strcmp(argv[0],"plasticDeformation") == 0) {
    
    output.tag("ResponseType","epsP");
    output.tag("ResponseType","thetaP_1");
    output.tag("ResponseType","thetaP_2");
    
    theResponse =  new ElementResponse(this, 6, Vector(3));

}
  output.endTag(); // ElementOutput
  
  return theResponse;
}

int
BeamColumnwLHNMYS::getResponse (int responseID, Information &eleInfo)
{
  double N, M1, M2, V;
  double L = theCoordTransf->getInitialLength();
  this->getResistingForce();

  switch (responseID) {
  case 1: // stiffness
    return eleInfo.setMatrix(this->getTangentStiff());
    
  case 2: // global forces
    return eleInfo.setVector(this->getResistingForce());
    
  case 3: // local forces
    // Axial
    N = q(0);
    P(3) =  N;
    P(0) = -N+p0[0];
    // Moment
    M1 = q(1);
    M2 = q(2);
    P(2) = M1;
    P(5) = M2;
    // Shear
    V = (M1+M2)/L;
    P(1) =  V+p0[1];
    P(4) = -V+p0[2];
    return eleInfo.setVector(P);
    
  case 4: // basic forces
    return eleInfo.setVector(q);

  case 5: // basic deformations
    return eleInfo.setVector(theCoordTransf->getBasicTrialDisp());
          
  case 6: // plastic deformations
    return eleInfo.setVector(vppres);
          
  default:
    return -1;
  }
}

int
BeamColumnwLHNMYS::setParameter(const char **argv, int argc, Parameter &param)
{
  if (argc < 1)
    return -1;

  // E of the beam interior
  if (strcmp(argv[0],"E") == 0)
    return param.addObject(1, this);

  // A of the beam interior
  if (strcmp(argv[0],"A") == 0)
    return param.addObject(2, this);
  
  // I of the beam interior
  if (strcmp(argv[0],"I") == 0)
    return param.addObject(3, this);
  
  return -1;
}

int
BeamColumnwLHNMYS::updateParameter (int parameterID, Information &info)
{
	switch (parameterID) {
	case -1:
		return -1;
	case 1:
		E = info.theDouble;
		return 0;
	case 2:
		A = info.theDouble;
		return 0;
	case 3:
		I = info.theDouble;
		return 0;
	default:
		return -1;
	}
}

