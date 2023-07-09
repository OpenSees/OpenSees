/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 1999, The Regents of the University of California    **
** All Rights Reserved.                                               **
**                                                                    **
** See file 'COPYRIGHT'  in main directory for information on usage   **
** and redistribution of OpenSees, and for a DISCLAIMER OF ALL        **
** WARRANTIES.                                                        **
**                                                                    **
** UpdatedLagrangianBeam2D.cpp: implementation of the                 **
**                             UpdatedLagrangianBeam2D class          **
** Developed by:                                                      **
**    Rohit Kaul       (rkaul@stanford.edu)                           **
**    Greg Deierlein   (ggd@stanford.edu)                             **
**                                                                    **
**           John A. Blume Earthquake Engineering Center              **
**                    Stanford University                             **
** ****************************************************************** **/

#include <Domain.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Renderer.h>
#include <Information.h>
#include <ElementResponse.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "UpdatedLagrangianBeam2D.h"

#define	  _debug  0
#define   _Kdebug 0

Matrix UpdatedLagrangianBeam2D::K(6,6);  
Matrix UpdatedLagrangianBeam2D::Kg(6,6);
Matrix UpdatedLagrangianBeam2D::Kt(6,6);
Matrix UpdatedLagrangianBeam2D::M(6,6);
Matrix UpdatedLagrangianBeam2D::D(6,6);
Matrix UpdatedLagrangianBeam2D::T(6,6);
Vector UpdatedLagrangianBeam2D::force(6);
Vector UpdatedLagrangianBeam2D::disp(6);
Vector UpdatedLagrangianBeam2D::ZeroVector(6);
Matrix UpdatedLagrangianBeam2D::ZeroMatrix(6,6);
Vector UpdatedLagrangianBeam2D::end1IncrDisp(3);
Vector UpdatedLagrangianBeam2D::end2IncrDisp(3);

Node *UpdatedLagrangianBeam2D::theNodes[2];

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

UpdatedLagrangianBeam2D::UpdatedLagrangianBeam2D(int classTag)
  :Element(0, classTag), isLinear(true), L(0.0), sn(0.0), cs(0),
   connectedExternalNodes(2), load(6),
   end1Ptr(0), end2Ptr(0),  eleForce(6), eleForce_hist(6), 
   nodeRecord(0), dofRecord(0), m_Iter(-1), Ki(0)
{
	numDof  = 6;
	massDof = -1;// assumes no lumped mass
    // will be overwritten by the sub-class massDof
}

UpdatedLagrangianBeam2D::UpdatedLagrangianBeam2D(int tag, int classTag, int nd1, int nd2, bool islinear)
  :Element(tag, classTag), isLinear(islinear), L(0.0), sn(0.0), cs(0),
   connectedExternalNodes(2), load(6),
   end1Ptr(0), end2Ptr(0),  eleForce(6), eleForce_hist(6), 
   nodeRecord(0), dofRecord(0), m_Iter(-1), Ki(0)
{
	connectedExternalNodes(0) = nd1;
	connectedExternalNodes(1) = nd2;
	numDof  = 6;
	massDof = -1;// assumes no lumped mass
    // will be overwritten by the sub-class massDof
}

UpdatedLagrangianBeam2D::~UpdatedLagrangianBeam2D()
{
  if (Ki != 0)
    delete Ki;
}

//////////////////////////////////////////////////////////////////////
// Add the 2D element to the domain
//////////////////////////////////////////////////////////////////////

void UpdatedLagrangianBeam2D::setDomain(Domain *theDomain)
{
	// check Domain is not null - invoked when object removed from a domain
    if (theDomain == 0) {
	end1Ptr = 0;
	end2Ptr = 0;
	L = 0;
    }
    // first set the node pointers
    //cout <<  "Element " << getTag() << ", Nodes = " << connectedExternalNodes; cin.get();

    int Nd1 = connectedExternalNodes(0);
    int Nd2 = connectedExternalNodes(1);
    end1Ptr = theDomain->getNode(Nd1);
    end2Ptr = theDomain->getNode(Nd2);	
    if (end1Ptr == 0) {
	opserr << "WARNING (W_C_10) - UpdatedLagrangianBeam2D::setDomain(..) [" << getTag() << "]\n";
	opserr << Nd1 << "Nd1 does not exist in model for element \n" <<" Tag = " << getTag();
	return;
    }
    if (end2Ptr == 0) {
	opserr << "WARNING (W_C_20) - UpdatedLagrangianBeam2D::setDomain(..) [" << getTag() << "]\n";
	opserr << Nd2 << "Nd2 does not exist in model for element\n" << " Tag = " << getTag();
	return;
    }

    // now verify the number of dof at node ends
    int dofNd1 = end1Ptr->getNumberDOF();
    int dofNd2 = end2Ptr->getNumberDOF();
    if (dofNd1 != 3 && dofNd2 != 3) {
	opserr << "WARNING (W_C_30) - UpdatedLagrangianBeam2D::setDomain() [" << getTag() <<"]\n";
	opserr << "node and/or node " << Nd1 << Nd2 << " have/has incorrect number ";
	opserr << "of dof's at end for element\n " << *this;
	return;
    }

    // call the base class method
    this->DomainComponent::setDomain(theDomain);

    // determine length and direction cosines
	{
		double dx,dy;
		const Vector &end1Crd = end1Ptr->getCrds();
		const Vector &end2Crd = end2Ptr->getCrds();
    
		dx = end2Crd(0)-end1Crd(0);
		dy = end2Crd(1)-end1Crd(1);	
  
		L = sqrt(dx*dx + dy*dy);
		L_hist = L;

		if(L==0)
		{
			opserr << "WARNING UpdatedLagrangianBeam2D::setDomain(): zero length\n";
			return;	
		}
		cs = dx/L;
		sn = dy/L;
		cs_hist = dx/L;
		sn_hist = dy/L;
	}


}//setDomain


int UpdatedLagrangianBeam2D::update()
{
	return 0;
}


int UpdatedLagrangianBeam2D::commitState()
{
  int success = 0 ;

  // call element commitState to do any base class stuff
  if ((success = this->Element::commitState()) != 0) {
    opserr << "UpdatedLagrangianBeam2D::commitState () - failed in base class";
  }    

#ifdef _G3DEBUG
	opserr << m_Iter << " "; // cin.get();
#endif
		m_Iter = 0;
	
	if(_debug || _Kdebug)
	{
	    opserr << "\n Beam N0: "<< this->getTag()
	    << " ----- Inside Commit State -----\n";
	    //cin.get();
	}

    if(!isLinear)
    {
        this->updateState();
        cs_hist = cs;
        sn_hist = sn;
		L_hist = L;
    }

	// save the prev. step element forces
	eleForce_hist = eleForce;
	return success;
}

void UpdatedLagrangianBeam2D::updateState()
{
	// Update the direction cosines to the new trial state
	/////////////////////////////////////////////////////////////////////
	//return;
	/////////////////////////////////////////////////////////////////////

	const Vector &end1Crd = end1Ptr->getCrds();
	const Vector &end2Crd = end2Ptr->getCrds();
	const Vector &end1Disp  = end1Ptr->getTrialDisp();
    const Vector &end2Disp  = end2Ptr->getTrialDisp();

	double x1 = end1Crd(0) + end1Disp(0);
	double y1 = end1Crd(1) + end1Disp(1);
	double x2 = end2Crd(0) + end2Disp(0);
	double y2 = end2Crd(1) + end2Disp(1);

	double dx = x2 - x1;
	double dy = y2 - y1;

	L = sqrt(dx*dx + dy*dy);
	if(L==0 )
	{
		opserr << "WARNING (W_B_40) - UpdatedLagrangianBeam2D::updateState() [" << getTag() << "\n";
		opserr << "L = 0\n";
		return;
	}
	cs = dx/L;
	sn = dy/L;
}


int UpdatedLagrangianBeam2D::revertToLastCommit()
{
	cs = cs_hist;
	sn = sn_hist;
	L = L_hist;
	eleForce = eleForce_hist;
	return 0;
}

//////////////////////////////////////////////////////////////////////
// protected methods for local stiffness and force formulation
//////////////////////////////////////////////////////////////////////

// set default for 2D beam-column elements

void UpdatedLagrangianBeam2D::addInternalGeomStiff(Matrix &K)
{
   if (isLinear)
    	return;
    	
 double P = eleForce_hist(3);  //Last committed
 double l = L_hist;

    K(0,0) += P/l;
	K(0,3) += -P/l;
	K(3,0) += -P/l;
	K(3,3) += P/l;

	K(1,1) += 1.2*P/l;
	K(1,4) += -1.2*P/l;
	K(4,1) += -1.2*P/l;
	K(4,4) += 1.2*P/l;

	K(1,2) += P/10;
	K(1,5) += P/10;
	K(2,1) += P/10;
	K(5,1) += P/10;

	K(2,2) += 2*P*l/15;
	K(2,5) += -P*l/30;
	K(5,2) += -P*l/30;
	K(5,5) += 2*P*l/15;


	K(2,4) += -P/10;
	K(4,2) += -P/10;
	K(4,5) += -P/10;
	K(5,4) += -P/10;

/////////////////////////////old ///////
/*    K(1, 1) +=  0.2*P/l;
    K(4, 4) +=  0.2*P/l;
	K(1, 4) += -0.2*P/l;
    K(4, 1) += -0.2*P/l;
	K(1, 2) += P/10;
    K(1, 5) += P/10;
    K(2, 1) += P/10;
    K(5, 1) += P/10;
	K(2, 4) += -P/10;
    K(4, 2) += -P/10;
    K(4, 5) += -P/10;
    K(5, 4) += -P/10;
	K(2, 2) += 2*P*l/15;
    K(5, 5) += 2*P*l/15;
	K(2, 5) += -P*l/30;
    K(5, 2) += -P*l/30;
*/

}

void UpdatedLagrangianBeam2D::addExternalGeomStiff(Matrix &K)
{
//     if (isLinear)
    	return;
// everything is included above for rf2
 
double P = eleForce_hist(3);  //Last committed
double V = eleForce_hist(1);
double l = L_hist;

    K(0, 1) += V/l;
    K(1, 0) += V/l;
	K(0, 4) += -V/l;
    K(4, 0) += -V/l;
	K(1, 1) += P/l;
    K(4, 4) += P/l;
	K(1, 3) += -V/l;
    K(3, 1) += -V/l;
	K(1, 4) += -P/l;
    K(4, 1) += -P/l;
	K(3, 4) += P/l;
    K(4, 3) += P/l;

}


//////////////////////////////////////////////////////////////////////
// Public methods called, taken care of for 2D element subclasses
//////////////////////////////////////////////////////////////////////

int UpdatedLagrangianBeam2D::getNumExternalNodes(void) const
{
    return 2;
}

const ID &UpdatedLagrangianBeam2D::getExternalNodes(void) 
{
    return connectedExternalNodes;
}

Node **
UpdatedLagrangianBeam2D::getNodePtrs(void) 
{
  theNodes[0] = end1Ptr;
  theNodes[1] = end2Ptr;
  
  return theNodes;
}

int UpdatedLagrangianBeam2D::getNumDOF(void) 
{
    return numDof;
}

const Matrix &UpdatedLagrangianBeam2D::getTangentStiff(void)
{
    // Get the local elastic stiffness matrix, store in Kt
    getLocalStiff(Kt);

    // Add internal geometric stiffness matrix
    addInternalGeomStiff(Kt);

    // Add external geometric stiffness matrix
    addExternalGeomStiff(Kt);

	if(_Kdebug){
	    opserr << "UpdatedLagrangianBeam2D::getTangentStiff(void) tag = " << getTag() << "\n";
	    opserr << Kt;
	}

    //Kt = T^(Kt*T);
    transformToGlobal(Kt);

#ifdef _G3DEBUG
	for(int i=0; i< 6; i++)
	{
		double aii = Kt(i, i);
		if(aii < 1e-6)
		{
			opserr << " WARNING (W_B_50) - UpdatedLagrangianBeam2D::getTangentStiff(..) [" << getTag() << "]\n";
			opserr << " aii = " << aii << ", i = " << i << "\n";
		}
	}
#endif

    return Kt;
}

const Matrix &UpdatedLagrangianBeam2D::getInitialStiff(void)
{
  if (Ki == 0)
    Ki = new Matrix(this->getTangentStiff());
  
  return *Ki;
}


const Matrix &UpdatedLagrangianBeam2D::getMass(void)
{
  if(massDof==0)
    return ZeroMatrix;

  getLocalMass(M);
  transformToGlobal(M);

  return M;
}

//////////////////////////////////////////////////////////////////////
// methods for applying and returning loads
//////////////////////////////////////////////////////////////////////
// Add uniformly varying load to the element
Vector &UpdatedLagrangianBeam2D::getUVLoadVector(double q1, double q2)
{
 load(0) = 0;
 load(1) = (7*q1 + 3*q2)*(L/20);
 load(2) = (3*q1 + 2*q2)*(L*L/60);
 load(3) = 0;
 load(4) = (3*q1 + 7*q2)*(L/20);
 load(5) = -(2*q1 + 3*q2)*(L*L/60);
 return load;
}

void UpdatedLagrangianBeam2D::zeroLoad(void)
{
    load.Zero();
}

int UpdatedLagrangianBeam2D::addLoad(const Vector &moreLoad)
{
    if (moreLoad.Size() != numDof) {
	opserr << "WARNING (W_C_80) - UpdatedLagrangianBeam2D::addLoad(..) [" << getTag() << "]\n";
	opserr << "vector not of correct size\n";
	return -1;
    }
    load += moreLoad;
    return 0;
}

void UpdatedLagrangianBeam2D::getTrialLocalForce(Vector &lforce)
{
	// Get the local elastic stiffness matrix, store in Kt
    getLocalStiff(Kt);

    // Add internal geometric stiffness matrix
    addInternalGeomStiff(Kt);

    // Get incremental local displacements  trial-conv.
	
	if(isLinear)
		getIncrLocalDisp(disp);
	else
		getIncrNaturalDisp(disp);
	
	//~ cout << disp << endln;
/* 
//////////////////////////////////////////
// using the natural deformation approach
//////////////////////////////////////////

	double l  = L_hist;
	double ua = disp(0);
	double ub = disp(3);
	double va = disp(1);
	double vb = disp(4);
	double ra = disp(2);
	double rb = disp(5);

	double un =   (ub - ua) 
	            + ( (ub - ua)*(ub - ua) + (vb - va)*(vb - va) )/(2*l);

	double rr  = atan( (vb - va)/(l + ub - ua) );
	double ran = ra - rr;
	double rbn = rb - rr;

	disp(0) = 0;
	disp(1) = 0;
	disp(2) = ran;
	disp(3) = un;
	disp(4) = 0;
	disp(5) = rbn;
///////////////////////////////////////////////
// end natural defo
///////////////////////////////////////////////
*/
    // Compute local incremental force
    force = Kt*disp;

    // Compute total local force
    lforce = eleForce_hist + force;

}

const Vector &UpdatedLagrangianBeam2D::getResistingForce()
{
    // check for quick return
    if (L == 0)
	    return ZeroVector;
	
	m_Iter++;

	if(!isLinear)
		this->updateState();
	
	this->getTrialLocalForce(eleForce);

///////////////////////////////////////////////
// start elementary
///////////////////////////////////////////////

/*
double f0 = eleForce(0);
double f1 = eleForce(1);
double f2 = eleForce(2);
double f3 = eleForce(3);
double f4 = eleForce(4);
double f5 = eleForce(5);

    eleForce(0) = (cs*cs_hist+sn*sn_hist)*f0+(-cs*sn_hist+sn*cs_hist)*f1;
 	eleForce(1) = (-sn*cs_hist+cs*sn_hist)*f0+(cs*cs_hist+sn*sn_hist)*f1;
  	eleForce(2) = f2;
 	eleForce(3) = (cs*cs_hist+sn*sn_hist)*f3+(-cs*sn_hist+sn*cs_hist)*f4;
 	eleForce(4) = (-sn*cs_hist+cs*sn_hist)*f3+(cs*cs_hist+sn*sn_hist)*f4;
    eleForce(5) = f5;
*/
//////////////////////////////////////////////////
// end elementary
//////////////////////////////////////////////////

	double cos = cs;
	double sin = sn;
    // determine the ele end forces in global coords - want -F into rForce
    force(0) =  cos*eleForce(0) - sin*eleForce(1);
    force(1) =  sin*eleForce(0) + cos*eleForce(1);
    force(2) =  eleForce(2);
    force(3) =  cos*eleForce(3) - sin*eleForce(4);
    force(4) =  sin*eleForce(3) + cos*eleForce(4);
    force(5) =  eleForce(5);

    if(_debug)
      { opserr << "Global forces:\n " << force; 
      }

    return force;
}

const Vector &
UpdatedLagrangianBeam2D::getResistingForceIncInertia()
{
    // check for quick return
    if (L == 0)
        return ZeroVector;

    // form the stiffness matrix - we will do a P-kU-MA
    force = this->getResistingForce(); // now we have P-Ku

    // massDof <  0 for distributed mass
    // massDof =  0 if mass is specified at nodes
    // massDof >  0 for lumped mass


    // determine -Ma
    // optimize for lumped mass
    if (massDof != 0) {
      if(massDof > 0) {
	const Vector &end1Accel = end1Ptr->getTrialAccel();
	const Vector &end2Accel = end2Ptr->getTrialAccel();
	force(0) -= massDof*end1Accel(0);
	force(1) -= massDof*end1Accel(1);
	force(3) -= massDof*end2Accel(0);
	force(4) -= massDof*end2Accel(1);
      } else if(massDof < 0) {
	M = this->getMass();
	
	const Vector &end1Accel = end1Ptr->getTrialAccel();
	const Vector &end2Accel = end2Ptr->getTrialAccel();
	Vector Accel(6), f(6);
	int i=0;
	for(i=0; i<3; i++)
	  {
	    Accel(i)   = end1Accel(i);
	    Accel(i+3) = end2Accel(i);
	  }
	f = M*Accel;
	for(i=0; i<6; i++) force(i) -= f(i);
      }

      if (alphaM != 0.0 || betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0)
	force += this->getRayleighDampingForces();
    } else {

      if (betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0)
	force += this->getRayleighDampingForces();
    }    
    
    return force;
}

void UpdatedLagrangianBeam2D::transformToGlobal(Matrix &K)
{
double k00=K(0,0), k01=K(0,1), k02=K(0,2), k03=K(0,3), k04=K(0,4), k05=K(0,5);
double k11=K(1,1), k12=K(1,2), k13=K(1,3), k14=K(1,4), k15=K(1,5);
double k22=K(2,2), k23=K(2,3), k24=K(2,4), k25=K(2,5);
double k33=K(3,3), k34=K(3,4), k35=K(3,5);
double k44=K(4,4), k45=K(4,5);
double k55=K(5,5);

 double cos = cs;
 double sin = sn;

    K(0,0) = (cos*k00-sin*k01)*cos-(cos*k01-sin*k11)*sin;
    K(0,1) = (cos*k00-sin*k01)*sin+(cos*k01-sin*k11)*cos;
    K(0,2) =  cos*k02-sin*k12;
    K(0,3) = (cos*k03-sin*k13)*cos-(cos*k04-sin*k14)*sin;
    K(0,4) = (cos*k03-sin*k13)*sin+(cos*k04-sin*k14)*cos;
    K(0,5) =  cos*k05-sin*k15;

    K(1,1) = (sin*k00+cos*k01)*sin+(sin*k01+cos*k11)*cos;
    K(1,2) =  sin*k02+cos*k12;
    K(1,3) = (sin*k03+cos*k13)*cos-(sin*k04+cos*k14)*sin;
    K(1,4) = (sin*k03+cos*k13)*sin+(sin*k04+cos*k14)*cos;
    K(1,5) =  sin*k05+cos*k15;

    K(2,2) =  k22;
    K(2,3) =  k23*cos-k24*sin;
    K(2,4) =  k23*sin+k24*cos;
    K(2,5) =  k25;

    K(3,3) = (cos*k33-sin*k34)*cos-(cos*k34-sin*k44)*sin;
    K(3,4) = (cos*k33-sin*k34)*sin+(cos*k34-sin*k44)*cos;
    K(3,5) =  cos*k35-sin*k45;

    K(4,4) = (sin*k33+cos*k34)*sin+(sin*k34+cos*k44)*cos;
    K(4,5) =  sin*k35+cos*k45;

    K(5,5) =  k55;

    for(int i=1; i<6; i++)
    {
        for(int j=0; j<i; j++)
        {
            K(i,j) = K(j,i);
        }
    }

}

//////////////////////////////////////////////////////////////////////
// methods for getting local displacements
//////////////////////////////////////////////////////////////////////
void UpdatedLagrangianBeam2D::getTrialNaturalDisp(Vector &nDisp)
{
    // getIncrLocalDisp(disp);
	getTrialLocalDisp(disp);
//////////////////////////////////////////
// using the natural deformation approach
//////////////////////////////////////////

	double l  = L_hist;
	double ua = disp(0);
	double ub = disp(3);
	double va = disp(1);
	double vb = disp(4);
	double ra = disp(2);
	double rb = disp(5);

	double un =   (ub - ua) 
	            + ( (ub - ua)*(ub - ua) + (vb - va)*(vb - va) )/(2*l);

	double rr  = atan( (vb - va)/(l + ub - ua) );
	double ran = ra - rr;
	double rbn = rb - rr;

	nDisp(0) = 0;
	nDisp(1) = 0;
	nDisp(2) = ran;
	nDisp(3) = un;
	nDisp(4) = 0;
	nDisp(5) = rbn;
}

void UpdatedLagrangianBeam2D::getIncrNaturalDisp(Vector &nDisp)
{
    getIncrLocalDisp(disp);

//////////////////////////////////////////
// using the natural deformation approach
//////////////////////////////////////////

	double l  = L_hist;
	double ua = disp(0);
	double ub = disp(3);
	double va = disp(1);
	double vb = disp(4);
	double ra = disp(2);
	double rb = disp(5);

	double un =   (ub - ua) 
	            + ( (ub - ua)*(ub - ua) + (vb - va)*(vb - va) )/(2*l);

	double rr  = atan( (vb - va)/(l + ub - ua) );
	double ran = ra - rr;
	double rbn = rb - rr;

	nDisp(0) = 0;
	nDisp(1) = 0;
	nDisp(2) = ran;
	nDisp(3) = un;
	nDisp(4) = 0;
	nDisp(5) = rbn;


}


void UpdatedLagrangianBeam2D::getIncrLocalDisp(Vector &lDisp)
{
    if (L == 0.0)
	return;

    const Vector &end1TrialDisp  = end1Ptr->getTrialDisp();
    const Vector &end2TrialDisp  = end2Ptr->getTrialDisp();    
	const Vector &end1CommitDisp = end1Ptr->getDisp();
	const Vector &end2CommitDisp = end2Ptr->getDisp();
	

	for(int i=0; i<3; i++)
	{
		end1IncrDisp(i) = end1TrialDisp(i) - end1CommitDisp(i);
		end2IncrDisp(i) = end2TrialDisp(i) - end2CommitDisp(i);
	}

    lDisp(0) = cs_hist * end1IncrDisp(0) + sn_hist * end1IncrDisp(1);
    lDisp(1) = cs_hist * end1IncrDisp(1) - sn_hist * end1IncrDisp(0);
    lDisp(2) = end1IncrDisp(2);
    lDisp(3) = cs_hist * end2IncrDisp(0) + sn_hist * end2IncrDisp(1);
    lDisp(4) = cs_hist * end2IncrDisp(1) - sn_hist * end2IncrDisp(0);
    lDisp(5) = end2IncrDisp(2);

    return;
}

void UpdatedLagrangianBeam2D::getTrialLocalDisp(Vector &lDisp)
{   
    if(L == 0.0)
	return;

    const Vector &end1Disp  = end1Ptr->getTrialDisp();
    const Vector &end2Disp  = end2Ptr->getTrialDisp();

    lDisp(0) = cs * end1Disp(0) + sn * end1Disp(1);
    lDisp(1) = cs * end1Disp(1) - sn * end1Disp(0);
    lDisp(2) = end1Disp(2);    
    lDisp(3) = cs * end2Disp(0) + sn * end2Disp(1);
    lDisp(4) = cs * end2Disp(1) - sn * end2Disp(0);
    lDisp(5) = end2Disp(2);  
	
    return;
}


void UpdatedLagrangianBeam2D::getConvLocalDisp(Vector &lDisp)
{   
    if (L == 0.0)
	return;
    
	const Vector &end1Disp = end1Ptr->getDisp();
	const Vector &end2Disp = end2Ptr->getDisp();

    lDisp(0) = cs_hist * end1Disp(0) + sn_hist * end1Disp(1);
    lDisp(1) = cs_hist * end1Disp(1) - sn_hist * end1Disp(0);
    lDisp(2) = end1Disp(2);    
    lDisp(3) = cs_hist * end2Disp(0) + sn_hist * end2Disp(1);
    lDisp(4) = cs_hist * end2Disp(1) - sn_hist * end2Disp(0);
    lDisp(5) = end2Disp(2);    
	
    return;
}


//////////////////////////////////////////////////////////////////////
// methods for display/recorders... may be overridden if required
//////////////////////////////////////////////////////////////////////

int UpdatedLagrangianBeam2D::displaySelf(Renderer &theViewer, int displayMode, float fact, const char **modes, int numMode)
{
    // first determine the two end points of the element based on
    // the display factor (a measure of the distorted image)
    // store this information in 2 3d vectors v1 and v2

    const Vector &end1Crd = end1Ptr->getCrds();
    const Vector &end2Crd = end2Ptr->getCrds();	
    const Vector &end1Disp = end1Ptr->getDisp();
    const Vector &end2Disp = end2Ptr->getDisp();    
    Vector rgb(3);
	rgb(0) = 0; rgb(1) = 0; rgb(2) = 1;

	Vector v1(3);
    Vector v2(3);

    for (int i=0; i<2; i++) {
	v1(i) = end1Crd(i)+end1Disp(i)*fact;
	v2(i) = end2Crd(i)+end2Disp(i)*fact;    	
    }	
	v1(2) = 0;
	//v1(2) = end1Disp(2)*fact;
	v2(2) = 0;
	//v2(2) = end2Disp(2)*fact;
	
	//opserr << v1 << v2;
	if (displayMode == 1) theViewer.drawLine(v1, v2, rgb, rgb);	
	// cout << "line drawn << " << getTag() << endln;
	// cin.get();
	return 0;
    
    if (displayMode == 2) // axial force
	return theViewer.drawLine(v1, v2, eleForce(0), eleForce(3));
    else if (displayMode == 3) // shear forces
	return theViewer.drawLine(v1, v2, eleForce(1), eleForce(4));
    else if (displayMode == 4) // bending moments
	return theViewer.drawLine(v1, v2, eleForce(2), eleForce(5));
	
	return 0;
}

Response* UpdatedLagrangianBeam2D::setResponse(const char **argv, int argc)
{
    // force (axialForce)
    if ((strcmp(argv[0],"force") == 0) ||
	(strcmp(argv[0],"forces") == 0) ||
	(strcmp(argv[0],"localForce") == 0))
	{
		return new ElementResponse(this, 1, Vector(6));
    }

	else if((strcmp(argv[0],"forceDisp")==0))
	{
		if(strcmp(argv[1],"1") == 0) nodeRecord = 1;
		else nodeRecord = 2;

		if(strcmp(argv[2],"0") == 0) dofRecord = 0;
		if(strcmp(argv[2],"1") == 0) dofRecord = 1;
		if(strcmp(argv[2],"2") == 0) dofRecord = 2;

		return  new ElementResponse(this, 4, Vector(7));

	}

	else if ((strcmp(argv[0],"globalForce") == 0))
	{
		return  new ElementResponse(this, 5, Vector(6));
    }

    else if ((strcmp(argv[0],"disp") == 0) ||
	(strcmp(argv[0],"displacements") == 0) ||
	(strcmp(argv[0],"displacement") == 0))
	{
		return  new ElementResponse(this, 2, Vector(6));
    }

    // tangent stiffness matrix
    else if (strcmp(argv[0],"stiffness") == 0) {
	return new ElementResponse(this, 3, Matrix(6,6));
    }

    // a material quantity: needs to be implemented in subclasses

    // otherwise response quantity is unknown for the UpdatedLagrangianBeam2D class
    /*else
    {
    	opserr << "WARNING (W_C_90) - UpdatedLagrangianBeam2D::setResponse(..) [" << getTag() << "]\n";
		opserr << "unknown response quantity\n";
		return 0;
	}*/
	
	return 0;
}

int UpdatedLagrangianBeam2D::getResponse(int responseID, Information &eleInformation)
{
  switch (responseID) {
    case -1:
      return -1;

    case 1:
		if(eleInformation.theVector!=0)
		{
			*(eleInformation.theVector) = eleForce;
		}

      return 0;

    case 2:
		if(eleInformation.theVector!=0)
		{
			this->getTrialLocalDisp(disp);
			*(eleInformation.theVector) = disp;
		}
      return 0;

    case 3:
      if (eleInformation.theMatrix != 0)
	  {
		*(eleInformation.theMatrix) = this->getTangentStiff();
	  }
      return 0;

	case 4:
		if(eleInformation.theVector!=0)
		{
			Vector disp(3);
			if(nodeRecord==1)	disp = end1Ptr->getDisp();
			else				disp = end2Ptr->getDisp();
			
			Vector temp(7);
			temp(0) = disp(dofRecord);
			for(int i=1; i < 7; i++) temp(i) = eleForce(i-1);
			
			eleInformation.theVector->addVector(0, temp, 1);
		}
		return 0;
	case 5:
		if(eleInformation.theVector!=0)
		{
			double cos = cs;
			double sin = sn;
			// determine the ele end forces in global coords - want -F into rForce
			force(0) =  cos*eleForce(0) - sin*eleForce(1);
			force(1) =  sin*eleForce(0) + cos*eleForce(1);
			force(2) =  eleForce(2);
			force(3) =  cos*eleForce(3) - sin*eleForce(4);
			force(4) =  sin*eleForce(3) + cos*eleForce(4);
			force(5) =  eleForce(5);
			
			*(eleInformation.theVector) = force;
		}
      return 0;



    default:      
	  return -1;
  }
}
