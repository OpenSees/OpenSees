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

// $Revision: 1.2
// $Date: 2010-11-10
// $Source: /OpenSees/SRC/material/nD/ContactMaterial3D.cpp,v $
                                                                        
// Written: Kathryn Petek
// Created: February 2004
// Modified: Chris McGann
//           November 2010 -> changes for incorporation into main source code
// Modified: Chris McGann
//           Jan 2011 -> added update for frictional state

// Description: This file contains the implementation for the ContactMaterial3D class.
//

#include <ContactMaterial3D.h>

#include <Information.h>
#include <MaterialResponse.h>
#include <Parameter.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <string.h>

#include <elementAPI.h>
#define OPS_Export
static int numContactMaterial3DMaterials = 0;
int ContactMaterial3D::mFrictFlag = 1;

OPS_Export void * OPS_ADD_RUNTIME_VPV(OPS_ContactMaterial3DMaterial)
{
  if (numContactMaterial3DMaterials == 0) {
    numContactMaterial3DMaterials++;
    opserr << "ContactMaterial3D nDmaterial - Written: K.Petek, P.Mackenzie-Helnwein, P.Arduino, U.Washington\n";
  }

  // Pointer to a uniaxial material that will be returned
  NDMaterial *theMaterial = 0;

  int numArgs = OPS_GetNumRemainingInputArgs();

  if (numArgs < 5) {
    opserr << "Want: nDMaterial ContactMaterial3D tag? mu? G? c? t?\n";
    return 0;  
  }
 
  int tag;
  double dData[4];

  int numData = 1;
  if (OPS_GetInt(&numData, &tag) != 0) {
    opserr << "WARNING invalid tag for  ContactMaterial3D material" << endln;
    return 0;
  }
  numData = 4;
  if (OPS_GetDouble(&numData, dData) != 0) {
    opserr << "WARNING invalid material data for nDMaterial ContactMaterial3D material  with tag: " << tag << endln;
    return 0;
  }

  theMaterial = new ContactMaterial3D(tag, dData[0], dData[1], dData[2], dData[3]);

  if (theMaterial == 0) {
    opserr << "WARNING ran out of memory for nDMaterial ContactMaterial3D material  with tag: " << tag << endln;
  }

  return theMaterial;
}


//full constructor
ContactMaterial3D::ContactMaterial3D (int tag, double mu, double Gmod, double c, double t)
 : NDMaterial(tag,ND_TAG_ContactMaterial3D),
   s_e_n(2),
   s_e_nplus1(2),
   r_nplus1(2),
   g(2,2),
   G(2,2),
   strain_vec(4),
   stress_vec(4),
   tangent_matrix(4,4)
{
#ifdef DEBUG
  opserr << "ContactMaterial3D::ContactMaterial3D(...)" << endln;
#endif

  frictionCoeff   = mu;
  mMu             = mu;
  stiffness       = Gmod;
  cohesion        = c;
  mCo             = c;
  tensileStrength = t;
  mTen            = t;
  
  this->zero();
}

//full constructor
ContactMaterial3D::ContactMaterial3D (const ContactMaterial3D &other)
  : NDMaterial(other.getTag(),ND_TAG_ContactMaterial3D),
    s_e_n(2),
    s_e_nplus1(2),
    r_nplus1(2),
    g(2,2),
    G(2,2),
    strain_vec(4),
    stress_vec(4),
    tangent_matrix(4,4)
{
#ifdef DEBUG
  opserr << "ContactMaterial3D::ContactMaterial3D(...)" << endln;
#endif

  this->frictionCoeff = other.frictionCoeff;
  this-> mMu = other. mMu;
  this->stiffness = other.stiffness;
  this->cohesion  = other.cohesion; 
  this->mCo = other.mCo;
  this->tensileStrength = other.tensileStrength;
  this->mTen = other.mTen;

  this->zero();
}
   
//null constructor
ContactMaterial3D::ContactMaterial3D ()
  : NDMaterial(0, ND_TAG_ContactMaterial3D),
    s_e_n(2),
    s_e_nplus1(2),
    r_nplus1(2),
    g(2,2),
    G(2,2),
    strain_vec(4),
    stress_vec(4),
    tangent_matrix(4,4)
{
#ifdef DEBUG
  opserr << "ContactMaterial3D::ContactMaterial3D()" << endln;
#endif

  frictionCoeff   = 0.0;
  mMu             = 0.0;
  stiffness       = 1.0;
  cohesion        = 0.0;
  mCo             = 0.0;
  tensileStrength = 0.0;
  mTen            = 0.0;
}

//destructor
ContactMaterial3D::~ContactMaterial3D ()
{
#ifdef DEBUG
        opserr << "ContactMaterial3D::~ContactMaterial3D()" << endln;
#endif
}

//zero internal variables
void ContactMaterial3D::zero( )
{
#ifdef DEBUG
  opserr << "ContactMaterial3D::zero( )" << endln;
#endif
  s_e_n.Zero();         // elastic slip from previous increment
  s_e_nplus1.Zero();    // elastic slip after current increment
       
  r_nplus1.Zero();      // direction of plastic slip

  inSlip = false;    
  mFlag = 1;

  stress_vec.Zero();
  strain_vec.Zero();
  tangent_matrix.Zero();
}


int ContactMaterial3D::setTrialStrain (const Vector &strain_from_element)
{
#ifdef DEBUG
  opserr << "ContactMaterial3D::setTrialStrain (const Vector &strain_from_element)" << endln;
#endif
  Vector t_s(2);          // tangential contact force
  double t_n;             // normal contact force
  double f_nplus1_trial;  // trial slip condition

  double gap;             // current gap
  Vector slip(2);         // incremental slip

  double t_s_norm;        // norm of tangential contact force
       
  strain_vec = strain_from_element;

  gap        = strain_vec(0);
  slip(0)    = strain_vec(1);
  slip(1)    = strain_vec(2);
  t_n        = strain_vec(3);    

  Vector zeroVec = slip;  
  zeroVec.Zero();

  // update frictional status
  this->UpdateFrictionalState();

  // trial state (elastic predictor step) -> assume sticking
  inSlip = false;
       
  s_e_nplus1 = (t_n > -tensileStrength) ?  s_e_n + slip : zeroVec; // ctv

  t_s        = stiffness * g * s_e_nplus1; // cov

  // Norm(s_e_nplus1) = sqrt( s_e_nplus1' * g * s_e_nplus1 )
  s_e_nplus1_norm = sqrt( s_e_nplus1(0) * g(0,0) * s_e_nplus1(0)
                          + s_e_nplus1(1) * g(1,0) * s_e_nplus1(0) * 2.0
                          + s_e_nplus1(1) * g(1,1) * s_e_nplus1(1) );

  // Norm(t_s) = sqrt( t_s' * g * t_s )
  //t_s_norm = sqrt( t_s(0) * G(0,0) * t_s(0)
  //                      + t_s(1) * G(1,0) * t_s(0) * 2.0
  //                      + t_s(1) * G(1,1) * t_s(1) );

  //Norm(t_s) = k*Norm(s_e_nplus1)  -- yields same result as above
  t_s_norm = stiffness * s_e_nplus1_norm;

  // slip condition
  f_nplus1_trial = t_s_norm - frictionCoeff*t_n - cohesion;

  // if (f_nplus1_trial > 0.0) {
  // if ( (f_nplus1_trial > 0.0) && (t_n > -cohesion/frictionCoeff) &&  (slip.Norm() > 1e-12) ) {
  if ( (f_nplus1_trial > 0.0) && (t_n > -tensileStrength) && (s_e_nplus1_norm > 1e-12) ) {
 
    // plastic corrector step -> sliding
    inSlip = true;

    gamma = f_nplus1_trial / stiffness * 0.999999999999 ;

    r_nplus1 = s_e_nplus1 / s_e_nplus1_norm; // ctv

    // s_p_nplus1 = s_p_n + gamma * r_nplus1
    // s_e_nplus1 = s_nplus1 - s_p_nplus1
    //        = (s_nplus1 - s_p_n) - gamma * r_nplus1
    //        = (s_n + slip - s_p_n) - gamma * r_nplus1
    //        = (s_e_n + slip) - gamma * r_nplus1
    //        = s_e_nplus1_trial - gamma * r_nplus1
    double scale = (1.0 - gamma/s_e_nplus1_norm);

    s_e_nplus1 = scale * s_e_nplus1; // ctv
    t_s        = scale * t_s;        // cov

  }

#ifdef DEBUG
  if (DEBUG_LEVEL > 1) {
    if (inSlip) {
      opserr << "   ** SLIDING (material)" << endln;
    } else {
      opserr << "   ** STICKING (material)" << endln;
	}
  }
#endif

  //update stress and strain values
  stress_vec(0) = t_n;
  stress_vec(1) = t_s(0);
  stress_vec(2) = t_s(1);
  stress_vec(3) = gap;
       
  return 0;
}


//unused trial strain functions
int ContactMaterial3D::setTrialStrain (const Vector &v, const Vector &r)
{
#ifdef DEBUG
  opserr << "ContactMaterial3D::setTrialStrain (const Vector &v, const Vector &r)" << endln;
#endif
  opserr << "YOU SHOULD NOT SEE THIS: ContactMaterial3D::setTrialStrain (const Vector &v, const Vector &r)" << endln;
  return this->setTrialStrain (v);
}


const Matrix & ContactMaterial3D::getTangent ()
{
#ifdef DEBUG
  opserr << "ContactMaterial3D::getTangent()\n";
#endif

  double C_nl;
  Matrix C_ss(2,2);
  Vector C_sl(2);

  double t_n = strain_vec(3);

  C_nl = 1.0;

  if (t_n < - tensileStrength) {
    C_ss.Zero();
    C_sl.Zero();
               
  } else if (inSlip) {
    // sliding coefficients
    Matrix r_dyadic_r(2,2);

    Vector R_nplus1 = g * r_nplus1;

    r_dyadic_r(0,0) = R_nplus1(0)*R_nplus1(0);
    r_dyadic_r(0,1) = R_nplus1(0)*R_nplus1(1);
    r_dyadic_r(1,1) = R_nplus1(1)*R_nplus1(1);
    r_dyadic_r(1,0) = r_dyadic_r(0,1);
           
    double scale = (1.0 - gamma/s_e_nplus1_norm);
    C_ss = stiffness * scale * (g - r_dyadic_r);
    C_sl = R_nplus1*frictionCoeff;

  } else {
    // sticking coefficients
    C_ss = stiffness * g;
    C_sl.Zero();
  }

#ifdef DEBUG
  if (DEBUG_LEVEL > 1) {
    opserr << "   strain_vec = " << strain_vec;
    opserr << "   is sliding? " << inSlip << endln;
  }
  if (DEBUG_LEVEL > 2) {
    opserr << "   C_nl = " << C_nl
           << "   C_ss = " << C_ss
           << "   C_sl = " << C_sl
           << endln;
    opserr << "   stiffness: " << stiffness
           << "   mu: " << frictionCoeff << endln;
 }
#endif

  //tangent matrix was zeroed initially
  tangent_matrix(0,3) = 1;
  tangent_matrix(1,1) = C_ss(0,0);
  tangent_matrix(1,2) = C_ss(0,1);
  tangent_matrix(2,1) = C_ss(1,0);
  tangent_matrix(2,2) = C_ss(1,1);

  tangent_matrix(1,3) = C_sl(0);
  tangent_matrix(2,3) = C_sl(1);
  tangent_matrix(3,0) = 1;
       
  return tangent_matrix;          
}


const Matrix & ContactMaterial3D::getInitialTangent ()
{
#ifdef DEBUG
  opserr << "ContactMaterial3D::getInitialTangent ()" << endln;
#endif
  return tangent_matrix;          //tangent is empty matrix
}


const Vector & ContactMaterial3D::getStress ()
{
#ifdef DEBUG
  opserr << "ContactMaterial3D::getStress ()" << endln;
#endif
  return stress_vec;
}


const Vector & ContactMaterial3D::getStrain ()
{
#ifdef DEBUG
  opserr << "ContactMaterial3D::getStrain ()" << endln;
#endif
  return strain_vec;
}


void ContactMaterial3D::setMetricTensor(Matrix &v)
{
#ifdef DEBUG
  opserr << "ContactMaterial3D::setMetricTensor(Matrix &v)" << endln;
#endif
  g = v;
 
  // dual basis metric tensor G = inv(g)
  double det = g(0,0)*g(1,1) - g(0,1)*g(1,0);
  G(0,0) =  g(1,1);
  G(1,0) = -g(1,0);
  G(0,1) = -g(0,1);
  G(1,1) =  g(0,0);
  G = G / det;
}


int ContactMaterial3D::commitState (void)
{
#ifdef DEBUG
  opserr << "ContactMaterial3D::commitState (void)" << endln;
#endif
  s_e_n = s_e_nplus1;
	
  return 0;
}
 

int ContactMaterial3D::revertToLastCommit (void)
{
#ifdef DEBUG
  opserr << "ContactMaterial3D::revertToLastCommit (void)" << endln;
#endif
  return 0;
}

int ContactMaterial3D::revertToStart(void)
{
#ifdef DEBUG
  opserr << "ContactMaterial3D::revertToStart(void)" << endln;
#endif

  this->zero();

  return 0;
}

NDMaterial * ContactMaterial3D::getCopy (void)
{
#ifdef DEBUG
  opserr << "ContactMaterial3D::getCopy (void)" << endln;
#endif
  ContactMaterial3D * copy = new ContactMaterial3D(*this);
  return copy;
}

NDMaterial * ContactMaterial3D::getCopy (const char *code)
{
#ifdef DEBUG
  opserr << "ContactMaterial3D::getCopy (const char *code)" << endln;
#endif
  if (strcmp(code,"ContactMaterial3D")==0) {
    ContactMaterial3D * copy = new ContactMaterial3D(*this);
    return copy;
  }

  return 0;
}

const char * ContactMaterial3D::getType (void) const
{
#ifdef DEBUG
  opserr << "ContactMaterial3D::getType (void) const" << endln;
#endif
  return "ThreeDimensional";

}

int ContactMaterial3D::getOrder (void) const
{
#ifdef DEBUG
  opserr << "ContactMaterial3D::getOrder (void) const" << endln;
#endif
  return 6;
}

int ContactMaterial3D::UpdateFrictionalState(void)
{
	if (mFrictFlag == 1 && mFlag == 1) {
		frictionCoeff = mMu;
		tensileStrength = mTen;
		cohesion = mCo;
		mFlag = 0;

		// ensure tensile strength is inbounds
		if (tensileStrength > cohesion / frictionCoeff ) {
			tensileStrength = cohesion / frictionCoeff;
		}
		
	} else if (mFrictFlag != 1) {
		frictionCoeff = 0.0;
		cohesion = 0.0;
		tensileStrength = 0.0;
		mFlag = 1;
	}

	return 0;
}

int ContactMaterial3D::setParameter(const char **argv, int argc, Parameter &param)
{
	return -1;
}

int ContactMaterial3D::sendSelf(int commitTag, Channel &theChannel)
{
#ifdef DEBUG
  opserr << "ContactMaterial3D::sendSelf(int commitTag, Channel &theChannel)" << endln;
#endif
  // we place all the data needed to define material and it's state
  // int a vector object
  static Vector data(29);
  data(0)  = this->getTag();
  data(1)  = mMu;
  data(2)  = mCo;
  data(3)  = mTen;
  data(4)  = mFrictFlag;
  data(5)  = mFlag;
  data(6)  = frictionCoeff;
  data(7)  = stiffness;
  data(8)  = cohesion;
  data(9)  = tensileStrength;
  data(10) = s_e_n(0);
  data(11) = s_e_n(1);
  data(12) = stress_vec(0);
  data(13) = stress_vec(1);
  data(14) = stress_vec(2);
  data(15) = stress_vec(3);
  data(16) = strain_vec(0);
  data(17) = strain_vec(1);
  data(18) = strain_vec(2);
  data(19) = strain_vec(3);
  data(20) = inSlip;
  data(21) = g(0,0);
  data(22) = g(0,1);
  data(23) = g(1,0);
  data(24) = g(1,1);
  data(25) = r_nplus1(0);
  data(26) = r_nplus1(1);
  data(27) = gamma;
  data(28) = s_e_nplus1_norm;

  // send the vector object to the channel
  if (theChannel.sendVector(this->getDbTag(), commitTag, data) < 0) {
    opserr << "ContactMaterial3D::sendSelf - failed to send vector to channel\n";
    return -1;
  }

  return 0;
 
}

int ContactMaterial3D::recvSelf(int commitTag, Channel &theChannel,FEM_ObjectBroker &theBroker)    
{
#ifdef DEBUG
  opserr << "ContactMaterial3D::recvSelf(...)" << endln;
#endif
  // recv the vector object from the channel which defines material param and state
  static Vector data(29);
  if (theChannel.recvVector(this->getDbTag(), commitTag, data) < 0) {
    opserr << "ContactMaterial3D::recvSelf - failed to recv vector from channel\n";
    return -1;
  }
  
  // set the material parameters and state variables
  this->setTag((int)data(0));
  mMu             = data(1);
  mCo             = data(2);
  mTen            = data(3);
  mFrictFlag      = data(4);
  mFlag           = data(5);
  frictionCoeff   = data(6);
  stiffness       = data(7);
  cohesion        = data(8);
  tensileStrength = data(9);
  s_e_n(0)        = data(10);
  s_e_n(1)        = data(11);
  stress_vec(0)   = data(12);
  stress_vec(1)   = data(13);
  stress_vec(2)   = data(14);
  stress_vec(3)   = data(15);
  strain_vec(0)   = data(16);
  strain_vec(1)   = data(17);
  strain_vec(2)   = data(18);
  strain_vec(3)   = data(19);
  inSlip          = (int)data(20);
  g(0,0)          = data(21);
  g(0,1)          = data(22);
  g(1,0)          = data(23);
  g(1,1)          = data(24);
  r_nplus1(0)     = data(25);
  r_nplus1(1)     = data(26);
  gamma           = data(27);
  s_e_nplus1_norm = data(28);

  s_e_nplus1 = s_e_n;

  return 0;

}


void ContactMaterial3D::Print(OPS_Stream &s, int flag )
{
#ifdef DEBUG
  opserr << "ContactMaterial3D::Print(OPS_Stream &s, int flag )" << endln;
#endif
  s << "ContactMaterial3D: " << endln;
}

int ContactMaterial3D::updateParameter(int responseID, Information &info)
{
  if (responseID == 1) {
    mFrictFlag = info.theDouble;
  }
  
  return 0;
}

double ContactMaterial3D::getcohesion(void)
{
    return cohesion;
}

void ContactMaterial3D::ScaleCohesion(const double len) 
{
    cohesion *= len;
}

double ContactMaterial3D::getTensileStrength(void)
{
    return tensileStrength;
}

void ContactMaterial3D::ScaleTensileStrength(const double len)
{
    tensileStrength *= len;
}
