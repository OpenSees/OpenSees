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
// $Source: /OpenSees/SRC/material/nD/ContactMaterial2D.cpp,v $
                                                                        
// Written: Kathryn Petek
// Created: February 2004
// Modified: Chris McGann
//           November 2010 -> changes for incorporation into main source code
// Modified: Chris McGann
//           Jan 2011 -> added update for frictional state

// Description: This file contains the implementation for the ContactMaterial2D class.
//				

#include <ContactMaterial2D.h>

#include <Information.h>
#include <MaterialResponse.h>
#include <Parameter.h>
#include <string.h>

#include <Channel.h>
#include <FEM_ObjectBroker.h>

#include <elementAPI.h>
#define OPS_Export
static int numContactMaterial2DMaterials = 0;
int ContactMaterial2D::mFrictFlag = 1;

void * OPS_ADD_RUNTIME_VPV(OPS_ContactMaterial2DMaterial)
{
  if (numContactMaterial2DMaterials == 0) {
    numContactMaterial2DMaterials++;
    opserr << "ContactMaterial2D nDmaterial - Written: K.Petek, P.Mackenzie-Helnwein, P.Arduino, U.Washington\n";
  }

  // Pointer to a nDmaterial that will be returned
  NDMaterial *theMaterial = 0;

  int numArgs = OPS_GetNumRemainingInputArgs();

  if (numArgs < 5) {
    opserr << "Want: nDMaterial ContactMaterial2D tag? mu? G? c? t?\n";
    return 0;	
  }
  
  int tag;
  double dData[4];

  int numData = 1;
  if (OPS_GetInt(&numData, &tag) != 0) {
    opserr << "WARNING invalid tag for  ContactMaterial2D material" << endln;
    return 0;
  }
  numData = 4;
  if (OPS_GetDouble(&numData, dData) != 0) {
    opserr << "WARNING invalid material data for nDMaterial ContactMaterial2D material  with tag: " << tag << endln;
    return 0;
  }

  theMaterial = new ContactMaterial2D(tag, dData[0], dData[1], dData[2], dData[3]);

  if (theMaterial == 0) {
    opserr << "WARNING ran out of memory for nDMaterial ContactMaterial2D material  with tag: " << tag << endln;
  }

  return theMaterial;
}

//full constructor
ContactMaterial2D::ContactMaterial2D (int tag, double mu, double G, double c, double t)
 : NDMaterial(tag,ND_TAG_ContactMaterial2D),
   strain_vec(3),
   stress_vec(3),
   tangent_matrix(3,3)
{
#ifdef DEBUG
        opserr << "ContactMaterial2D::ContactMaterial2D" << endln;
#endif
        frictionCoeff = mu;
		mMu = mu;
        stiffness = G;
        cohesion  = c;
		mCo = c;
		tensileStrength = t;
		mTen = t;

        this->zero();
}
   
//null constructor
ContactMaterial2D::ContactMaterial2D () 
 : NDMaterial(0, ND_TAG_ContactMaterial2D),
   strain_vec(3),
   stress_vec(3),
   tangent_matrix(3,3)
{
        frictionCoeff = 0.0;
        stiffness = 1.0;
        cohesion  = 0.0;
	tensileStrength = 0.0;
}

//destructor
ContactMaterial2D::~ContactMaterial2D ()
{
}

//zero internal variables
void ContactMaterial2D::zero( )
{
#ifdef DEBUG
        opserr << "ContactMaterial2D::zero()" << endln;
#endif

    s_e_n      = 0.0;       // elastic slip from previous increment
    s_e_nplus1 = 0.0;       // elastic slip after current increment
    r_nplus1   = 0.0;       // direction of plastic slip

    inSlip     = false;     
	mFlag = 1;

    stress_vec.Zero();
    strain_vec.Zero();
    tangent_matrix.Zero();

	

	// ensure that tensileStrength is within bounds
	if (frictionCoeff == 0.0) {
		tensileStrength = 0.0;
	}
	else if (tensileStrength > cohesion / frictionCoeff ) {
		tensileStrength = cohesion / frictionCoeff;
	}


}


int ContactMaterial2D::setTrialStrain (const Vector &strain_from_element)
{
#ifdef DEBUG
        opserr << "ContactMaterial2D::setTrialStrain()" << endln;
#endif

    double t_s;             // tangential contact force
    double t_n;             // normal contact force
    double f_nplus1_trial;  // trial slip condition
    double gamma;           // consistency parameter

    double gap;             // current gap
    double slip;            // incremental slip
   
    strain_vec = strain_from_element;

    gap    = strain_vec(0);
    slip   = strain_vec(1);
    t_n    = strain_vec(2);

	// update frictional status
	this->UpdateFrictionalState();

// trial state (elastic predictor step) -> assume sticking
        inSlip = false;

        s_e_nplus1 = (t_n > -tensileStrength) ? s_e_n + slip : 0.0;
        t_s        = stiffness * s_e_nplus1;

        // slip condition
        f_nplus1_trial = fabs(t_s) - frictionCoeff*t_n - cohesion;

        // if ((f_nplus1_trial > 0.0) ) {
        if ((f_nplus1_trial > 0.0) && (t_n > -tensileStrength) && (fabs(s_e_nplus1) > 1.0e-12))       {


// plastic corrector step -> sliding
            inSlip = true;

            gamma = f_nplus1_trial / stiffness ;

            r_nplus1 = (t_s < 0) ? -1 : 1;

            // s_p_nplus1 = s_p_n + gamma * r_nplus1
            // s_e_nplus1 = s_nplus1 - s_p_nplus1
            //        = (s_nplus1 - s_p_n) - gamma * r_nplus1
            //        = (s_n + slip - s_p_n) - gamma * r_nplus1
            //        = (s_e_n + slip) - gamma * r_nplus1
            //        = s_e_nplus1_trial - gamma * r_nplus1
            s_e_nplus1 = s_e_nplus1 - gamma * r_nplus1;

            t_s = stiffness * s_e_nplus1;
    }

#ifdef DEBUG
    if (DEBUG_LEVEL > 1) {
        if (inSlip) {
            opserr << "   ** SLIDING" << endln; }
        else {
            opserr << "   ** STICKING" << endln;}
    }
#endif

    //update stress and strain values
    stress_vec(0) = t_n;
    stress_vec(1) = t_s;
    stress_vec(2) = gap;
    
    return 0;

}


//unused trial strain functoins
int ContactMaterial2D::setTrialStrain (const Vector &v, const Vector &r)
{
  return this->setTrialStrain (v);
}


const Matrix & ContactMaterial2D::getTangent ()
{
#ifdef DEBUG
    opserr << "ContactMaterial2D::getTangent()\n";
#endif

    double C_nl;
    double C_ss;
    double C_sl;

    double t_n = strain_vec(2);

    C_nl = 1.0;

//    C_ss = (inSlip || t_n < - tensileStrength) ? 0.0 : stiffness;
//    C_sl = (inSlip) ? r_nplus1*frictionCoeff : 0.0;


	if (t_n < - tensileStrength) {
		C_ss = 0.0;
		C_sl = 0.0;
		
	} else if (inSlip) {
    // sliding coefficients
		C_ss = 0.0;
		C_sl = r_nplus1*frictionCoeff; 

	} else {
	// sticking coefficients
		C_ss = stiffness;
		C_sl = 0.0;
	} 
	

#ifdef DEBUG
    if (DEBUG_LEVEL > 1) {
        opserr << "   is sliding? " << inSlip << endln;
        opserr << "   C_nl = " << C_nl
               << "   C_ss = " << C_ss
               << "   C_sl = " << C_sl
               << endln;
        opserr << "   stiffness: " << stiffness 
               << "   mu: " << frictionCoeff << endln;
    }
#endif

//tangent matrix was zeroed initially
    tangent_matrix(0,2) = 1;
    tangent_matrix(1,1) = C_ss;
    tangent_matrix(1,2) = C_sl;
    tangent_matrix(2,0) = 1;
    
    return tangent_matrix;      
}


const Matrix & ContactMaterial2D::getInitialTangent ()
{
#ifdef DEBUG
        opserr << "ContactMaterial2D::getInitialTangent()" << endln;
#endif

    return tangent_matrix;      //tangent is empty matrix
}


const Vector & ContactMaterial2D::getStress()
{
#ifdef DEBUG
        opserr << "ContactMaterial2D::getStress()" << endln;
#endif

    return stress_vec;
}


const Vector & ContactMaterial2D::getStrain ()
{
#ifdef DEBUG
        opserr << "ContactMaterial2D::setStrain()" << endln;
#endif

    return strain_vec;
}


int ContactMaterial2D::commitState (void)
{
#ifdef DEBUG
        opserr << "ContactMaterial2D::commitState" << endln;
#endif

        s_e_n = s_e_nplus1;
  
        return 0;
}
 

int ContactMaterial2D::revertToLastCommit (void)
{
#ifdef DEBUG
        opserr << "ContactMaterial2D::revertToLastCommit()" << endln;
#endif

        return 0;
}

int ContactMaterial2D::revertToStart(void)
{
#ifdef DEBUG
        opserr << "ContactMaterial2D::revertToStart()" << endln;
#endif
	
    this->zero();

    return 0;
}


NDMaterial * ContactMaterial2D::getCopy (void)
{
#ifdef DEBUG
        opserr << "ContactMaterial2D::getCopy()" << endln;
#endif

  ContactMaterial2D * copy = new ContactMaterial2D(*this);
  return copy;
}


NDMaterial * ContactMaterial2D::getCopy (const char *code)
{
#ifdef DEBUG
        opserr << "ContactMaterial2D::getCopy()" << endln;
#endif

  if (strcmp(code,"ContactMaterial2D")==0) {
    ContactMaterial2D * copy = new ContactMaterial2D(*this);
    return copy;
  }

  return 0;
}



const char * ContactMaterial2D::getType (void) const
{
#ifdef DEBUG
        opserr << "ContactMaterial2D::getType()" << endln;
#endif
    return "Two_Dimensional";

}


int ContactMaterial2D::getOrder (void) const
{
#ifdef DEBUG
        opserr << "ContactMaterial2D::getOrder()" << endln;
#endif
    return 3;
}

int ContactMaterial2D::UpdateFrictionalState(void)
{
	if (mFrictFlag == 1 && mFlag == 1) {
		frictionCoeff = mMu;
		cohesion = mCo;
		tensileStrength = mTen;
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

int ContactMaterial2D::setParameter(const char **argv, int argc, Parameter &param)
{
	return -1;
}

int ContactMaterial2D::sendSelf(int commitTag, Channel &theChannel)
{
  // we place all the data needed to define material and it's state
  // int a vector object
  static Vector data(6);
  int cnt = 0;
  data(cnt++) = this->getTag();
  data(cnt++) = frictionCoeff;
  data(cnt++) = stiffness;
  data(cnt++) = cohesion;
  data(cnt++) = tensileStrength;
  data(cnt++) = s_e_n;

  // send the vector object to the channel
  if (theChannel.sendVector(this->getDbTag(), commitTag, data) < 0) {
    opserr << "ContactMaterial2D::sendSelf - failed to send vector to channel\n";
    return -1;
  }


  return 0;
 
}


int ContactMaterial2D::recvSelf(int commitTag, Channel &theChannel, 
                     FEM_ObjectBroker &theBroker)    
{
  // recv the vector object from the channel which defines material param and state
  static Vector data(5);
  if (theChannel.recvVector(this->getDbTag(), commitTag, data) < 0) {
    opserr << "ContactMaterial2D::recvSelf - failed to recv vector from channel\n";
    return -1;
  }

  // set the material parameters and state variables
  int cnt = 0;
  this->setTag(data(cnt++));
  frictionCoeff = data(cnt++);
  stiffness = data(cnt++);
  cohesion = data(cnt++);
  tensileStrength = data(cnt++);
  s_e_n = data(cnt++);

  s_e_nplus1 = s_e_n;

  return 0;

}


void ContactMaterial2D::Print(OPS_Stream &s, int flag )
{
  s << "ContactMaterial2D" << endln;
}


int ContactMaterial2D::updateParameter(int responseID, Information &info)
{
    if (responseID==20) frictionCoeff=info.theDouble;
	if (responseID==21) stiffness=info.theDouble;

	if (responseID == 1) {
		mFrictFlag = info.theDouble;
	}

  return 0;
}

double ContactMaterial2D::getcohesion(void)
{
    return cohesion;
}

void ContactMaterial2D::ScaleCohesion(const double len) 
{
    cohesion *= len;
    //return 0;
}

double ContactMaterial2D::getTensileStrength(void)
{
    return tensileStrength;
}

void ContactMaterial2D::ScaleTensileStrength(const double len)
{
    tensileStrength *= len;
    //return 0;
}
