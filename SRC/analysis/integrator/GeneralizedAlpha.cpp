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

// $Revision$
// $Date$
// $URL$

// Written: fmk
// Created: 1/10
// Revision: A
//
// Description: This file contains the implementation of the GeneralizedAlpha class.
// J.Chung, G.M.Hulbert "A Time Integration Algorithm for Structural Dynamics With
// Improved Numerical Dissipation: The Generalized-alpha Method" ASME Journal of Applied
// Mechanics, Vol 60, 371-375, 1993.

#include <GeneralizedAlpha.h>
#include <FE_Element.h>
#include <LinearSOE.h>
#include <AnalysisModel.h>
#include <Vector.h>
#include <DOF_Group.h>
#include <DOF_GrpIter.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>

#include <elementAPI.h>
#define OPS_Export 

void *OPS_GeneralizedAlpha(void)
{
  // Pointer to an integrator that will be returned
  TransientIntegrator *theIntegrator = 0;

  int argc = OPS_GetNumRemainingInputArgs();
  if (argc != 2 && argc != 4) {
    opserr << "WARNING - incorrect number of args want GeneralizedAlpha $alphaM $alphaF <$gamma $beta>\n";
    return 0;
  }

  double dData[4];
  if (OPS_GetDouble(&argc, dData) != 0) {
    opserr << "WARNING - invalid args want GeneralizedAlpha $alphaM $alphaF <$gamma $beta>\n";
    return 0;
  }
  
  if (argc == 2)
    theIntegrator = new GeneralizedAlpha(dData[0], dData[1]);
  else
    theIntegrator = new GeneralizedAlpha(dData[0], dData[1], dData[2], dData[3]);

  if (theIntegrator == 0)
    opserr << "WARNING - out of memory creating GeneralizedAlpha integrator\n";

  return theIntegrator;
}


GeneralizedAlpha::GeneralizedAlpha()
  : TransientIntegrator(INTEGRATOR_TAGS_GeneralizedAlpha),
    alphaM(0.0), alphaF(0.0), gamma(0.0),
    deltaT(0.0), 
    c1(0.0), c2(0.0), c3(0.0), 
    Ut(0), Utdot(0), Utdotdot(0),
    U(0), Udot(0), Udotdot(0),
    Ualpha(0), Ualphadot(0), Ualphadotdot(0)
{
    
}


GeneralizedAlpha::GeneralizedAlpha(double _alphaM, double _alphaF)
  : TransientIntegrator(INTEGRATOR_TAGS_GeneralizedAlpha),
    alphaM(_alphaM), alphaF(_alphaF),
    beta((1+_alphaM-_alphaF)*(1+_alphaM-_alphaF)*0.25), gamma(0.5+_alphaM-_alphaF),
    deltaT(0.0), 
    c1(0.0), c2(0.0), c3(0.0), 
    Ut(0), Utdot(0), Utdotdot(0), 
    U(0), Udot(0), Udotdot(0),
    Ualpha(0), Ualphadot(0), Ualphadotdot(0)
{
    
}


GeneralizedAlpha::GeneralizedAlpha(double _alphaM, double _alphaF, double _beta, double _gamma)
    : TransientIntegrator(INTEGRATOR_TAGS_GeneralizedAlpha),
      alphaM(_alphaM), alphaF(_alphaF),
      beta(_beta), gamma(_gamma),
      deltaT(0.0), 
      c1(0.0), c2(0.0), c3(0.0), 
      Ut(0), Utdot(0), Utdotdot(0), 
      U(0), Udot(0), Udotdot(0),
      Ualpha(0), Ualphadot(0), Ualphadotdot(0)
{
    
}


GeneralizedAlpha::~GeneralizedAlpha()
{
    // clean up the memory created
    if (Ut != 0)
        delete Ut;
    if (Utdot != 0)
        delete Utdot;
    if (Utdotdot != 0)
        delete Utdotdot;
    if (U != 0)
        delete U;
    if (Udot != 0)
        delete Udot;
    if (Udotdot != 0)
        delete Udotdot;
    if (Ualpha != 0)
        delete Ualpha;
    if (Ualphadot != 0)
        delete Ualphadot;
    if (Ualphadotdot != 0)
        delete Ualphadotdot;
}


int GeneralizedAlpha::newStep(double _deltaT)
{
    deltaT = _deltaT;
    if (beta == 0 || gamma == 0 )  {
        opserr << "GeneralizedAlpha::newStep() - error in variable\n";
        opserr << "gamma = " << gamma << " beta = " << beta << endln;
        return -1;
    }
    
    if (deltaT <= 0.0)  {
        opserr << "GeneralizedAlpha::newStep() - error in variable\n";
        opserr << "dT = " << deltaT << endln;
        return -2;	
    }

    // get a pointer to the AnalysisModel
    AnalysisModel *theModel = this->getAnalysisModel();
    
    // set the constants
    c1 = 1.0;
    c2 = gamma/(beta*deltaT);
    c3 = 1.0/(beta*deltaT*deltaT);
       
    if (U == 0)  {
        opserr << "GeneralizedAlpha::newStep() - domainChange() failed or hasn't been called\n";
        return -3;
    }
    
    // set response at t to be that at t+deltaT of previous step
    (*Ut) = *U;
    (*Utdot) = *Udot;
    (*Utdotdot) = *Udotdot;

    // determine new velocities and accelerations at t+deltaT
    double a1 = (1.0 - gamma/beta);
    double a2 = deltaT*(1.0 - 0.5*gamma/beta);
    Udot->addVector(a1, *Utdotdot, a2);
    
    double a3 = -1.0/(beta*deltaT);
    double a4 = 1.0 - 0.5/beta;
    Udotdot->addVector(a4, *Utdot, a3);

    // determine the velocities at t+alphaF*deltaT
    (*Ualphadot) = *Utdot;
    Ualphadot->addVector((1.0-alphaF), *Udot, alphaF);

    // determine the velocities at t+alphaM*deltaT
    (*Ualphadotdot) = *Utdotdot;
    Ualphadotdot->addVector((1.0-alphaM), *Udotdot, alphaM);
    
    // set the trial response quantities
    theModel->setVel(*Ualphadot);
    theModel->setAccel(*Ualphadotdot);
        
    // increment the time to t+alpha*deltaT and apply the load
    double time = theModel->getCurrentDomainTime();
    time += alphaF*deltaT;

    if (theModel->updateDomain(time, deltaT) < 0)  {
        opserr << "GeneralizedAlpha::newStep() - failed to update the domain\n";
        return -4;
    }

    return 0;
}


int GeneralizedAlpha::revertToLastStep()
{
    // set response at t+deltaT to be that at t .. for next step
    if (U != 0)  {
        (*U) = *Ut;
        (*Udot) = *Utdot;
        (*Udotdot) = *Utdotdot;
    }

    return 0;
}


int GeneralizedAlpha::formEleTangent(FE_Element *theEle)
{
    theEle->zeroTangent();
    if (statusFlag == CURRENT_TANGENT)  {
        theEle->addKtToTang(alphaF*c1);
        theEle->addCtoTang(alphaF*c2);
        theEle->addMtoTang(alphaM*c3);
    } else if (statusFlag == INITIAL_TANGENT)  {
        theEle->addKiToTang(alphaF*c1);
        theEle->addCtoTang(alphaF*c2);
        theEle->addMtoTang(alphaM*c3);
    } else if (statusFlag == HALL_TANGENT)  {
        theEle->addKtToTang(c1*cFactor);
        theEle->addKiToTang(c1*iFactor);
        theEle->addCtoTang(c2);
        theEle->addMtoTang(c3);
    }    
    
    return 0;
}   
 

int GeneralizedAlpha::formNodTangent(DOF_Group *theDof)
{
    theDof->zeroTangent();

    theDof->addCtoTang(alphaF*c2);
    theDof->addMtoTang(alphaM*c3);
    
    return 0;
}


int GeneralizedAlpha::domainChanged()
{
    AnalysisModel *myModel = this->getAnalysisModel();
    LinearSOE *theLinSOE = this->getLinearSOE();
    const Vector &x = theLinSOE->getX();
    int size = x.Size();
    
    // create the new Vector objects
    if (Ut == 0 || Ut->Size() != size)  {
      
        // delete the old
      if (Ut != 0) {
	delete Ut;
	delete Utdot;
	delete Utdotdot;
	delete U;
	delete Udot;
	delete Udotdot;
	delete Ualpha;
	delete Ualphadot;
	delete Ualphadotdot;
      }        
      // create the new
      Ut = new Vector(size);
      Utdot = new Vector(size);
      Utdotdot = new Vector(size);
      U = new Vector(size);
      Udot = new Vector(size);
      Udotdot = new Vector(size);
      Ualpha = new Vector(size);
      Ualphadot = new Vector(size);
      Ualphadotdot = new Vector(size);
      
      // check we obtained the new
      if (Ut == 0 || Ut->Size() != size ||
	  Utdot == 0 || Utdot->Size() != size ||
	  Utdotdot == 0 || Utdotdot->Size() != size ||
	  U == 0 || U->Size() != size ||
	  Udot == 0 || Udot->Size() != size ||
	  Udotdot == 0 || Udotdot->Size() != size ||
	  Ualpha == 0 || Ualpha->Size() != size ||
	  Ualphadot == 0 || Ualphadot->Size() != size ||
	  Ualphadotdot == 0 || Ualphadotdot->Size() != size)  {
	
	opserr << "GeneralizedAlpha::domainChanged - ran out of memory\n";
        
	// delete the old
	if (Ut != 0)
	  delete Ut;
	if (Utdot != 0)
	  delete Utdot;
	if (Utdotdot != 0)
	  delete Utdotdot;
	if (U != 0)
	  delete U;
	if (Udot != 0)
	  delete Udot;
	if (Udotdot != 0)
	  delete Udotdot;
	if (Ualpha != 0)
	  delete Ualpha;
	if (Ualphadot != 0)
	  delete Ualphadot;
	if (Ualphadotdot != 0)
	  delete Ualphadotdot;
	
	Ut = 0; Utdot = 0; Utdotdot = 0;
	U = 0; Udot = 0; Udotdot = 0;
	Ualpha = 0; Ualphadot = 0; Ualphadotdot = 0;
	
	return -1;
      }
    }        
    
    // now go through and populate U, Udot and Udotdot by iterating through
    // the DOF_Groups and getting the last committed velocity and accel
    DOF_GrpIter &theDOFs = myModel->getDOFs();
    DOF_Group *dofPtr;
    while ((dofPtr = theDOFs()) != 0)  {
      const ID &id = dofPtr->getID();
      int idSize = id.Size();
      
      int i;
      const Vector &disp = dofPtr->getCommittedDisp();	
      for (i=0; i < idSize; i++)  {
            int loc = id(i);
            if (loc >= 0)  {
                (*U)(loc) = disp(i);		
            }
        }
        
        const Vector &vel = dofPtr->getCommittedVel();
        for (i=0; i < idSize; i++)  {
            int loc = id(i);
            if (loc >= 0)  {
                (*Udot)(loc) = vel(i);
            }
        }
        
        const Vector &accel = dofPtr->getCommittedAccel();	
        for (i=0; i < idSize; i++)  {
            int loc = id(i);
            if (loc >= 0)  {
                (*Udotdot)(loc) = accel(i);
            }
        }        
    }    
    
    return 0;
}


int GeneralizedAlpha::update(const Vector &deltaU)
{
    AnalysisModel *theModel = this->getAnalysisModel();
    if (theModel == 0)  {
        opserr << "WARNING GeneralizedAlpha::update() - no AnalysisModel set\n";
        return -1;
    }
    
    // check domainChanged() has been called, i.e. Ut will not be zero
    if (Ut == 0)  {
        opserr << "WARNING GeneralizedAlpha::update() - domainChange() failed or not called\n";
        return -2;
    }
    
    // check deltaU is of correct size
    if (deltaU.Size() != U->Size())  {
        opserr << "WARNING GeneralizedAlpha::update() - Vectors of incompatible size ";
        opserr << " expecting " << U->Size() << " obtained " << deltaU.Size() << endln;
        return -3;
    }
    
    //  determine the response at t+deltaT
    (*U) += deltaU;
    Udot->addVector(1.0, deltaU, c2);
    Udotdot->addVector(1.0, deltaU, c3);

    // determine displacement and velocity at t+alphaF*deltaT
    (*Ualpha) = *Ut;
    Ualpha->addVector((1.0-alphaF), *U, alphaF);

    (*Ualphadot) = *Utdot;
    Ualphadot->addVector((1.0-alphaF), *Udot, alphaF);

    // determine the velocities at t+alphaM*deltaT
    (*Ualphadotdot) = *Utdotdot;
    Ualphadotdot->addVector((1.0-alphaM), *Udotdot, alphaM);

    
    // update the response at the DOFs
    theModel->setResponse(*Ualpha,*Ualphadot,*Udotdot);        
    if (theModel->updateDomain() < 0)  {
        opserr << "GeneralizedAlpha::update() - failed to update the domain\n";
        return -4;
    }
    
    return 0;
}


int GeneralizedAlpha::commit(void)
{
    AnalysisModel *theModel = this->getAnalysisModel();
    if (theModel == 0)  {
        opserr << "WARNING GeneralizedAlpha::commit() - no AnalysisModel set\n";
        return -1;
    }	  
    
    // update the response at the DOFs
    theModel->setResponse(*U,*Udot,*Udotdot);
    if (theModel->updateDomain() < 0)  {
        opserr << "GeneralizedAlpha::commit() - failed to update the domain\n";
        return -4;
    }
    
    // set the time to be t+deltaT
    double time = theModel->getCurrentDomainTime();
    time += (1.0-alphaF)*deltaT;
    theModel->setCurrentDomainTime(time);

    return theModel->commitDomain();
}

const Vector &
GeneralizedAlpha::getVel()
{
  return *Udot;
}

int GeneralizedAlpha::sendSelf(int cTag, Channel &theChannel)
{
    Vector data(4);
    data(0) = alphaF;
    data(1) = alphaM;
    data(2) = beta;
    data(3) = gamma;

    
    if (theChannel.sendVector(this->getDbTag(), cTag, data) < 0)  {
        opserr << "WARNING GeneralizedAlpha::sendSelf() - could not send data\n";
        return -1;
    }

    return 0;
}


int GeneralizedAlpha::recvSelf(int cTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
    Vector data(4);

    if (theChannel.recvVector(this->getDbTag(), cTag, data) < 0)  {
        opserr << "WARNING GeneralizedAlpha::recvSelf() - could not receive data\n";
        return -1;
    }
    
    alphaF  = data(0);
    alphaM  = data(1);
    beta   = data(2);
    gamma  = data(3);
    
    return 0;
}


void GeneralizedAlpha::Print(OPS_Stream &s, int flag)
{
  AnalysisModel *theModel = this->getAnalysisModel();
  if (theModel != 0)  {
    double currentTime = theModel->getCurrentDomainTime();
    s << "\t GeneralizedAlpha - currentTime: " << currentTime << endln;
    s << "  alphaF: " << alphaF << "  alphaM: " << alphaM << "  beta: " << beta  << "  gamma: " << gamma << endln;
    s << "  c1: " << c1 << "  c2: " << c2 << "  c3: " << c3 << endln;
  } else 
    s << "\t GeneralizedAlpha - no associated AnalysisModel\n";
}
