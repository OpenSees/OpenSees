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
                                                                        
// $Revision: 1.8 $
// $Date: 2007-04-02 23:42:26 $
// $Source: /usr/local/cvs/OpenSees/SRC/analysis/integrator/ArcLength.cpp,v $
                                                                        
                                                                        
// File: ~/analysis/integrator/ArcLength.C
// 
// Written: fmk 
// Created: 07/98
// Revision: A
//
// Description: This file contains the class definition for ArcLength.
// ArcLength is an algorithmic class for perfroming a static analysis
// using the arc length scheme, that is within a load step the follwing
// constraint is enforced: dU^TdU + alpha^2*dLambda^2 = arcLength^2
// where dU is change in nodal displacements for step, dLambda is
// change in applied load and arcLength is a control parameter.
//
// What: "@(#) ArcLength.C, revA"


#include <ArcLength.h>
#include <AnalysisModel.h>
#include <LinearSOE.h>
#include <Vector.h>
#include <Channel.h>
#include <math.h>
#include <stdlib.h>
#include <elementAPI.h>
#include<Domain.h>
#include<ID.h>
#include<FE_Element.h>
#include<FE_EleIter.h>
#include<LoadPattern.h>
#include<LoadPatternIter.h>
#include<Parameter.h>
#include<ParameterIter.h>
#include<Node.h>
#include<DOF_Group.h>
#include<DOF_GrpIter.h>
#include<TaggedObjectStorage.h>
#include<EquiSolnAlgo.h>
void* OPS_ArcLength()
{
    double arcLength;
    double alpha;
    if (OPS_GetNumRemainingInputArgs() < 2) {
	opserr << "WARNING integrator ArcLength arcLength alpha \n";
	return 0;
    }

    int numdata = 1;
    if (OPS_GetDoubleInput(&numdata, &arcLength) < 0) {
	opserr << "WARNING integrator ArcLength failed to read arc lenght\n";
	return 0;
    }
    if (OPS_GetDoubleInput(&numdata, &alpha) < 0) {
	opserr << "WARNING integrator ArcLength failed to read alpha\n";
	return 0;
    }
    return new ArcLength(arcLength,alpha); 
}

ArcLength::ArcLength(double arcLength, double alpha)
:StaticIntegrator(INTEGRATOR_TAGS_ArcLength),
 arcLength2(arcLength*arcLength), alpha2(alpha*alpha),
 deltaUhat(0), deltaUbar(0), deltaU(0), deltaUstep(0),deltaUstep2(0),dDeltaUstepdh(0), 
 phat(0), deltaLambdaStep(0.0),dDeltaLambdaStepdh(0.0), currentLambda(0.0), dLAMBDA(0.0),dLAMBDA2(0.0),dlambda1dh(0.0),dLAMBDAdh(0),Residual(0),sensU(0),sensitivityFlag(0),
 signLastDeltaLambdaStep(1), dUhatdh(0),dphatdh(0),dUIJdh(0),dlambdaJdh(0.0),gradNumber(0), a(0.0),b(0.0),c(0.0),b24ac(0.0)
{

}

ArcLength::~ArcLength()
{
    // delete any vector object created
    if (deltaUhat != 0)
	delete deltaUhat;
    if (deltaU != 0)
	delete deltaU;
    if (deltaUstep != 0)
	delete deltaUstep;
    if(deltaUstep2 !=0)
       delete deltaUstep2;
    if (deltaUbar != 0)
	delete deltaUbar;
    if (phat != 0)
	delete phat;
    if(dUhatdh !=0)
       delete dUhatdh;
    if(dphatdh !=0)
       delete dphatdh;
    if(dLAMBDAdh !=0)
       delete dLAMBDAdh;
    if(dUIJdh !=0)
      delete dUIJdh;
    if(dDeltaUstepdh !=0)
       delete dDeltaUstepdh;
    if(Residual !=0)
       delete Residual;
    if(sensU !=0)
       delete sensU;
}

int
ArcLength::newStep(void)
{
  // opserr<<"sewStep: start"<<endln;
//opserr<<"alpha2 is "<<alpha2<<endln;
    // get pointers to AnalysisModel and LinearSOE
    AnalysisModel *theModel = this->getAnalysisModel();
    LinearSOE *theLinSOE = this->getLinearSOE();    
    if (theModel == 0 || theLinSOE == 0) {
	opserr << "WARNING ArcLength::newStep() ";
	opserr << "No AnalysisModel or LinearSOE has been set\n";
	return -1;
    }

    // get the current load factor
    currentLambda = theModel->getCurrentDomainTime();

    if (deltaLambdaStep < 0)
	signLastDeltaLambdaStep = -1;
    else
	signLastDeltaLambdaStep = +1;

    // determine dUhat
    this->formTangent();
    theLinSOE->setB(*phat);
   // opserr<<" ArcLength: Phat is "<<*phat<<endln;
    if (theLinSOE->solve() < 0) {
      opserr << "ArcLength::newStep(void) - failed in solver\n";
      return -1;
    }

    (*deltaUhat) = theLinSOE->getX();
    Vector &dUhat = *deltaUhat;
    
    // determine delta lambda(1) == dlambda
    double dLambda = sqrt(arcLength2/((dUhat^dUhat)+alpha2));
    dLambda *= signLastDeltaLambdaStep; // base sign of load change
   //    opserr<<"newStep:   the sign is "<<signLastDeltaLambdaStep<<endln;     // on what was happening last step
    deltaLambdaStep = dLambda;
    dLAMBDA=dLambda;
  //  opserr<<"newStep:   dLAMBDA= "<<dLAMBDA<<endln;
    currentLambda += dLambda;

    // determine delta U(1) == dU
    (*deltaU) = dUhat;
    (*deltaU) *= dLambda;
    (*deltaUstep) = (*deltaU);
    (*deltaUstep2)=(*deltaU);

  // update model with delta lambda and delta U
    theModel->incrDisp(*deltaU);    
    //////

   if(this->activateSensitivity()==true) { 
          Domain *theDomain=theModel->getDomainPtr();
      ParameterIter &paramIter = theDomain->getParameters();
      Parameter *theParam;

      // De-activate all parameters

      // Now, compute sensitivity wrt each parameter
      //int numGrads = theDomain->getNumParameters();

      while ((theParam = paramIter()) != 0)
	 theParam->activate(false);

      paramIter = theDomain->getParameters();
      while ((theParam = paramIter()) != 0) {
	 // Activate this parameter
	 theParam->activate(true);
	 // Get the grad index for this parameter
	 gradNumber = theParam->getGradIndex();

	 this->formTangDispSensitivity(gradNumber);

	 this->formdLambdaDh(gradNumber);
	// sensU->addVector(1.0,*dUhatdh,dLambda);
	// sensU->addVector(1.0,*deltaUhat,dlambda1dh);

	 dDeltaUstepdh->addVector(0.0,*dUhatdh,dLambda);
	 dDeltaUstepdh->addVector(1.0,*deltaUhat,dlambda1dh);
	 dDeltaLambdaStepdh =dlambda1dh; 
	 theParam->activate(false);
      } 
   }
   ///////////////Abbas/////////////////////////////


         /////
    theModel->applyLoadDomain(currentLambda);    
    theModel->updateDomain();
//opserr<<" newStep:   end"<<endln;

    return 0;
}

int
ArcLength::update(const Vector &dU)
{
  /// opserr<<" update: start"<<endln;

    AnalysisModel *theModel = this->getAnalysisModel();
    LinearSOE *theLinSOE = this->getLinearSOE();    
    if (theModel == 0 || theLinSOE == 0) {
	opserr << "WARNING ArcLength::update() ";
	opserr << "No AnalysisModel or LinearSOE has been set\n";
	return -1;
    }

    (*deltaUbar) = dU; // have to do this as the SOE is gonna change
//opserr<<"deltaUbar= "<<*deltaUbar<<endln;

//opserr<<"Update:   phat is = "<<*phat<<"////////////////////////////"<<endln;
    // determine dUhat    
    theLinSOE->setB(*phat);
    theLinSOE->solve();

    (*deltaUhat) = theLinSOE->getX();    

    // determine the coeeficients of our quadratic equation
           a = alpha2 + ((*deltaUhat)^(*deltaUhat));
           b = alpha2*deltaLambdaStep 
      + ((*deltaUhat)^(*deltaUbar))
      + ((*deltaUstep)^(*deltaUhat));
      b *= 2.0;
      c = 2*((*deltaUstep)^(*deltaUbar)) + ((*deltaUbar)^(*deltaUbar));
    // check for a solution to quadratic
         b24ac = b*b - 4.0*a*c;
    if (b24ac < 0) {
      opserr << "ArcLength::update() - imaginary roots due to multiple instability";
      opserr << " directions - initial load increment was too large\n";
      opserr << "a: " << a << " b: " << b << " c: " << c << " b24ac: " << b24ac << endln;
      return -1;
    }			       
    double a2 = 2.0*a;
    if (a2 == 0.0) {
      opserr << "ArcLength::update() - zero denominator";
      opserr << " alpha was set to 0.0 and zero reference load\n";
      return -2;
    }			       

    // determine the roots of the quadratic
    double sqrtb24ac = sqrt(b24ac);
    double dlambda1 = (-b + sqrtb24ac)/a2;
    double dlambda2 = (-b - sqrtb24ac)/a2;
//opserr<<"squareRoot of b24ac= "<<sqrtb24ac<<endln;
    double val = (*deltaUhat)^(*deltaUstep);
    double theta1 = ((*deltaUstep)^(*deltaUstep)) + ((*deltaUbar)^(*deltaUstep));
    //    double theta2 = theta1 + dlambda2*val;
    theta1 += dlambda1*val;

    // choose dLambda based on angle between incremental displacement before
    // and after this step -- want positive
    double dLambda;
    if (theta1 > 0)
      dLambda = dlambda1;
    else
      dLambda = dlambda2;

dLAMBDA2=dLambda;
    // determine delta U(i)
    (*deltaU) = (*deltaUbar);    
    deltaU->addVector(1.0, *deltaUhat,dLambda);
    
    // update dU and dlambda
    (*deltaUstep) += *deltaU;
    deltaLambdaStep += dLambda;
    currentLambda += dLambda;
///
//This is just a check
//double Result=a*(dLambda*dLambda)+(b*dLambda)+c;
    // update the model
    theModel->incrDisp(*deltaU);    
    theModel->applyLoadDomain(currentLambda);    


    theModel->updateDomain();
    
    // set the X soln in linearSOE to be deltaU for convergence Test
    theLinSOE->setX(*deltaU);
//opserr<<" update:  end"<<endln;
    return 0;
}



int 
ArcLength::domainChanged(void)
{
 //  opserr<<"domainChanged: start"<<endln;

    // we first create the Vectors needed
    AnalysisModel *theModel = this->getAnalysisModel();
    LinearSOE *theLinSOE = this->getLinearSOE();    
    if (theModel == 0 || theLinSOE == 0) {
	opserr << "WARNING ArcLength::update() ";
	opserr << "No AnalysisModel or LinearSOE has been set\n";
	return -1;
    }    
    int size = theModel->getNumEqn(); // ask model in case N+1 space

    if (deltaUhat == 0 || deltaUhat->Size() != size) { // create new Vector
	if (deltaUhat != 0)
	    delete deltaUhat;   // delete the old
	deltaUhat = new Vector(size);
	if (deltaUhat == 0 || deltaUhat->Size() != size) { // check got it
	    opserr << "FATAL ArcLength::domainChanged() - ran out of memory for";
	    opserr << " deltaUhat Vector of size " << size << endln;
	    exit(-1);
	}
    }

    if (deltaUbar == 0 || deltaUbar->Size() != size) { // create new Vector
	if (deltaUbar != 0)
	    delete deltaUbar;   // delete the old
	deltaUbar = new Vector(size);
	if (deltaUbar == 0 || deltaUbar->Size() != size) { // check got it
	    opserr << "FATAL ArcLength::domainChanged() - ran out of memory for";
	    opserr << " deltaUbar Vector of size " << size << endln;
	    exit(-1);
	}
    }

    
    if (deltaU == 0 || deltaU->Size() != size) { // create new Vector
	if (deltaU != 0)
	    delete deltaU;   // delete the old
	deltaU = new Vector(size);
	if (deltaU == 0 || deltaU->Size() != size) { // check got it
	    opserr << "FATAL ArcLength::domainChanged() - ran out of memory for";
	    opserr << " deltaU Vector of size " << size << endln;
	    exit(-1);
	}
    }

    if (deltaUstep == 0 || deltaUstep->Size() != size) { 
	if (deltaUstep != 0)
	    delete deltaUstep;  
	deltaUstep = new Vector(size);
	if (deltaUstep == 0 || deltaUstep->Size() != size) { 
	    opserr << "FATAL ArcLength::domainChanged() - ran out of memory for";
	    opserr << " deltaUstep Vector of size " << size << endln;
	    exit(-1);
	}
    }

if (deltaUstep2 == 0 || deltaUstep2->Size() != size) { 
	if (deltaUstep2 != 0)
	    delete deltaUstep2;  
	deltaUstep2 = new Vector(size);
	if (deltaUstep2 == 0 || deltaUstep2->Size() != size) { 
	    opserr << "FATAL ArcLength::domainChanged() - ran out of memory for";
	    opserr << " deltaUstep2 Vector of size " << size << endln;
	    exit(-1);
	}
    }


    if (dDeltaUstepdh == 0 || dDeltaUstepdh->Size() != size) { 
	if (dDeltaUstepdh != 0)
	    delete dDeltaUstepdh;  
	dDeltaUstepdh = new Vector(size);
	if (dDeltaUstepdh == 0 || dDeltaUstepdh->Size() != size) { 
	    opserr << "FATAL ArcLength::domainChanged() - ran out of memory for";
	    opserr << " dDeltaUstepdh Vector of size " << size << endln;
	    exit(-1);
	}
    }

    if (phat == 0 || phat->Size() != size) { 
	if (phat != 0)
	    delete phat;  
	phat = new Vector(size);
	if (phat == 0 || phat->Size() != size) { 
	    opserr << "FATAL ArcLength::domainChanged() - ran out of memory for";
	    opserr << " phat Vector of size " << size << endln;
	    exit(-1);
	}
    }    
 
   if (dphatdh == 0 || dphatdh->Size() != size) { 
      if (dphatdh != 0)
	 delete dphatdh;  
      dphatdh = new Vector(size);
      if (dphatdh == 0 || dphatdh->Size() != size) { 
	 opserr << "FATAL DisplacementControl::domainChanged() - ran out of memory for";
	 opserr << " dphatdh Vector of size " << size << endln;
	 exit(-1);
      }
   }    

 if (dUhatdh == 0 || dUhatdh->Size() != size) { 
      if (dUhatdh != 0)
	 delete dUhatdh;  
      dUhatdh = new Vector(size);
      if (dUhatdh == 0 || dUhatdh->Size() != size) { 
	 opserr << "FATAL DisplacementControl::domainChanged() - ran out of memory for";
	 opserr << " dUhatdh Vector of size " << size << endln;
	 exit(-1);
      }
   } 
  if (dUIJdh == 0 || dUIJdh->Size() != size) { 
      if (dUIJdh != 0)
	 delete dUIJdh;  
      dUIJdh = new Vector(size);
      if (dUIJdh == 0 || dUIJdh->Size() != size) { 
	 opserr << "FATAL DisplacementControl::domainChanged() - ran out of memory for";
	 opserr << " dUIJdh Vector of size " << size << endln;
	 exit(-1);
      }
   }
 if (Residual == 0 || Residual->Size() != size) { 
      if (Residual != 0)
	 delete Residual;  
      Residual = new Vector(size);
      if (Residual == 0 || Residual->Size() != size) { 
	 opserr << "FATAL DisplacementControl::domainChanged() - ran out of memory for";
	 opserr << " Residual Vector of size " << size << endln;
	 exit(-1);
      }
   } 


   if (sensU == 0 || sensU->Size() != size) { 
      if (sensU != 0)
	 delete sensU;  
      sensU = new Vector(size);
      if (sensU == 0 || sensU->Size() != size) { 
	 opserr << "FATAL DisplacementControl::domainChanged() - ran out of memory for";
	 opserr << " sensU Vector of size " << size << endln;
	 exit(-1);
      }
   } 



 Domain *theDomain=theModel->getDomainPtr();
   int numGrads = theDomain->getNumParameters();

   if (dLAMBDAdh == 0 || dLAMBDAdh->Size() != (numGrads)) { 
      if (dLAMBDAdh != 0 )  
	 delete dLAMBDAdh;
      dLAMBDAdh = new Vector(numGrads);
      if (dLAMBDAdh== 0 || dLAMBDAdh->Size() != (numGrads)) { 
	 opserr << "FATAL DisplacementControl::domainChanged() - ran out of memory for";
	 opserr << " dLAMBDAdh Vector of size " << numGrads << endln;
	 exit(-1);
      }
   } 

    // now we have to determine phat
    // do this by incrementing lambda by 1, applying load
    // and getting phat from unbalance.
    currentLambda = theModel->getCurrentDomainTime();
    currentLambda += 1.0;
    theModel->applyLoadDomain(currentLambda);    
    this->formUnbalance(); // NOTE: this assumes unbalance at last was 0
    (*phat) = theLinSOE->getB();
    currentLambda -= 1.0;
    theModel->setCurrentDomainTime(currentLambda);    
    

    // check there is a reference load
    int haveLoad = 0;
    for (int i=0; i<size; i++)
      if ( (*phat)(i) != 0.0 ) {
	haveLoad = 1;
	i = size;
      }

    if (haveLoad == 0) {
      opserr << "WARNING ArcLength::domainChanged() - zero reference load";
      return -1;
    }
  //  opserr<<" domainChanged: end"<<endln;
  return 0;  
}

int
ArcLength::sendSelf(int cTag,
		    Channel &theChannel)
{
 //  opserr<<"sendSelf: start"<<endln;

  Vector data(5);
  data(0) = arcLength2;
  data(1) = alpha2;
  data(2) = deltaLambdaStep;
  data(3) = currentLambda;
  data(4)  = signLastDeltaLambdaStep;

  if (theChannel.sendVector(this->getDbTag(), cTag, data) < 0) {
      opserr << "ArcLength::sendSelf() - failed to send the data\n";
      return -1;
  }
 // opserr<<"sendSelf: end"<<endln;
  return 0;
}


int
ArcLength::recvSelf(int cTag,
		    Channel &theChannel, FEM_ObjectBroker &theBroker)
{
 //  opserr<<"ArcLength:: recSelf: start"<<endln;

  Vector data(5);
  if (theChannel.recvVector(this->getDbTag(), cTag, data) < 0) {
      opserr << "ArcLength::sendSelf() - failed to send the data\n";
      return -1;
  }      

  // set the data
  arcLength2 = data(0);
  alpha2 = data(1);
  deltaLambdaStep = data(2);
  currentLambda = data(3);
  signLastDeltaLambdaStep = data(4);
 // opserr<<"recSelf: end"<<endln;

  return 0;
}

void
ArcLength::Print(OPS_Stream &s, int flag)
{
    AnalysisModel *theModel = this->getAnalysisModel();
    if (theModel != 0) {
	double cLambda = theModel->getCurrentDomainTime();
	s << "\t ArcLength - currentLambda: " << cLambda;
	s << "  arcLength: " << sqrt(arcLength2) <<  "  alpha: ";
	s << sqrt(alpha2) << endln;
    } else 
	s << "\t ArcLength - no associated AnalysisModel\n";
}


//////////////////////////////////// Sensitivity Begin/////////////////////////////
//Added by Abbas
//obtain the derivative of the tangent displacement (dUhatdh)
   void
ArcLength::formTangDispSensitivity(int gradNumber)
{
   AnalysisModel *theModel=this->getAnalysisModel();
   int size=theModel->getNumEqn();
   LinearSOE *theLinSOE = this->getLinearSOE(); 
// To get the structural stiffness Matrix
//...............................................................
dphatdh->Zero();
//static Matrix K(size,size);
//K.Zero();
//this->formTangent();
//opserr<<" before getK"<<endln;
//K=this->getK();
//opserr<<"after getK"<<endln;

/*
opserr<<"the tangent printed from the DisplacementControl.cpp"<<endln;
for(int i=0;i<size;i++)
{
for(int j=0;j<size;j++)
{
opserr<<K(i,j)<<"      ";

}
opserr<<endln;

}
*/
// ................................................................


//this->formTangentSensitivity(CURRENT_TANGENT);
//static Matrix dKdh(size,size);
//dKdh.Zero();
//dKdh=this->getdKdh(gradNumber);
// To print dKdh to the Screen 
//................................................................

//opserr<<"dKdh from the DisplacementControl.cpp"<<endln;
//for(int i=0;i<size;i++)
//{
//for(int j=0;j<size;j++)
//{
//opserr<<dKdh(i,j)<<"      ";

//}
//opserr<<endln;
//}

//form dKdh*Uft
//...............................................................

//dphatdh->addMatrixVector(1.0,dKdh,*deltaUhat,-1.0);
//call the tangent (K)
this->formTangent();
theLinSOE->setB(*dphatdh);
if(theLinSOE->solve()<0) {
opserr<<"SOE failed to obtained dUhatdh ";
exit(-1);
}
(*dUhatdh)=theLinSOE->getX();
//opserr<<"final dUhatdh is "<<*dUhatdh<<endln;

  
  
/*
   //call the tangent (K)
   theLinSOE->setB(*dphatdh);
   if(theLinSOE->solve()<0) {
      opserr<<"SOE failed to obtained dUhatdh ";
      exit(-1);
   }
   (*dUhatdh)=theLinSOE->getX();
*/

   // if the parameter is a load parameter.
   ////////////////////////////////////////////////////////
   // Loop through the loadPatterns and add the dPext/dh contributions

   static Vector oneDimVectorWithOne(1);
   oneDimVectorWithOne(0) = 1.0;
   static ID oneDimID(1);
   Node *aNode;
   DOF_Group *aDofGroup;

   int nodeNumber, dofNumber, relevantID, i, sizeRandomLoads, numRandomLoads;

   LoadPattern *loadPatternPtr;
   //AnalysisModel *theModel = this->getAnalysisModel();
   Domain *theDomain = theModel->getDomainPtr();
   LoadPatternIter &thePatterns = theDomain->getLoadPatterns();

   while((loadPatternPtr = thePatterns()) != 0) {
      const Vector &randomLoads = loadPatternPtr->getExternalForceSensitivity(gradNumber);
      sizeRandomLoads = randomLoads.Size();
      if (sizeRandomLoads == 1) {
	 // No random loads in this load pattern

      }
      else {
	 // opserr<<"there is sensitivity load parameter"<<endln;//Abbas.............
	 // Random loads: add contributions to the 'B' vector
	 numRandomLoads = (int)(sizeRandomLoads/2);
	 for (i=0; i<numRandomLoads*2; i=i+2) {
	    nodeNumber = (int)randomLoads(i);
	    dofNumber = (int)randomLoads(i+1);
	    aNode = theDomain->getNode(nodeNumber);
	    aDofGroup = aNode->getDOF_GroupPtr();
	    const ID &anID = aDofGroup->getID();
	    relevantID = anID(dofNumber-1);
	    oneDimID(0) = relevantID;
	    theLinSOE->addB(oneDimVectorWithOne, oneDimID);
	    (*dphatdh)=theLinSOE->getB();


	 }
      }
   }

   if(theLinSOE->solve()<0) {
      opserr<<"SOE failed to obtained dUhatdh ";
      exit(-1);
   }
//if(dphatdh !=0) {
 //  this->formTangent();
//   (*dUhatdh)=theLinSOE->getX();
//}

}

// form dLambda for each time step dLambda
   double 
ArcLength::formdLambdaDh(int gradNumber)
{
  
double dUhatTdUhat=((*deltaUhat)^(*deltaUhat));
double dUhatTUhatdh=(*deltaUhat)^(*dUhatdh);
//opserr<<"deltaUhat="<<*deltaUhat<<endln;
//opserr<<"dUhatdh= "<<*dUhatdh<<endln;
if(dLAMBDA==0.0 )
{
dlambda1dh=0.0;
//opserr<<"ArcLength: dLAMBDA=0 !!!!!!!!!!!!!!!!!!!"<<endln;
} else {
double ALPHA2=alpha2*alpha2;
double denomerator=pow((dUhatTdUhat+alpha2),2.0);
      dlambda1dh=signLastDeltaLambdaStep *1.0/(dLAMBDA)*(-arcLength2*dUhatTUhatdh/(denomerator));
    //  opserr<<"first dLambda1dh= "<<dlambda1dh<<endln;
 //dlambda1dh=1.0/(2.0)*dLAMBDA*(-sqrt(arcLength2)*dUhatTUhatdh/(dUhatTdUhat+alpha2));
//opserr<<"second dlambdadh: "<<dlambda1dh<<".......................,,,,,,,,,,,,////////////"<<endln;

 //dlambda1dh *=signLastDeltaLambdaStep;
//opserr<<"formdLAmbdaDh:   dLAMBDA= "<<dLAMBDA<<endln;
//opserr<<"denomerator= "<<denomerator<<endln;
//opserr<<"ArcLength2= "<<arcLength2<<endln;
//opserr<<"dUhatTUhatdh= "<<dUhatTUhatdh<<endln;
//opserr<<"dUhatTdUhat= "<<dUhatTdUhat<<endln;
}
//opserr<<"formdLambdaDh= "<<signLastDeltaLambdaStep<<endln;
   if(dLAMBDAdh != 0) {
      (*dLAMBDAdh)(gradNumber) = (*dLAMBDAdh)(gradNumber) + dlambda1dh;
  //    opserr<<"1st iteration "<<(*dLAMBDAdh)(gradNumber)<<endln;
      return (*dLAMBDAdh)(gradNumber);
   } else {
      return 0.0;
   }
}

// dLambdadh of the subsequent iterations (J>1)
   double 
ArcLength::getLambdaSensitivity(int gradNumber)
{

    //...................
 // determine the coeeficients of our quadratic equation
  //  double a = alpha2 + ((*deltaUhat)^(*deltaUhat));
 //   double b = alpha2*deltaLambdaStep 
  //    + ((*deltaUhat)^(*deltaUbar))
 //     + ((*deltaUstep)^(*deltaUhat));
 //   b *= 2.0;
 //   double c = 2*((*deltaUstep)^(*deltaUbar)) + ((*deltaUbar)^(*deltaUbar));
    // check for a solution to quadratic
  //  double b24ac = b*b - 4.0*a*c;
     if (b24ac < 0) {
      opserr << "ArcLength::update() - imaginary roots due to multiple instability";
      opserr << " directions - initial load increment was too large\n";
      opserr << "a: " << a << " b: " << b << " c: " << c << " b24ac: " << b24ac << endln;
      return -1;
    }			       
    double a2 = 2.0*a;
    if (a2 == 0.0) {
      opserr << "ArcLength::update() - zero denominator";
      opserr << " alpha was set to 0.0 and zero reference load\n";
      return -2;
    }			       
 double dAdh;
 double dBdh;
 double dCdh;
 dAdh =2.0*((*deltaUhat)^(*dUhatdh));
 dBdh = 2.0*(((*dUIJdh)^(*deltaUhat))+((*deltaUbar)^(*dUhatdh))+((*deltaUstep2)^(*dUhatdh))+((*dDeltaUstepdh)^(*deltaUhat))+alpha2*dDeltaLambdaStepdh);
 dCdh=2.0*(((*deltaUstep2)^(*dUIJdh))+((*dDeltaUstepdh)^(*deltaUbar))+((*deltaUbar)^(*dUIJdh)));//+2.0*((*deltaUstep)^(*dDeltaUstepdh));//+alpha2*dDeltaLambdaStepdh*((*phat)^(*phat)));

   
// opserr<<":a is "<<a<<endln;
    double sqrtb24ac = sqrt(b24ac);
    double dSqrtb24acdh=(2.0*b*dBdh-4.0*(a*dCdh+dAdh*c))/(2.0*sqrtb24ac);
    double dlambda1 = (-b + sqrtb24ac)/a2;
    double dlambdaj1dh=(a2*(-dBdh+dSqrtb24acdh)-((-b+sqrtb24ac)*2.0*dAdh))/(4.0*a*a);

    double dlambda2 = (-b - sqrtb24ac)/a2;
    double dlambdaj2dh= (a2*(-dBdh-dSqrtb24acdh)-((-b-sqrtb24ac)*2.0*dAdh))/(4.0*a*a);
   
    double val = (*deltaUhat)^(*deltaUstep2);
    double theta1 = ((*deltaUstep2)^(*deltaUstep2)) + ((*deltaUbar)^(*deltaUstep2));
    theta1 += dlambda1*val;

    double dTheta1dh=2*((*deltaUstep2)^(*dDeltaUstepdh))+((*deltaUbar)^(*dDeltaUstepdh))+((*dUIJdh)^(*deltaUstep2));
    //    double theta2 = theta1 + dlambda2*val;
      double dvaldh= ((*deltaUhat)^(*dDeltaUstepdh))+((*dUhatdh)^(*deltaUstep2));
     dTheta1dh +=dlambdaj1dh*val+dlambda1*dvaldh;
    // choose dLambda based on angle between incremental displacement before
    // and after this step -- want positive
    
      if (dTheta1dh > 0)
      dlambdaJdh = dlambdaj1dh;
    else
      dlambdaJdh = dlambdaj2dh;
  //  double dResultdh=a*2*dlambdaJdh*dLAMBDA2+dAdh*(dLAMBDA2*dLAMBDA2)+b*dlambdaJdh+dBdh*dLAMBDA2+dCdh;
//    opserr<<"dRdh= "<<dResultdh<<",,.,.,.,.,.,.,.,.,.,.,.,.,,.,.,.,.,.,.,.,.,///////"<<endln;
///////////////
   // determine delta U(i)
  ////////////////
   //  sensU->addVector(1.0,*dUhatdh,dLAMBDA2);
  (*deltaU) = (*deltaUbar);    
    deltaU->addVector(1.0, *deltaUhat,dLAMBDA2);
    
    // update dU and dlambda
    (*deltaUstep2) += *deltaU;

    dDeltaUstepdh->addVector(1.0,*dUhatdh,dLAMBDA2);
   dDeltaUstepdh->addVector(1.0,*deltaUhat,dlambdaJdh);
   (*dDeltaUstepdh) +=(*dUIJdh);
  



    dDeltaLambdaStepdh +=dlambdaJdh;

//opserr<<"Ok....."<<2.0*((*deltaUstep)^(*dDeltaUstepdh))<<"......//////////"<<endln;
//..............

  // dDeltaLambdaStepdh +=dlambdaJdh;
  // Now update Lambda_ij
   if(dLAMBDAdh !=0) {
      (*dLAMBDAdh)(gradNumber) = (*dLAMBDAdh)(gradNumber)+ dlambdaJdh;
    
      return (*dLAMBDAdh)(gradNumber);
   } else {
      
      return 0.0;
   }
}

 void
ArcLength::formResidualDispSensitivity( int gradNumber)
{
AnalysisModel *theModel=this->getAnalysisModel();
int size=theModel->getNumEqn();
//static Matrix dKdh(size,size);
//dKdh=this->getdKdh(gradNumber);
//dUIJdh-> addMatrixVector(1.0,dKdh,*deltaUbar,-1.0);
  
}



 int
ArcLength::formEleResidual(FE_Element* theEle)
{
   if(sensitivityFlag == 0) {  // no sensitivity
      this->StaticIntegrator::formEleResidual(theEle);
   } else {
      theEle->zeroResidual();
      theEle->addResistingForceSensitivity(gradNumber);
   }
   return 0;
}


 int
ArcLength::formIndependentSensitivityRHS()
{
   return 0;
}

  int
ArcLength::formSensitivityRHS(int passedGradNumber)
{
   sensitivityFlag = 1;

   // Set a couple of data members
   gradNumber = passedGradNumber;

   // get model
   AnalysisModel* theAnalysisModel = this->getAnalysisModel();
   LinearSOE* theSOE = this->getLinearSOE();

   // Loop through elements
   FE_Element *elePtr;
   FE_EleIter &theEles = theAnalysisModel->getFEs(); 

   while((elePtr = theEles()) != 0) {
      theSOE->addB(elePtr->getResidual(this) ,elePtr->getID()  );
   }

  
   (*Residual)=theSOE->getB();
//opserr<<"Residual= "<<*Residual<<endln;
   double CallDlambda1dh=(*dLAMBDAdh)(gradNumber);
  Residual->addVector(1.0,*phat, CallDlambda1dh ); //needed to calculate dLambdadh
 // Residual->addVector(1.0,*dphatdh, currentLambda ); //needed to calculate dLambdadh

/////////////
//int size=theAnalysisModel->getNumEqn();
//static Matrix dKdh(size,size);
//dKdh=this->getdKdh(gradNumber);
//Residual-> addMatrixVector(1.0,dKdh,*deltaUbar,-1.0);

///////

   theSOE->setB(*Residual);


   // Loop through the loadPatterns and add the dPext/dh contributions
   static Vector oneDimVectorWithOne(1);
   oneDimVectorWithOne(0) = 1.0;
   static ID oneDimID(1);

   Node *aNode;
   DOF_Group *aDofGroup;
   int nodeNumber, dofNumber, relevantID, i, sizeRandomLoads, numRandomLoads;
   LoadPattern *loadPatternPtr;
   Domain *theDomain = theAnalysisModel->getDomainPtr();
   LoadPatternIter &thePatterns = theDomain->getLoadPatterns();
   while((loadPatternPtr = thePatterns()) != 0) {
      const Vector &randomLoads = loadPatternPtr->getExternalForceSensitivity(gradNumber);
      sizeRandomLoads = randomLoads.Size();
      if (sizeRandomLoads == 1) {
	 // No random loads in this load pattern

      }
      else {

	 // Random loads: add contributions to the 'B' vector
	 numRandomLoads = (int)(sizeRandomLoads/2);
	 for (i=0; i<numRandomLoads*2; i=i+2) {
	    nodeNumber = (int)randomLoads(i);
	    dofNumber = (int)randomLoads(i+1);
	    aNode = theDomain->getNode(nodeNumber);
	    aDofGroup = aNode->getDOF_GroupPtr();
	    const ID &anID = aDofGroup->getID();
	    relevantID = anID(dofNumber-1);
	    oneDimID(0) = relevantID;
	    theSOE->addB(oneDimVectorWithOne, oneDimID);


	 }
      }

   }


   theSOE->setB(*Residual);




   //reset sensitivity flag
   sensitivityFlag=0;

   return 0;
}


   int
ArcLength::saveSensitivity(const Vector &v, int gradNum, int numGrads)
{
   AnalysisModel* theAnalysisModel = this->getAnalysisModel();

   DOF_GrpIter &theDOFGrps = theAnalysisModel->getDOFs();
   DOF_Group 	*dofPtr;

   while ( (dofPtr = theDOFGrps() ) != 0)  {
      //	dofPtr->saveSensitivity(v,0,0,gradNum,numGrads);
      dofPtr->saveDispSensitivity(v,gradNum,numGrads);

   }

   return 0;
}

   int
ArcLength::saveLambdaSensitivity(double dlambdadh, int gradNum, int numGrads)
{
   AnalysisModel* theAnalysisModel = this->getAnalysisModel();
   Domain *theDomain = theAnalysisModel->getDomainPtr();

   LoadPattern *lpPtr;
   LoadPatternIter &thePatterns = theDomain->getLoadPatterns();
  
   while ( (lpPtr = thePatterns() ) != 0)
            lpPtr->saveLoadFactorSensitivity(dlambdadh, gradNum, numGrads);

   return 0;
}

   int 
ArcLength::commitSensitivity(int gradNum, int numGrads)
{

   AnalysisModel* theAnalysisModel = this->getAnalysisModel();


   // Loop through the FE_Elements and set unconditional sensitivities
   FE_Element *elePtr;
   FE_EleIter &theEles = theAnalysisModel->getFEs();    
   while((elePtr = theEles()) != 0) {
      elePtr->commitSensitivity(gradNum, numGrads);
   }
   return 0;
}



   bool 
ArcLength::computeSensitivityAtEachIteration()
{

   return true;     
}


   int 
ArcLength::computeSensitivities(void)
{
   LinearSOE *theSOE = this->getLinearSOE();
//opserr<<"computeSensitivities : start"<<endln;
   // Zero out the old right-hand side of the SOE
   theSOE->zeroB();
   if (this == 0) {
      opserr << "ERROR SensitivityAlgorithm::computeSensitivities() -";
      opserr << "the SensitivityIntegrator is NULL\n";
      return -1;
   }


   // Form the part of the RHS which are indepent of parameter
   this->formIndependentSensitivityRHS();

   AnalysisModel *theModel = this->getAnalysisModel();   
   Domain *theDomain=theModel->getDomainPtr();//Abbas
   ParameterIter &paramIter = theDomain->getParameters();
   Parameter *theParam;

   // De-activate all parameters
   while ((theParam = paramIter()) != 0)
      theParam->activate(false);

   // Now, compute sensitivity wrt each parameter
   int  numGrads = theDomain->getNumParameters();
   paramIter = theDomain->getParameters();


   while ((theParam = paramIter()) != 0 ) {
      // Activate this parameter
      theParam->activate(true);

      // Zero the RHS vector
      theSOE->zeroB();

      // Get the grad index for this parameter
      int gradIndex = theParam->getGradIndex();
     
      // Form the RHS

        this->formTangDispSensitivity(gradIndex);
      this->formSensitivityRHS(gradIndex);

      this->formTangent();
      theSOE->solve();
      *dUIJdh=theSOE->getX();// sensitivity of the residual displacement

     //    opserr<<"deltaUbar is "<<*deltaUbar<<endln;   
     // this->formTangDispSensitivity(gradIndex);
      double dlamdh = this->getLambdaSensitivity(gradIndex);
//opserr<<" dLAMBDAdh======================================================="<<dlamdh<<endln;

///////////////////
      // To obtain the response sensitivity 
      theSOE->setB(*Residual);
      theSOE->solve();
      (*sensU) = theSOE->getX();
    //  this->formResidualDispSensitivity(gradNumber);
// (*sensU)=(*dUIJdh);
     //..............
   //  (*sensU) +=(*dUIJdh);
  
  //  sensU->addVector(1.0,*dUhatdh,dLAMBDA2);
 // sensU->addVector(1.0,*deltaUhat,dlamdh);
 //(*sensU) =(*dUIJdh);

  //(*dDeltaUstepdh) +=(*dUIJdh);

   
    //  (*sensU) +=(*dDeltaUstepdh);

    //.......
      // Save sensitivity to nodes
      this->saveSensitivity( (*sensU), gradIndex, numGrads );
     
      
      this->saveLambdaSensitivity(dlamdh, gradIndex, numGrads);

      // Commit unconditional history variables (also for elastic problems; strain sens may be needed anyway)
      this->commitSensitivity(gradIndex, numGrads);
      // De-activate this parameter for next sensitivity calc
    //  opserr<<"dUIJdh "<<*dUIJdh<<"...................Abbas."<<endln;
    //  opserr<<"deltaUbar "<<*deltaUbar<<"ABBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBASSSSSSSSSS"<<endln;
      theParam->activate(false);

    //  theSOE->zeroB();//reset the SOE to zero ;Abbas

   } 
   // end of if statment to be run only one time during the iteration process.
//opserr<<"computeSensitivities : end"<<endln;
   return 0;
}

