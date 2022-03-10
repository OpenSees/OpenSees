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
                                                                        
// $Revision: 1.9 $
// $Date: 2007-04-02 23:42:26 $
// $Source: /usr/local/cvs/OpenSees/SRC/analysis/integrator/MinUnbalDispNorm.cpp,v $
                                                                        
                                                                        
// File: ~/analysis/integrator/MinUnbalDispNorm.C
// 
// Written: fmk 
// Created: 07/99
// Revision: A
//
// What: "@(#) MinUnbalDispNorm.C, revA"


#include <MinUnbalDispNorm.h>
#include <AnalysisModel.h>
#include <LinearSOE.h>
#include <Vector.h>
#include <Channel.h>
#include <math.h>
#include <Domain.h>
#include<Node.h>
#include<DOF_Group.h>
#include<DOF_GrpIter.h>
#include<ID.h>
#include <stdlib.h>
#include<FE_Element.h>
#include<FE_EleIter.h>
#include<LoadPattern.h>
#include<LoadPatternIter.h>
#include<Parameter.h>
#include<ParameterIter.h>
#include<EquiSolnAlgo.h>
#include<TaggedObjectStorage.h>
#include <elementAPI.h>
#include <Matrix.h>
void* OPS_MinUnbalDispNorm()
{
    double lambda11, minlambda, maxlambda;
    int numIter;
    if (OPS_GetNumRemainingInputArgs() < 1) {
	opserr << "WARNING integrator MinUnbalDispNorm lambda11 <Jd minLambda1j maxLambda1j>\n";
	return 0;
    }

    int numdata = 1;
    if (OPS_GetDoubleInput(&numdata, &lambda11) < 0) {
	opserr << "WARNING integrator MinUnbalDispNorm invalid lambda11\n";
	return 0;
    }

    if (OPS_GetNumRemainingInputArgs() >= 3) {
	if (OPS_GetIntInput(&numdata, &numIter) < 0) {
	    opserr << "WARNING integrator MinUnbalDispNorm invalid numIter\n";
	    return 0;
	}
	if (OPS_GetDoubleInput(&numdata, &minlambda) < 0) {
	    opserr << "WARNING integrator MinUnbalDispNorm invalid minlambda\n";
	    return 0;
	}
	if (OPS_GetDoubleInput(&numdata, &maxlambda) < 0) {
	    opserr << "WARNING integrator MinUnbalDispNorm invalid maxlambda\n";
	    return 0;
	}
    }
    else {
	minlambda = lambda11;
	maxlambda = lambda11;
	numIter = 1;
    }

    int signFirstStepMethod = SIGN_LAST_STEP;
    if (OPS_GetNumRemainingInputArgs() > 0) {
	const char* flag = OPS_GetString();
	if ((strcmp(flag,"-determinant") == 0) ||
	    (strcmp(flag,"-det") == 0)) {
	    signFirstStepMethod = CHANGE_DETERMINANT;
	}
    }

    return new MinUnbalDispNorm(lambda11,numIter,minlambda,maxlambda,signFirstStepMethod);

}

MinUnbalDispNorm::MinUnbalDispNorm(double lambda1, int specNumIter,
		     double min, double max, int signFirstStep)
:StaticIntegrator(INTEGRATOR_TAGS_MinUnbalDispNorm),
 dLambda1LastStep(lambda1), 
 specNumIncrStep(specNumIter), numIncrLastStep(specNumIter),
 deltaUhat(0), deltaUbar(0), deltaU(0), deltaUstep(0), dUhatdh(0), dLambdaj(0.0),
 phat(0), deltaLambdaStep(0.0), currentLambda(0.0), dLambdaStepDh(0.0),dUIJdh(0),Dlambdadh(0.0),dphatdh(0),Residual2(0),
 signLastDeltaLambdaStep(1), sensitivityFlag(0),Residual(0), dlambdadh(0.0),dLambda(0.0), sensU(0),d_deltaU_dh(0),gradNumber(0),dLAMBDAdh(0),
 dLambda1min(min), dLambda1max(max), signLastDeterminant(1), signFirstStepMethod(signFirstStep)
{
  // to avoid divide-by-zero error on first update() ensure numIncr != 0
  if (specNumIncrStep == 0) {
    opserr << "WARNING LoadControl::LoadControl() - numIncr set to 0, 1 assumed\n";
    specNumIncrStep = 1.0;
    numIncrLastStep = 1.0;
  }
}

MinUnbalDispNorm::~MinUnbalDispNorm()
{
    // delete any vector object created
    if (deltaUhat != 0)
	delete deltaUhat;
    if (deltaU != 0)
	delete deltaU;
    if (deltaUstep != 0)
	delete deltaUstep;
    if (deltaUbar != 0)
	delete deltaUbar;
    if (phat != 0)
	delete phat;
    if(dUhatdh !=0)
     delete dUhatdh;
   if(dUIJdh !=0)
     delete dUIJdh; 
   if(Residual !=0)
     delete Residual;
   if(sensU !=0)
     delete sensU;
   if(Residual2 !=0)
     delete Residual2;
   if(dLAMBDAdh !=0) 
     delete dLAMBDAdh;
   if(dphatdh !=0)
      delete dphatdh;

   dLAMBDAdh=0;
   dUhatdh=0;


}

int
MinUnbalDispNorm::newStep(void)
{
//opserr<<"   New Step................................"<<endln;    // get pointers to AnalysisModel and LinearSOE
    AnalysisModel *theModel = this->getAnalysisModel();
    LinearSOE *theLinSOE = this->getLinearSOE();    
    if (theModel == 0 || theLinSOE == 0) {
	opserr << "WARNING MinUnbalDispNorm::newStep() ";
	opserr << "No AnalysisModel or LinearSOE has been set\n";
	return -1;
    }

    // get the current load factor
    currentLambda = theModel->getCurrentDomainTime();

//opserr<<" NewStep=      "<<*phat<<endln;
    // determine dUhat
    this->formTangent();
    theLinSOE->setB(*phat);
    if (theLinSOE->solve() < 0) {
      opserr << "MinUnbalanceDispNorm::newStep(void) - failed in solver\n";
      return -1;
    }
    (*deltaUhat) = theLinSOE->getX();
    Vector &dUhat = *deltaUhat;

    // determine delta lambda(1) == dlambda
    double factor = specNumIncrStep/numIncrLastStep;
    double dLambda = dLambda1LastStep*factor;

    // check aaint min and max values specified in constructor
    if (dLambda < dLambda1min)
      dLambda = dLambda1min;
    else if (dLambda > dLambda1max)
      dLambda = dLambda1max;

    dLambda1LastStep = dLambda;


    if (signFirstStepMethod == SIGN_LAST_STEP) {
      if (deltaLambdaStep < 0)
	signLastDeltaLambdaStep = -1;
      else
	signLastDeltaLambdaStep = +1;
      
      dLambda *= signLastDeltaLambdaStep; // base sign of load change
                                        // on what was happening last step
    } else {
    
      double det = theLinSOE->getDeterminant();
      int signDeterminant = 1;
      if (det < 0)
	signDeterminant = -1;
    
      dLambda *= signDeterminant * signLastDeterminant;
    
      signLastDeterminant = signDeterminant;
    }

    /*
    double work = (*phat)^(dUhat);
    int signCurrentWork = 1;
    if (work < 0) signCurrentWork = -1;

    if (signCurrentWork != signLastDeltaStep)
    */

    deltaLambdaStep = dLambda;
    currentLambda += dLambda;
    numIncrLastStep = 0;

    // determine delta U(1) == dU
    (*deltaU) = dUhat;
    (*deltaU) *= dLambda;
    (*deltaUstep) = (*deltaU);

//////////////////


   ////////////////Abbas////////////////////////////

  if(this->activateSensitivity()==true) { 
    Domain *theDomain=theModel->getDomainPtr();
    ParameterIter &paramIter = theDomain->getParameters();
    Parameter *theParam;

    // De-activate all parameters
     
    // Now, compute sensitivity wrt each parameter
    int numGrads = theDomain->getNumParameters();
    
    while ((theParam = paramIter()) != 0)
      theParam->activate(false);
    
    paramIter = theDomain->getParameters();
    while ((theParam = paramIter()) != 0) {
      // Activate this parameter
      theParam->activate(true);
      // Get the grad index for this parameter
      gradNumber = theParam->getGradIndex();
      
      this->formTangDispSensitivity(dUhatdh,gradNumber);
//opserr<<"NewStep: dUhatdh=   "<<*dUhatdh<<endln;
//opserr<<"NewStep: dLambda= "<< dLambda<<endln;
      this->formdLambdaDh(gradNumber);

      
           sensU->addVector(1.0, *dUhatdh ,dLambda);
 //opserr<<"NewStep:    dUfdh...= "<<*sensU<<endln;

      theParam->activate(false);
 //     opserr<<"NewStep:  Uft= "<<*deltaUhat<<endln;
    } 
  }
  ///////////////Abbas/////////////////////////////




/////////////////    




    // update model with delta lambda and delta U
    theModel->incrDisp(*deltaU);    
    theModel->applyLoadDomain(currentLambda);    
    if (theModel->updateDomain() < 0) {
      opserr << "MinUnbalDispNorm::newStep - model failed to update for new dU\n";
      return -1;
    }

    return 0;
}

int
MinUnbalDispNorm::update(const Vector &dU)
{
 //  opserr<<"Update Function.............................."<<endln;
 
  
   AnalysisModel *theModel = this->getAnalysisModel();
    LinearSOE *theLinSOE = this->getLinearSOE();    
    if (theModel == 0 || theLinSOE == 0) {
	opserr << "WARNING MinUnbalDispNorm::update() ";
	opserr << "No AnalysisModel or LinearSOE has been set\n";
	return -1;
    }

    (*deltaUbar) = dU; // have to do this as the SOE is gonna change

    // determine dUhat    
    theLinSOE->setB(*phat);
    theLinSOE->solve();
    (*deltaUhat) = theLinSOE->getX();    

    // determine delta lambda(i)
    double a = (*deltaUhat)^(*deltaUbar);
    double b = (*deltaUhat)^(*deltaUhat);
    if (b == 0) {
      opserr << "MinUnbalDispNorm::update() - zero denominator\n";
      return -1;
    }

    double dLambda = -a/b;
    dLambdaj=dLambda;//Abbas
   //opserr<<"update function.... "<< " dLambda= "<<dLambda<<" deltaUhat= "<<(*deltaUhat)<<" deltaUbar: "<<*deltaUbar<<endln; 
    // determine delta U(i)
    (*deltaU) = (*deltaUbar);    
    deltaU->addVector(1.0, *deltaUhat,dLambda);
    
    // update dU and dlambda
    (*deltaUstep) += *deltaU;
    deltaLambdaStep += dLambda;
    currentLambda += dLambda;
 //   opserr<<"dLambda= "<<dLambda<<endln;
//opserr<<"CURRENTLAMBDA= ..................................................../////////......."<<currentLambda<<endln;
    // update the model
    theModel->incrDisp(*deltaU);    
    theModel->applyLoadDomain(currentLambda);    

    if (theModel->updateDomain() < 0) {
      opserr << "MinUnbalDispNorm::update - model failed to update for new dU\n";
      return -1;
    }
    
    // set the X soln in linearSOE to be deltaU for convergence Test
    theLinSOE->setX(*deltaU);

    numIncrLastStep++;

// opserr<<" UpdateFunction.....Uft= "<<*deltaUhat<<endln;
// opserr<<"UpdateFunction.... deltaUbar= "<<*deltaUbar<<endln;

    return 0;
}



int 
MinUnbalDispNorm::domainChanged(void)
{
    // we first create the Vectors needed
    AnalysisModel *theModel = this->getAnalysisModel();
    LinearSOE *theLinSOE = this->getLinearSOE();    
    if (theModel == 0 || theLinSOE == 0) {
	opserr << "WARNING MinUnbalDispNorm::update() ";
	opserr << "No AnalysisModel or LinearSOE has been set\n";
	return -1;
    }    
    int size = theModel->getNumEqn(); // ask model in case N+1 space

    if (deltaUhat == 0 || deltaUhat->Size() != size) { // create new Vector
	if (deltaUhat != 0)
	    delete deltaUhat;   // delete the old
	deltaUhat = new Vector(size);
	if (deltaUhat == 0 || deltaUhat->Size() != size) { // check got it
	    opserr << "FATAL MinUnbalDispNorm::domainChanged() - ran out of memory for";
	    opserr << " deltaUhat Vector of size " << size << endln;
	    exit(-1);
	}
    }

    if (deltaUbar == 0 || deltaUbar->Size() != size) { // create new Vector
	if (deltaUbar != 0)
	    delete deltaUbar;   // delete the old
	deltaUbar = new Vector(size);
	if (deltaUbar == 0 || deltaUbar->Size() != size) { // check got it
	    opserr << "FATAL MinUnbalDispNorm::domainChanged() - ran out of memory for";
	    opserr << " deltaUbar Vector of size " << size << endln;
	    exit(-1);
	}
    }

    
    if (deltaU == 0 || deltaU->Size() != size) { // create new Vector
	if (deltaU != 0)
	    delete deltaU;   // delete the old
	deltaU = new Vector(size);
	if (deltaU == 0 || deltaU->Size() != size) { // check got it
	    opserr << "FATAL MinUnbalDispNorm::domainChanged() - ran out of memory for";
	    opserr << " deltaU Vector of size " << size << endln;
	    exit(-1);
	}
    }

    if (deltaUstep == 0 || deltaUstep->Size() != size) { 
	if (deltaUstep != 0)
	    delete deltaUstep;  
	deltaUstep = new Vector(size);
	if (deltaUstep == 0 || deltaUstep->Size() != size) { 
	    opserr << "FATAL MinUnbalDispNorm::domainChanged() - ran out of memory for";
	    opserr << " deltaUstep Vector of size " << size << endln;
	    exit(-1);
	}
    }

    if (phat == 0 || phat->Size() != size) { 
	if (phat != 0)
	    delete phat;  
	phat = new Vector(size);
	if (phat == 0 || phat->Size() != size) { 
	    opserr << "FATAL MinUnbalDispNorm::domainChanged() - ran out of memory for";
	    opserr << " phat Vector of size " << size << endln;
	    exit(-1);
	}
    }    
  if (dphatdh == 0 || dphatdh->Size() != size) { 
     if (dphatdh != 0)
       delete dphatdh;  
     dphatdh = new Vector(size);
     if (dphatdh == 0 || dphatdh->Size() != size) { 
       opserr << "FATAL MinUnbalDispNorm::domainChanged() - ran out of memory for";
       opserr << " dphatdh Vector of size " << size << endln;
       exit(-1);
     }
   }    


   if (dUhatdh == 0 || dUhatdh->Size() != size) { 
     if (dUhatdh != 0)
       delete dUhatdh;  
     dUhatdh = new Vector(size);
     if (dUhatdh == 0 || dUhatdh->Size() != size) { 
       opserr << "FATAL MinUnbalDisporm::domainChanged() - ran out of memory for";
       opserr << " dUhatdh Vector of size " << size << endln;
       exit(-1);
     }
   } 
   
   if (dUIJdh == 0 || dUIJdh->Size() != size) { 
      if (dUIJdh != 0)
	 delete dUIJdh;  
      dUIJdh = new Vector(size);
      if (dUIJdh == 0 || dUIJdh->Size() != size) { 
	 opserr << "FATAL MinUnbalDispNorm::domainChanged() - ran out of memory for";
	 opserr << " dUIJdh Vector of size " << size << endln;
	 exit(-1);
      }
   }

   if (Residual == 0 || Residual->Size() != size) { 
      if (Residual != 0)
	 delete Residual;  
      Residual = new Vector(size);
      if (Residual == 0 || Residual->Size() != size) { 
	 opserr << "FATAL MinUnbalDispNorm::domainChanged() - ran out of memory for";
	 opserr << " Residual Vector of size " << size << endln;
	 exit(-1);
      }
   } 
   
   if (Residual2 == 0 || Residual2->Size() != size) { 
     if (Residual2 != 0)
       delete Residual2;  
     Residual2 = new Vector(size);
     if (Residual2== 0 || Residual2->Size() != size) { 
       opserr << "FATAL MinUnbalDispNorm::domainChanged() - ran out of memory for";
       opserr << " N Vector of size " << size << endln;
       exit(-1);
     }
   } 

   if (sensU == 0 || sensU->Size() != size) { 
      if (sensU != 0)
	 delete sensU;  
      sensU = new Vector(size);
      if (sensU == 0 || sensU->Size() != size) { 
	 opserr << "FATAL MinUnbalDispNorm::domainChanged() - ran out of memory for";
	 opserr << " sensU Vector of size " << size << endln;
	 exit(-1);
      }
   } 

   Domain *theDomain=theModel->getDomainPtr();//Abbas
   int numGrads = theDomain->getNumParameters();

   if (dLAMBDAdh == 0 || dLAMBDAdh->Size() != (numGrads)) { 
     if (dLAMBDAdh != 0 )  
       delete dLAMBDAdh;
     dLAMBDAdh = new Vector(numGrads);
     if (dLAMBDAdh== 0 || dLAMBDAdh->Size() != (numGrads)) { 
       opserr << "FATAL MinUnbalDispNorm::domainChanged() - ran out of memory for";
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

    
    
    return 0;
}

int
MinUnbalDispNorm::sendSelf(int cTag,
		    Channel &theChannel)
{
  Vector data(8);
  data(0) = dLambda1LastStep;
  data(1) = specNumIncrStep;
  data(2) = numIncrLastStep;
  data(3) = deltaLambdaStep;
  data(4) = currentLambda;
  if (signLastDeltaLambdaStep == 1)
    data(5)  = 1.0;
  else
    data(5) = 0.0;
  data(6) = dLambda1min;
  data(7) = dLambda1max;

  if (theChannel.sendVector(this->getDbTag(), cTag, data) < 0) {
      opserr << "MinUnbalDispNorm::sendSelf() - failed to send the data\n";
      return -1;
  }
  return 0;
}


int
MinUnbalDispNorm::recvSelf(int cTag,
		    Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  Vector data(8);
  if (theChannel.recvVector(this->getDbTag(), cTag, data) < 0) {
      opserr << "MinUnbalDispNorm::sendSelf() - failed to send the data\n";
      return -1;
  }      

  // set the data

  dLambda1LastStep = data(0);
  specNumIncrStep = data(1);
  numIncrLastStep = data(2);
  deltaLambdaStep = data(3);
  currentLambda = data(4);
  if (data(5)== 1.0)
    signLastDeltaLambdaStep = 1;
  else
    signLastDeltaLambdaStep = -1;
  dLambda1min = data(6);
  dLambda1max = data(7);

  return 0;
}

void
MinUnbalDispNorm::Print(OPS_Stream &s, int flag)
{
    AnalysisModel *theModel = this->getAnalysisModel();
    if (theModel != 0) {
	double cLambda = theModel->getCurrentDomainTime();
	s << "\t MinUnbalDispNorm - currentLambda: " << cLambda;
    } else 
	s << "\t MinUnbalDispNorm - no associated AnalysisModel\n";
}


///////////////////////////Sensitivity Begin///////////////////////

////obtain the derivative of the tangent displacement (dUhatdh)
   Vector *
MinUnbalDispNorm::formTangDispSensitivity(Vector *dUhatdh,int gradNumber)
{
   LinearSOE *theLinSOE = this->getLinearSOE(); 
   dUhatdh->Zero();
   dphatdh->Zero();

   //call the tangent (K)
   this->formTangent();
   theLinSOE->setB(*dphatdh);
   if(theLinSOE->solve()<0) {
      opserr<<"SOE failed to obtained dUhatdh ";
      exit(-1);
   }
   (*dUhatdh)=theLinSOE->getX();
 //  opserr<<"final dUhatdh is "<<*dUhatdh<<endln;


   
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
   AnalysisModel *theModel = this->getAnalysisModel();   
   Domain *theDomain = theModel->getDomainPtr();
   LoadPatternIter &thePatterns = theDomain->getLoadPatterns();
  
while((loadPatternPtr = thePatterns()) != 0) {
  // opserr<<"Inside While loop"<<endln;
     const Vector &randomLoads = loadPatternPtr->getExternalForceSensitivity(gradNumber);
      sizeRandomLoads = randomLoads.Size();
      if (sizeRandomLoads == 1) {
	 // No random loads in this load pattern

      }
      else {
//	 opserr<<"there is sensitivity load parameter"<<endln;//Abbas.............
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
// dphatdh->addMatrixVector(1.0,dKdh,*deltaUhat,1.0);


	 }
      }
   }

 if(theLinSOE->solve()<0) {
   opserr<<"SOE failed to obtained dUhatdh ";
   exit(-1);
 }

    (*dUhatdh)=theLinSOE->getX();

// opserr<<"MinUnbalDispNorm::dUhatdh is  "<<*dUhatdh<<endln;//

/////////////////////////////////////////////////////////
   return dUhatdh;
}

// form dLambda for each time step dLambda
double 
MinUnbalDispNorm::formdLambdaDh(int gradNumber)
{
  // Here dLambda1dh=0 because its constant.
  //return 0.0;
  if(dLAMBDAdh != 0)
    return (*dLAMBDAdh)(gradNumber);
}


 // dLambdadh of the subsequent iterations (J>1)
double 
MinUnbalDispNorm::getLambdaSensitivity(int gradNumber)
{

   //opserr<<"getLambdaSensitivityFunction:   deltaUhat= "<<(*deltaUhat)<<" Pref= "<<*phat<<endln; 
 
 //  Vector &dufRdh=*dUIJdh;// component of the dUfrDh: derivative of the residual displacement
 

  //  opserr<<" LambdaFunction.....Uft= "<<*deltaUhat<<endln;
// opserr<<"LambdaFunation.... deltaUbar= "<<*deltaUbar<<endln;
    double temp= (*deltaUhat)^(*deltaUhat);
   double denomerator= pow(temp,2.0);
   double a= (*deltaUhat)^(*dUIJdh);
   double b= (*dUhatdh)^(*deltaUbar);
   double c= (*deltaUhat)^(*deltaUbar);
   double d= (*deltaUhat)^(*dUhatdh);
   double Numerator= -(temp*(a+b)-(c*2.0*d));
     Dlambdadh=Numerator/denomerator;  //   

 //opserr<<"dLambdadh=    "<<Dlambdadh<<" a= "<<a<< " b= "<<b <<" c= "<<c << " d= "<<d<< "temp= "<<temp<< " denomera= "<<denomerator<<endln;
// opserr<<"Numerator= "<<Numerator<<endln;

   // Now update Lambda_ij
   if(dLAMBDAdh !=0) {
     (*dLAMBDAdh)(gradNumber) = (*dLAMBDAdh)(gradNumber)+ Dlambdadh;
     //opserr << gradNumber << ' ' << deltaLambdaStep << ' ' << Dlambdadh << endln;
     return (*dLAMBDAdh)(gradNumber);
   } else {
     return 0.0;
   }
}



int
MinUnbalDispNorm::formEleResidual(FE_Element* theEle)
{
   if(sensitivityFlag == 0) {  // no sensitivity
     this->StaticIntegrator::formEleResidual(theEle);
   } else {
  //    opserr<<"formEleResidual:    else condition lllllllllllllllllllllllllllllllllllllllllllllllllllll"<<endln;
     theEle->zeroResidual();
     theEle->addResistingForceSensitivity(gradNumber);
   }
   return 0;
}

int
MinUnbalDispNorm::formIndependentSensitivityRHS()
{
   return 0;
}


   
int
MinUnbalDispNorm::formSensitivityRHS(int passedGradNumber)
{

  sensitivityFlag = 1;
  //this->Activatesensitivity();//Abbas
  // Set a couple of data members
  gradNumber = passedGradNumber;
//opserr<<" RHS..."<<"gradNumber= "<<gradNumber<<endln;
  // get model
  AnalysisModel* theAnalysisModel = this->getAnalysisModel();
  LinearSOE* theSOE = this->getLinearSOE();
  //theSOE->zeroB();

  // Loop through elements
  FE_Element *elePtr;
  FE_EleIter &theEles = theAnalysisModel->getFEs(); 

  while((elePtr = theEles()) != 0) {
    theSOE->addB(elePtr->getResidual(this) ,elePtr->getID()  );
  }

  //  opserr<<" dKdh*deltaUbar after while loop is "<<*dUIJdh<<endln;
  (*Residual)=theSOE->getB();
//opserr<<"RHS.....Residual= ..... "<<*Residual<<endln;
//////////////////////////
 int size=theAnalysisModel->getNumEqn();
  Matrix dKdh(size,size);
  dKdh.Zero();
  //dKdh=this->getdKdh(gradNumber);
    
//////////////////////////
/*
  opserr<<"RHS sensitivity function:    dKdh= "<<endln;
 for(int i=0;i<size;i++)
 {
 for(int j=0; j<size;j++)
 {
 opserr<<dKdh(i,j)<<"   ";

 
  }
 opserr<<endln;
 
 } 
*/ 
  



//////////////////
   //call the tangent (K)
 //  this->formTangent();


//   Residual->addMatrixVector(1.0,dKdh,*deltaUbar,-1.0);

//opserr<<"RHS.....Residual22222222= ..... "<<*Residual<<endln;

///////////////////////

  // double CallDlambda1dh=this->getLambdaSensitivity(gradNumber);
     
  double CallDlambda1dh=(*dLAMBDAdh)(gradNumber);
  // opserr<<"RHS:: dphatdh is ....................................//"<<*dphatdh<<endln;
  Residual->addVector(1.0,*phat, CallDlambda1dh ); 
//opserr<<"RHS..............Pref= "<<*phat<<endln;
//opserr<<"RHS............ ...........................= "<<CallDlambda1dh<<endln;
// opserr<<"RHS.....Residual333333333333= ..... "<<*Residual<<endln;
  Residual->addVector(1.0,*dphatdh,currentLambda);

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
      	 //opserr<<"No sensitivity Load Parameter is involved"<<endln;
    }
    else {
      //	 opserr<<"there is sensitivity load parameter"<<endln;//Abbas.............
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
    //  (*Residual) =theSOE->getB();
  }
  
  
  theSOE->setB(*Residual);
  
    
  
  //reset sensitivity flag
  sensitivityFlag=0;
  //  opserr<<"RHS :Ends"<<endln;
  return 0;
}


int
MinUnbalDispNorm::saveSensitivity(const Vector &v, int gradNum, int numGrads)
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
MinUnbalDispNorm::saveLambdaSensitivity(double dlambdadh, int gradNum, int numGrads)
{
   AnalysisModel* theAnalysisModel = this->getAnalysisModel();
   Domain *theDomain = theAnalysisModel->getDomainPtr();

   LoadPattern *lpPtr;
   LoadPatternIter &thePatterns = theDomain->getLoadPatterns();
   while ( (lpPtr = thePatterns() ) != 0)
     lpPtr->saveLoadFactorSensitivity(dlambdadh, gradNum, numGrads);
//opserr<<"saveLambdaSensitivity....."<<dlambdadh<<endln;
   return 0;
}

   int 
MinUnbalDispNorm::commitSensitivity(int gradNum, int numGrads)
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
MinUnbalDispNorm::computeSensitivityAtEachIteration()
{

   return true;     
}


   int 
MinUnbalDispNorm::computeSensitivities(void)
{

  // opserr<<"Compute Sensitivity................////////////////////"<<endln;
  LinearSOE *theSOE = this->getLinearSOE();
  
  // Zero out the old right-hand side of the SOE
  theSOE->zeroB();

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

 
  while ((theParam = paramIter()) != 0) {
    // Activate this parameter
    theParam->activate(true);
    
    // Zero the RHS vector
    theSOE->zeroB();
    
    // Get the grad index for this parameter
    int gradIndex = theParam->getGradIndex();
               // this->formResidualDispSensitivity( dUIJdh, gradIndex);//./.
// opserr<<"B.........BBBBBBBBBBB"<<endln;
    // Form the RHS

//	opserr<<"CCCCCCCCCCCCC"<<endln;
  //   opserr<<"DDDDDDDD"<<endln;

    this->formSensitivityRHS(gradIndex);




         this->formTangent();
	 //theSOE->setX(*dUIJdh);
    theSOE->solve();
    *dUIJdh=theSOE->getX();// sensitivity of the residual displacement

    this->formTangDispSensitivity(dUhatdh,gradIndex);
    double dlamdh = this->getLambdaSensitivity(gradIndex);
    //double dlamdh=(*dLAMBDAdh)(gradIndex);
//opserr<<" EEEEEEEEEEEEEEEE  dLambdaijdh= "<<dlamdh<<endln;

    //  opserr<<"computeSensitivities: dUfRdh is  "<<*dUIJdh<<endln;//
       //   this->formTangDispSensitivity(dUhatdh,gradIndex);
   
    // To obtain the response sensitivity 
    theSOE->setB(*Residual);





    theSOE->solve();
    
//  opserr<<"FFFFFFFFFFFFFFFFFFF"<<endln;

  //  dUIJdh->addVector(1.0, *dUhatdh,dLambdaj);
//    dUIJdh->addVector(1.0, *deltaUhat, Dlambdadh);
  //  opserr<<" GGGGGGGGGGGGGGGGGG"<<endln;
//(*sensU)=(*dUIJdh); //theSOE->getX();
  //sensU->addVector(1.0,*dUIJdh,1.0);
    (*sensU) = theSOE->getX();//.................................
 // opserr<<"JJJJJJJJJJJJJ"<<endln; 
    //(*sensU) =(*dUIJdh);


   // opserr<<" gradNumber= "<<gradIndex<<endln;
    // Save sensitivity to nodes
    this->saveSensitivity( (*sensU), gradIndex, numGrads );
    this->saveLambdaSensitivity(dlamdh, gradIndex, numGrads);
    //dUIJdh->Zero();
    // Commit unconditional history variables (also for elastic problems; strain sens may be needed anyway)
    this->commitSensitivity(gradIndex, numGrads);
    // De-activate this parameter for next sensitivity calc
    theParam->activate(false);
//opserr<<"KKKKKKKKKKKKKKKKKKKKKKK:   Lambda= "<<currentLambda<<endln;
    theSOE->zeroB();//reset the SOE to zero ;Abbas
 
  } 
 
  return 0;
}
