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

// $Revision: 1.13 $
// $Date: 2010-06-01 23:47:54 $
// $Source: /usr/local/cvs/OpenSees/SRC/analysis/integrator/DisplacementControl.cpp,v $


// File: ~/analysis/integrator/DisplacementControl.C
// 
// Written: fmk 
// Created: 07/98
// Revision: A
//
// Description: This file contains the class definition for DisplacementControl.
// DisplacementControl is an algorithmic class for performing a static analysis
// using the arc length scheme, that is within a load step the following
// constraint is enforced: dU^TdU + alpha^2*dLambda^2 = DisplacementControl^2
// where dU is change in nodal displacements for step, dLambda is
// change in applied load and DisplacementControl is a control parameter.
//
// What: "@(#) DisplacementControl.C, revA"


#include <DisplacementControl.h>
#include <AnalysisModel.h>
#include <LinearSOE.h>
#include <Vector.h>
#include <Channel.h>
#include <math.h>
#include <Domain.h>
#include <Node.h>
#include <DOF_Group.h>
#include <DOF_GrpIter.h>
#include <ID.h>
#include <stdlib.h>
#include <FE_Element.h>
#include <FE_EleIter.h>
#include <LoadPattern.h>
#include <LoadPatternIter.h>
#include <Parameter.h>
#include <ParameterIter.h>
#include <EquiSolnAlgo.h>
#include <Matrix.h>
#include <TaggedObjectStorage.h>
#include <elementAPI.h>

void* OPS_DisplacementControlIntegrator()
{
    if(OPS_GetNumRemainingInputArgs() < 3) {
       opserr<<"insufficient arguments for DisplacementControl\n";
       return 0;
    }

    // node, dof
    int iData[2];
    int numData = 2;
    if(OPS_GetIntInput(&numData,&iData[0]) < 0) {
	opserr << "WARNING failed to read node tag and ndf\n";
	return 0;
    }

    double incr;
    numData = 1;
    if(OPS_GetDoubleInput(&numData,&incr) < 0) {
	opserr << "WARNING failed to read incr\n";
	return 0;
    }

    // numIter,dumin,dumax
    int numIter = 1;
    int formTangent = 0;
    double data[2] = {incr,incr};
    if(OPS_GetNumRemainingInputArgs() > 2) {
       numData = 1;
       if(OPS_GetIntInput(&numData,&numIter) < 0) {
	   opserr << "WARNING failed to read numIter\n";
	   return 0;
       }
       numData = 2;
       if(OPS_GetDoubleInput(&numData,&data[0]) < 0) {
	   opserr << "WARNING failed to read dumin and dumax\n";
	   return 0;
       }
    }

    if (OPS_GetNumRemainingInputArgs() == 1) {
      	std::string type = OPS_GetString();
	if(type=="-initial" || type=="-Initial") {
	    formTangent = 1;
	} 
    }

    // check node
    Domain* theDomain = OPS_GetDomain();
    Node *theNode = theDomain->getNode(iData[0]);
    if(theNode == 0) {
       opserr << "WARNING integrator DisplacementControl node dof dU : Node does not exist\n";
       return 0;
    }

    int numDOF = theNode->getNumberDOF();
    if(iData[1] <= 0 || iData[1] > numDOF) {
       opserr << "WARNING integrator DisplacementControl node dof dU : invalid dof given\n";
       return 0;
    }

    return new DisplacementControl(iData[0],iData[1]-1,
				   incr,theDomain,
				   numIter,data[0],data[1], 
				   formTangent);
}


DisplacementControl::DisplacementControl(int node, int dof, 
					 double increment, 
					 Domain *domain,
					 int numIncr,
					 double min, double max, int tang) 
:StaticIntegrator(INTEGRATOR_TAGS_DisplacementControl),
   theNode(node), theDof(dof), theIncrement(increment), theDomain(domain),
   theDofID(-1),
   deltaUhat(0), deltaUbar(0), deltaU(0), deltaUstep(0),dUhatdh(0),
   phat(0), deltaLambdaStep(0.0), currentLambda(0.0), dLambdaStepDh(0.0),dUIJdh(0),Dlambdadh(0.0),
   specNumIncrStep(numIncr), numIncrLastStep(numIncr),
   minIncrement(min), maxIncrement(max),sensitivityFlag(0),Residual(0),dlambdadh(0.0),
   dLambda(0.0),  sensU(0),d_deltaU_dh(0),Residual2(0),gradNumber(0),dLAMBDAdh(0),dphatdh(0)
{
  tangFlag = tang;

   // to avoid divide-by-zero error on first update() ensure numIncr != 0
   if (numIncr == 0) {
      opserr << "WARNING DisplacementControl::DisplacementControl() -";
      opserr << " numIncr set to 0, 1 assumed\n";
      specNumIncrStep = 1.0;
      numIncrLastStep = 1.0;
           
   }

}

DisplacementControl::~DisplacementControl()
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
DisplacementControl::newStep(void)
{

  if (theDofID == -1) {
    opserr << "DisplacementControl::newStep() - dof is fixed or constrained (or domainChanged has not been called!)\n";
    return -1;
   }

   // get pointers to AnalysisModel and LinearSOE
   AnalysisModel *theModel = this->getAnalysisModel();
   LinearSOE *theLinSOE = this->getLinearSOE();    
   if (theModel == 0 || theLinSOE == 0) {
      opserr << "WARNING DisplacementControl::newStep() ";
      opserr << "No AnalysisModel or LinearSOE has been set\n";
      return -1;
   }

   // determine increment for this iteration
   double factor = specNumIncrStep/numIncrLastStep;
   theIncrement *=factor;

   if (theIncrement < minIncrement)
      theIncrement = minIncrement;
   else if (theIncrement > maxIncrement)
      theIncrement = maxIncrement;


   // get the current load factor
   currentLambda = theModel->getCurrentDomainTime();

   // determine dUhat
   this->formTangent(tangFlag);
   theLinSOE->setB(*phat);
   if (theLinSOE->solve() < 0) {
      opserr << "DisplacementControl::newStep(void) - failed in solver\n";
      return -1;
   }

   (*deltaUhat) = theLinSOE->getX();
   Vector &dUhat = *deltaUhat;// this is the Uft in the nonlinear lecture notes
   double dUahat = dUhat(theDofID);// this is the component of the Uft in our nonlinear lecture notes

   if (dUahat == 0.0) {
      opserr << "WARNING DisplacementControl::newStep() ";
      opserr << "dUahat is zero -- zero reference displacement at control node DOF\n";
      return -1;
   }

   // determine delta lambda(1) == dlambda    
   double dlambda = theIncrement/dUahat; // this is the dlambda of the 1st step

  // calldLambda1dh=theIncrement;
   deltaLambdaStep = dlambda;
   currentLambda += dlambda;

   // determine delta U(1) == dU
   (*deltaU) = dUhat;
   (*deltaU) *= dlambda;// this is eq(4) in the paper {dU}_1=dLAmbda1*Uft.
   (*deltaUstep) = (*deltaU);


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

      this->formdLambdaDh(gradNumber);
      theParam->activate(false);
    } 
  }
  ///////////////Abbas/////////////////////////////

   // update model with delta lambda and delta U
   theModel->incrDisp(*deltaU); 
   theModel->applyLoadDomain(currentLambda);
   if (theModel->updateDomain() < 0) {
      opserr << "DisplacementControl::newStep - model failed to update for new dU\n";
      return -1;
   }

   numIncrLastStep = 0;
   return 0;
}

int DisplacementControl::update(const Vector &dU)
{
 //  opserr<<"Update: Start"<<endln;

   if (theDofID == -1) {
      opserr << "DisplacementControl::newStep() - domainChanged has not been called\n";
      return -1;
   } 
   AnalysisModel *theModel = this->getAnalysisModel();
   LinearSOE *theLinSOE = this->getLinearSOE();    
   if (theModel == 0 || theLinSOE == 0) {
      opserr << "WARNING DisplacementControl::update() ";
      opserr << "No AnalysisModel or LinearSOE has been set\n";

      return -1;
   }

   (*deltaUbar) = dU; // have to do this as the SOE is gonna change
   double dUabar = (*deltaUbar)(theDofID);//dUbar is the vector of residual displacement and dUabar is its component
    // opserr<<" DisplacementControl:: deltaUbar = "<<*deltaUbar<<endln; 

   // determine dUhat    
   theLinSOE->setB(*phat);
   theLinSOE->solve();
   (*deltaUhat) = theLinSOE->getX();    

   double dUahat = (*deltaUhat)(theDofID);
   if (dUahat == 0.0) {
      opserr << "WARNING DisplacementControl::update() ";
      opserr << "dUahat is zero -- zero reference displacement at control node DOF\n";
      return -1;
   }

   // determine delta lambda(1) == dlambda    
   dLambda = -dUabar/dUahat;// this dLambda i,j

   // determine delta U(i)
   (*deltaU) = (*deltaUbar);   
   deltaU->addVector(1.0, *deltaUhat,dLambda);


   // update dU and dlambda
   (*deltaUstep) += *deltaU;
   deltaLambdaStep += dLambda;
   currentLambda += dLambda;

   // update the model
   theModel->incrDisp(*deltaU);    
   theModel->applyLoadDomain(currentLambda);    
   if (theModel->updateDomain() < 0) {
      opserr << "DisplacementControl::update - model failed to update for new dU\n";
      return -1;
   }

   // set the X soln in linearSOE to be deltaU for convergence Test
   theLinSOE->setX(*deltaU);

   numIncrLastStep++;

   return 0;
}



int 
DisplacementControl::domainChanged(void)
{
   // we first create the Vectors needed
  //opserr<<" this is the domain change function"<<endln;//Abbas
   AnalysisModel *theModel = this->getAnalysisModel();
   LinearSOE *theLinSOE = this->getLinearSOE(); 
    if (theModel == 0 || theLinSOE == 0) {
      opserr << "WARNING DisplacementControl::update() ";
      opserr << "No AnalysisModel or LinearSOE has been set\n";
      return -1;
   }    
   int size = theModel->getNumEqn(); // ask model in case N+1 space

   if (deltaUhat == 0 || deltaUhat->Size() != size) { // create new Vector
      if (deltaUhat != 0)
	 delete deltaUhat;   // delete the old
      deltaUhat = new Vector(size);
      if (deltaUhat == 0 || deltaUhat->Size() != size) { // check got it
	 opserr << "FATAL DisplacementControl::domainChanged() - ran out of memory for";
	 opserr << " deltaUhat Vector of size " << size << endln;
	 exit(-1);
      }
   }
  

   if (deltaUbar == 0 || deltaUbar->Size() != size) { // create new Vector
      if (deltaUbar != 0)
	 delete deltaUbar;   // delete the old
      deltaUbar = new Vector(size);
      if (deltaUbar == 0 || deltaUbar->Size() != size) { // check got it
	opserr << "FATAL DisplacementControl::domainChanged() - ran out of memory for";
	opserr << " deltaUbar Vector of size " << size << endln;
	exit(-1);
      }
   }

   if (deltaU == 0 || deltaU->Size() != size) { // create new Vector
     if (deltaU != 0)
       delete deltaU;   // delete the old
     deltaU = new Vector(size);
     if (deltaU == 0 || deltaU->Size() != size) { // check got it
       opserr << "FATAL DisplacementControl::domainChanged() - ran out of memory for";
       opserr << " deltaU Vector of size " << size << endln;
       exit(-1);
      }
   }

   if (deltaUstep == 0 || deltaUstep->Size() != size) { 
      if (deltaUstep != 0)
	delete deltaUstep;  
      deltaUstep = new Vector(size);
      if (deltaUstep == 0 || deltaUstep->Size() != size) { 
	 opserr << "FATAL DisplacementControl::domainChanged() - ran out of memory for";
	 opserr << " deltaUstep Vector of size " << size << endln;
	 exit(-1);
      }
   }

   if (phat == 0 || phat->Size() != size) { 
      if (phat != 0)
	 delete phat;  
      phat = new Vector(size);
      if (phat == 0 || phat->Size() != size) { 
	 opserr << "FATAL DisplacementControl::domainChanged() - ran out of memory for";
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
   
   if (Residual2 == 0 || Residual2->Size() != size) { 
     if (Residual2 != 0)
       delete Residual2;  
     Residual2 = new Vector(size);
     if (Residual2== 0 || Residual2->Size() != size) { 
       opserr << "FATAL DisplacementControl::domainChanged() - ran out of memory for";
       opserr << " N Vector of size " << size << endln;
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

   Domain *theDomain=theModel->getDomainPtr();//Abbas
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
      opserr << "WARNING DisplacementControl::domainChanged() - zero reference load";
      return -1;
   }

   // lastly we determine the id of the nodal dof
   // EXTRA CODE TO DO SOME ERROR CHECKING REQUIRED

   Node *theNodePtr = theDomain->getNode(theNode);
   if (theNodePtr == 0) {
      opserr << "DisplacementControl::domainChanged - no node\n";
      return -1;
   }

   DOF_Group *theGroup = theNodePtr->getDOF_GroupPtr();
   if (theGroup == 0) {
      return 0;
   }
   const ID &theID = theGroup->getID();
   theDofID = theID(theDof);
   return 0;
}

int
DisplacementControl::sendSelf(int cTag,
      Channel &theChannel)
{
   // TO FINISH
   return 0;
}


int
DisplacementControl::recvSelf(int cTag,
      Channel &theChannel, FEM_ObjectBroker &theBroker)
{
   // TO FINISH
   return 0;
}

void
DisplacementControl::Print(OPS_Stream &s, int flag)
{
   // TO FINISH    
}




///////////////////////////Sensitivity Begin///////////////////////////////////
//Added by Abbas
//obtain the derivative of the tangent displacement (dUhatdh)
   Vector *
DisplacementControl::formTangDispSensitivity(Vector *dUhatdh,int gradNumber)
{
  // opserr<<"DisplacementControl:; formTangDispsensitivity:::::::::::Start"<<endln;
   LinearSOE *theLinSOE = this->getLinearSOE(); 
   dUhatdh->Zero();
   dphatdh->Zero();
//   opserr<<"DisplacementControl::dUhatdh is  "<<*dUhatdh<<endln;//

   // To get the structural stiffness Matrix using the full general system of equations
   //...............................................................
 //  this->formTangent();
 //  Matrix K(size,size);
 //   K.Zero();
  // K=this->ActivateSensitivity();
  // // to print K to the screen
   //..........................
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

   // form dKdh

 // Matrix dKdh(size,size);
 //  dKdh.Zero();
//  dKdh=this->getdKdh();
//
   // To print dKdh to the Screen 
   //................................................................
   
//  opserr<<"dKdh from the DisplacementControl.cpp"<<endln;
//      for(int i=0;i<size;i++)
 //     {
 //     for(int j=0;j<size;j++)
 //     {
 //    opserr<<dKdh(i,j)<<"      ";

  //    }
 //    opserr<<endln;
  //    }
      
   //...............................................................


  // form dKdh*deltaUhat
  // dUhatdh->addMatrixVector(0.0,dKdh,*deltaUhat,-1.0);
   
   
   //call the tangent (K)
   this->formTangent(tangFlag);
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

    (*dUhatdh)=theLinSOE->getX();


/////////////////////////////////////////////////////////
   return dUhatdh;
}

// form dLambda for each time step dLambda
double 
DisplacementControl::formdLambdaDh(int gradNumber)
{
  Vector &duHatdh=*dUhatdh;

  //the component of the duHatdh vector
  double duHatdh_Comp=duHatdh(theDofID);

  // call the deltaUhat component
  Vector &UFT=(*deltaUhat);

  double UFT_Comp=UFT(theDofID);// cll the component of the dUhat again
  if(UFT_Comp == 0.0)
    dlambdadh=0.0;// to avoid dividing by zero
  else
    dlambdadh=-(theIncrement*duHatdh_Comp)/(UFT_Comp*UFT_Comp);

  if(dLAMBDAdh != 0) {
    (*dLAMBDAdh)(gradNumber) = (*dLAMBDAdh)(gradNumber) + dlambdadh;
    return (*dLAMBDAdh)(gradNumber);
  } else {
    return 0.0;
  }
}

/*
 // The following function is needed if dKdh is considered in the derivation (dKdh #0)
   Vector *
DisplacementControl::formResidualDispSensitivity( Vector *dUIJdh,int gradNumber)
{
   // this is just the product of the dKdh*deltaUbar. The complete solution to the residual displacement  
   // is inside the computeSensitivities function.
 
 
   Matrix dKdh(size,size);
   dKdh.Zero();
   dKdh=this->getdKdh();
   dUIJdh-> addMatrixVector(0.0,dKdh,*deltaUbar,-1.0);
   opserr<<" the dKdh*deltaUbar=  "<<*dUIJdh<<endln;

   return dUIJdh;

}

*/

 // dLambdadh of the subsequent iterations (J>1)
double 
DisplacementControl::getLambdaSensitivity(int gradNumber)
{

  //call the tangent displacement component (deltaUhat)
  // LinearSOE *theLinSOE = this->getLinearSOE(); 

  //  this->formTangent();
  //  theLinSOE->setB(*phat);
  //   if (theLinSOE->solve() < 0) {
  //      opserr << "DisplacementControl::newStep(void) - failed in solver\n";
  //      return -1;
  //   }

  //  (*deltaUhat) = theLinSOE->getX();
 
   Vector &UFT=*deltaUhat;
   double UFT_Comp=UFT(theDofID);

   //call the derivative of the tangent displacement component (dUftdh).
   Vector &duHatdh=*dUhatdh;
   double duHatdh_Comp=duHatdh(theDofID);

   //Call the residual displacement component
   Vector &dufR=*deltaUbar;
   double dufR_Comp=dufR(theDofID);

   Vector &dufRdh=*dUIJdh;// component of the dUfrDh: derivative of the residual displacement
   double dufRdh_Comp=dufRdh(theDofID);

   if(UFT_Comp==0.0 )
     Dlambdadh=0.0;
   else
     // dLambdadh_ij( the sensitivity of the load component to the parameter h)
     Dlambdadh=((-dufRdh_Comp*UFT_Comp  +(dufR_Comp*duHatdh_Comp)))/(UFT_Comp*UFT_Comp);  //   
 
   // Now update Lambda_ij
   if(dLAMBDAdh !=0) {
     (*dLAMBDAdh)(gradNumber) = (*dLAMBDAdh)(gradNumber)+ Dlambdadh;
     return (*dLAMBDAdh)(gradNumber);
   } else {
     return 0.0;
   }
}



int
DisplacementControl::formEleResidual(FE_Element* theEle)
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
DisplacementControl::formIndependentSensitivityRHS()
{
   return 0;
}


   
int
DisplacementControl::formSensitivityRHS(int passedGradNumber)
{
  sensitivityFlag = 1;
  //this->Activatesensitivity();//Abbas
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

  //  opserr<<" dKdh*deltaUbar after while loop is "<<*dUIJdh<<endln;
  (*Residual)=theSOE->getB();
 
  // opserr<<"RHS():getB()............"<<theSOE->getB()<<endln;
  // (*Residual) += (*dUIJdh);

  //   (*Residual2)=theSOE->getB();

  //  opserr<<" the new RHS is  "<<*Residual<<endln;
  //  opserr<<"where dUIJ is "<<*dUIJdh<<endln;

  
  //    if(CallParam==1) {
  //  opserr<<"inside if statement of gradIndex # "<<gradNumber<<endln;
  //  this->formTangDispSensitivity(dUhatdh,gradNumber);
  //  this->formdLambdaDh(gradNumber);
  // }
  

  //   opserr<<"RHS: dLambdaStepdh is "<<(*dLAMBDAdh)(gradNumber)<<endln;

  double CallDlambda1dh=(*dLAMBDAdh)(gradNumber);
  // opserr<<"RHS:: dphatdh is ....................................//"<<*dphatdh<<endln;
  Residual->addVector(1.0,*phat, CallDlambda1dh ); //needed to calculate dLambdadh
  //opserr<<"The residual sensitivity... dRdh= "<<*Residual<<endln;
  //    opserr<<"Pref=  "<<*phat<<endln;
  //    opserr<<"dLambdadh= "<<CallDlambda1dh<<endln;
  //opserr<<" dPrdh after adding dLamdh= "<<*Residual<<endln;
  Residual->addVector(1.0,*dphatdh, currentLambda ); //needed to calculate dLambdadh
  
  //  Matrix dKdh(size,size);
  //	 dKdh.Zero();
  //	 dKdh=this->getdKdh();
  //	 Residual->addMatrixVector(1.0,dKdh,*deltaUbar,-1.0);
  
  
  /*
    Matrix Kt(size,size);
    Kt.Zero();
    Kt=this->ActivateSensitivity();
    if(sensU==0) {
    opserr<<"sensU is zero ........."<<endln;
    sensU = new Vector(size);
    
    sensU->Zero();
    }
    Residual->addMatrixVector(1.0,Kt,(*sensU),-1.0);
    
    opserr<<"Resid= "<<*Residual<<endln;
  */
  




  Residual2->addVector(1.0,*phat,CallDlambda1dh);// needed to calculate dUdh
  //  Residual2->addVector(1.0,*dphatdh, currentLambda );
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
      //	 opserr<<"No sensitivity Load Parameter is involved"<<endln;
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
  
  //  (*Residual2)=theSOE->getB();
  
  
  //reset sensitivity flag
  sensitivityFlag=0;
  //  opserr<<"RHS :Ends"<<endln;
  return 0;
}


int
DisplacementControl::saveSensitivity(const Vector &v, int gradNum, int numGrads)
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
DisplacementControl::saveLambdaSensitivity(double dlambdadh, int gradNum, int numGrads)
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
DisplacementControl::commitSensitivity(int gradNum, int numGrads)
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
DisplacementControl::computeSensitivityAtEachIteration()
{

   return true;     
}


   int 
DisplacementControl::computeSensitivities(void)
{
  LinearSOE *theSOE = this->getLinearSOE();
  
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

 
  while ((theParam = paramIter()) != 0) {
    // Activate this parameter
    theParam->activate(true);
    
    // Zero the RHS vector
    theSOE->zeroB();
    
    // Get the grad index for this parameter
    int gradIndex = theParam->getGradIndex();
                  
    // this->formResidualDispSensitivity( dUIJdh, gradIndex);//./.
 
    // Form the RHS
    
    //  this->formTangDispSensitivity(dUhatdh,gradIndex);
    this->formSensitivityRHS(gradIndex);
     
    this->formTangent(tangFlag);
    theSOE->solve();
    *dUIJdh=theSOE->getX();// sensitivity of the residual displacement
 
    //  opserr<<"computeSensitivities: dUfRdh is  "<<*dUIJdh<<endln;//
    //    opserr<<"deltaUbar is "<<*deltaUbar<<endln;   
    this->formTangDispSensitivity(dUhatdh,gradIndex);
    double dlamdh = this->getLambdaSensitivity(gradIndex);

    // To obtain the response sensitivity 
    theSOE->setB(*Residual);
    theSOE->solve();
    //theSOE->getX();
       
    //Vector *x=new Vector(size);
    //(*x)=theSOE->getX();
    //  x->addVector(1.0,*deltaUhat,Dlambdadh);
    //  x->addVector(1.0,*dUhatdh,dLambda);
    (*sensU) = theSOE->getX();//.................................
    //   sensU->addVector(1.0,*deltaUhat,Dlambdadh);
    //    opserr<<" computeSensitivities::...... final dudh is "<<theSOE->getX()<<endln;
    //(*sensU) +=(*x);

    // Save sensitivity to nodes
    this->saveSensitivity( (*sensU), gradIndex, numGrads );
    this->saveLambdaSensitivity(dlamdh, gradIndex, numGrads);
    
    // Commit unconditional history variables (also for elastic problems; strain sens may be needed anyway)
    this->commitSensitivity(gradIndex, numGrads);
    // De-activate this parameter for next sensitivity calc
    theParam->activate(false);

    theSOE->zeroB();//reset the SOE to zero ;Abbas
 
  } 
  // end of if statement to be run only one time during the iteration process.
  //  opserr<<"computeSensitivity: End"<<endln;
  //  CallParam=0;
  // opserr<<"CallParam now is "<<CallParam<<endln;
  //opserr<<" computeSensitivities: end"<<endln;

  return 0;
}
