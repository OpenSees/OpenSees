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
                                                                        
// $Revision: 1.3 $
// $Date: 2003-02-14 23:00:44 $
// $Source: /usr/local/cvs/OpenSees/SRC/analysis/analysis/DomainDecompositionAnalysis.cpp,v $
                                                                        
                                                                        
// File: ~/analysis/Analysis/DomainDecompositionAnalysis.C
// 
// Written: fmk 
// Created: Tue Sept 17 16:34:47: 1996
// Revision: A
//
// Description: This file contains the class definition for 
// DomainDecompositionAnalysis. DomainDecompositionAnalysis is a subclass 
// of AnalysisAnalysis, it is used to perform the static condensation process
// on a subdomain.
//
// What: "@(#) DomainDecompositionAnalysis.C, revA"


#include <DomainDecompositionAnalysis.h>
#include <ConstraintHandler.h>
#include <DOF_Numberer.h>
#include <AnalysisModel.h>
#include <LinearSOE.h>
#include <DomainDecompAlgo.h>
#include <DomainSolver.h>

#include <IncrementalIntegrator.h>
#include <Subdomain.h>

#include <FE_Element.h>
#include <DOF_Group.h>
#include <FE_EleIter.h>
#include <DOF_GrpIter.h>
#include <Matrix.h>
#include <ID.h>
#include <Node.h>

#include <Channel.h>
#include <FEM_ObjectBroker.h>

DomainDecompositionAnalysis::DomainDecompositionAnalysis(Subdomain &the_Domain)
:Analysis(the_Domain),
 MovableObject(DomDecompANALYSIS_TAGS_DomainDecompositionAnalysis),
 theSubdomain(&the_Domain),
 theHandler(0),
 theNumberer(0),
 theModel(0),
 theAlgorithm(0),
 theIntegrator(0),
 theSOE(0),
 theSolver(0),
 theResidual(0),numEqn(0),numExtEqn(0),tangFormed(false),tangFormedCount(0),
 domainStamp(0)
{
    theSubdomain->setDomainDecompAnalysis(*this);
}


DomainDecompositionAnalysis::DomainDecompositionAnalysis(int clsTag,
							 Subdomain &the_Domain)
:Analysis(the_Domain),
 MovableObject(clsTag),
 theSubdomain(&the_Domain),
 theHandler(0),
 theNumberer(0),
 theModel(0),
 theAlgorithm(0),
 theIntegrator(0),
 theSOE(0),
 theSolver(0),
 theResidual(0),numEqn(0),numExtEqn(0),tangFormed(false),tangFormedCount(0),
 domainStamp(0)
{
    theSubdomain->setDomainDecompAnalysis(*this);
}

DomainDecompositionAnalysis::DomainDecompositionAnalysis(Subdomain &the_Domain,
			   ConstraintHandler &handler,
			   DOF_Numberer &numberer,
			   AnalysisModel &model,
			   DomainDecompAlgo &theSolnAlgo,
			   IncrementalIntegrator &integrator,
			   LinearSOE &theLinSOE,
			   DomainSolver &theDDSolver)


:Analysis(the_Domain),
 MovableObject(DomDecompANALYSIS_TAGS_DomainDecompositionAnalysis),
 theSubdomain( &the_Domain),
 theHandler( &handler),
 theNumberer( &numberer),
 theModel( &model),
 theAlgorithm( &theSolnAlgo),
 theIntegrator( &integrator),
 theSOE( &theLinSOE),
 theSolver( &theDDSolver),
 theResidual(0),numEqn(0),numExtEqn(0),tangFormed(false),tangFormedCount(0)
{
    theModel->setLinks(the_Domain);
    theHandler->setLinks(*theSubdomain,*theModel,*theIntegrator);
    theNumberer->setLinks(*theModel);
    theIntegrator->setLinks(*theModel,*theSOE);
    theAlgorithm->setLinks(*theModel,*theIntegrator,*theSOE,
			   *theSolver,*theSubdomain);


    theSubdomain->setDomainDecompAnalysis(*this);
}    


DomainDecompositionAnalysis::~DomainDecompositionAnalysis()
{

}    

int 
DomainDecompositionAnalysis::analyze(void)
{
    opserr << "DomainDecompositionAnalysis::analyze(void)";
    opserr << "does nothing and should not have been called\n";
    return -1;
}

int
DomainDecompositionAnalysis::domainChanged(void)
{
    // remove existing FE_elements and DOF_Groups from the Analysis
    theModel->clearAll();
    theHandler->clearAll();

    // now we invoke handle() on the constraint handler which
    // causes the creation of FE_Element and DOF_Group objects
    // and their addition to the AnalysisModel.

    numExtEqn = theHandler->handle(&(theSubdomain->getExternalNodes()));

    // we now get a node to number last

    const ID &theExtNodes = theSubdomain->getExternalNodes();
    int idSize = theExtNodes.Size();
    //    int theLastDOF = -1;

    ID theLastDOFs(1);
    int cnt = 0;

    // create an ID containing the tags of the DOF_Groups that are to
    // be numbered last
    for (int i=0; i<idSize; i++) {
	int nodeTag = theExtNodes(i);
	Node *nodePtr = theSubdomain->getNode(nodeTag);
	DOF_Group *dofGrpPtr = nodePtr->getDOF_GroupPtr();
	if (dofGrpPtr != 0) {
	    const ID theID = dofGrpPtr->getID();
	    int size = theID.Size();
	    for (int j=0; j<size; j++)
		if (theID(j) == -3) {
		    theLastDOFs[cnt]  = dofGrpPtr->getTag();
		    cnt++;
		    j = size;
		}
	}
    }

    // we now invoke number() on the numberer which causes
    // equation numbers to be assigned to all the DOFs in the
    // AnalysisModel.    

    theNumberer->numberDOF(theLastDOFs);

    /*************************
    for (int i=0; i<idSize; i++) {
	int nodeTag = theExtNodes(i);
	Node *nodePtr = theSubdomain->getNode(nodeTag);
	DOF_Group *dofPtr = nodePtr->getDOF_GroupPtr();
	if (dofPtr != 0) {
	    const ID theID = dofPtr->getID();
	    int size = theID.Size();
	    for (int j=0; j<size; j++)
		if (theID(j) == -3) {
		    theLastDOF = dofPtr->getTag();
	            i = idSize;
                    j=size;
		}
	}
    }
    theNumberer->numberDOF(theLastDOF);
    **********************/


    // we invoke setSize() on the LinearSOE which
    // causes that object to determine its size    
    
    theSOE->setSize(theModel->getDOFGraph());    
    numEqn = theSOE->getNumEqn();

    // we invoke domainChange() on the integrator and algorithm

    theIntegrator->domainChanged();
    theAlgorithm->domainChanged();        

    // now set the variables to indicate that tangent has not been formed

    tangFormed = false;
    tangFormedCount = 0;
    
    return 0;
}


int
DomainDecompositionAnalysis::getNumExternalEqn(void)
{
    return numExtEqn;
}

int
DomainDecompositionAnalysis::getNumInternalEqn(void)
{
  return numEqn-numExtEqn;
}




int  
DomainDecompositionAnalysis::newStep(double dT)
{
  return theIntegrator->newStep(dT);
}



int  
DomainDecompositionAnalysis::computeInternalResponse(void)
{
  return theAlgorithm->solveCurrentStep();
}




int  
DomainDecompositionAnalysis::formTangent(void)
{
    int result =0;

    Domain *the_Domain = this->getDomainPtr();

    // we check to see if the domain has changed 
    int stamp = the_Domain->hasDomainChanged();
    if (stamp != domainStamp) {
	domainStamp = stamp;
	this->domainChanged();
    }
    
    // if tangFormed == -1 then formTangent has already been
    // called for this state by formResidual() or formTangVectProduct()
    // so we won't be doing it again.

    if (tangFormedCount != -1) {
	result = theIntegrator->formTangent();
	if (result < 0)
	    return result;
	result = theSolver->condenseA(numEqn-numExtEqn);
	if (result < 0)
	    return result;
    }
	
    tangFormed = true;
    tangFormedCount++;
    
    return result;
}



int  
DomainDecompositionAnalysis::formResidual(void)
{
    int result =0;
    Domain *the_Domain = this->getDomainPtr();    
    
    // we check to see if the domain has changed 
    int stamp = the_Domain->hasDomainChanged();
    if (stamp != domainStamp) {
	domainStamp = stamp;
	this->domainChanged();
    }
    
    if (tangFormed == false) {
	result = this->formTangent();
	if (result < 0)
	    return result;
	tangFormedCount = -1; // set to minus number so tangent 
	                      // is not formed twice at same state
    }

    result = theIntegrator->formUnbalance();

    if (result < 0)
	return result;
    return theSolver->condenseRHS(numEqn-numExtEqn);
}



int  
DomainDecompositionAnalysis::formTangVectProduct(Vector &u)
{
    int result = 0;

    Domain *the_Domain = this->getDomainPtr();
    
    // we check to see if the domain has changed 
    int stamp = the_Domain->hasDomainChanged();
    if (stamp != domainStamp) {
	domainStamp = stamp;
	this->domainChanged();
    }
    
    if (tangFormed == false) {
	result = this->formTangent();
	if (result < 0)
	    return result;
	tangFormedCount = -1; // set to minus number so tangent 
	                      // is not formed twice at same state
    }    

    return theSolver->computeCondensedMatVect(numEqn-numExtEqn,u);
}



const Matrix &
DomainDecompositionAnalysis::getTangent()
{
    Domain *the_Domain = this->getDomainPtr();
    
    // we check to see if the domain has changed 
    int stamp = the_Domain->hasDomainChanged();
    if (stamp != domainStamp) {
	domainStamp = stamp;
	this->domainChanged();
    }

    if (tangFormed == false) {
	this->formTangent();
    }
    
    return (theSolver->getCondensedA());
}



const Vector &
DomainDecompositionAnalysis::getResidual()
{

    Domain *the_Domain = this->getDomainPtr();
    
    // we check to see if the domain has changed 
    // we check to see if the domain has changed 
    int stamp = the_Domain->hasDomainChanged();
    if (stamp != domainStamp) {
	domainStamp = stamp;
	this->domainChanged();
	this->formResidual();	
    }
    
    if (theResidual == 0) {
	theResidual = new Vector(theSolver->getCondensedRHS());
	return *theResidual;	
    }
    else if (theResidual->Size() != numExtEqn) {
	delete theResidual;
	theResidual = new Vector(theSolver->getCondensedRHS());
	return *theResidual;		    
    }
    else {
	(*theResidual) = theSolver->getCondensedRHS();
    }

    return *theResidual;
}




const Vector &
DomainDecompositionAnalysis::getTangVectProduct()
{
    Domain *the_Domain = this->getDomainPtr();
    
    // we check to see if the domain has changed 
    int stamp = the_Domain->hasDomainChanged();
    if (stamp != domainStamp) {
	domainStamp = stamp;
	this->domainChanged();
    }
    
    return theSolver->getCondensedMatVect();
}





Subdomain  *
DomainDecompositionAnalysis::getSubdomainPtr(void) const
{
    return theSubdomain;
}




ConstraintHandler *
DomainDecompositionAnalysis::getConstraintHandlerPtr(void) const
{
    return theHandler;
}



DOF_Numberer *
DomainDecompositionAnalysis::getDOF_NumbererPtr(void) const
{
    return theNumberer;
}



AnalysisModel  *
DomainDecompositionAnalysis::getAnalysisModelPtr(void) const
{
    return theModel;
}



DomainDecompAlgo  *
DomainDecompositionAnalysis::getDomainDecompAlgoPtr(void) const
{
    return theAlgorithm;
}



IncrementalIntegrator *
DomainDecompositionAnalysis::getIncrementalIntegratorPtr(void) const
{
    return theIntegrator;
}



LinearSOE *
DomainDecompositionAnalysis::getLinSOEPtr(void) const
{
    return theSOE;    
}



DomainSolver *
DomainDecompositionAnalysis::getDomainSolverPtr(void) const
{
    return theSolver;
}


int 
DomainDecompositionAnalysis::sendSelf(int commitTag,
				      Channel &theChannel)
{
    // determine the type of each object in the aggregation,
    // store it in an ID and send the info off.
    int dataTag = this->getDbTag();
    ID data(14);
    data(0) = theHandler->getClassTag();
    data(1) = theNumberer->getClassTag();
    data(2) = theModel->getClassTag();
    data(3) = theAlgorithm->getClassTag();
    data(4) = theIntegrator->getClassTag();    
    data(5) = theSOE->getClassTag();    
    data(6) = theSolver->getClassTag();        

    data(7) = theHandler->getDbTag();
    data(8) = theNumberer->getDbTag();
    data(9) = theModel->getDbTag();
    data(10) = theAlgorithm->getDbTag();
    data(11) = theIntegrator->getDbTag();    
    data(12) = theSOE->getDbTag();    
    data(13) = theSolver->getDbTag();        

    theChannel.sendID(dataTag, commitTag, data);

    // invoke sendSelf on each object in the aggregation
    
    theHandler->sendSelf(commitTag, theChannel);
    theNumberer->sendSelf(commitTag, theChannel);
    theModel->sendSelf(commitTag, theChannel);
    theAlgorithm->sendSelf(commitTag, theChannel);
    theIntegrator->sendSelf(commitTag, theChannel);    
    theSOE->sendSelf(commitTag, theChannel);    
    theSolver->sendSelf(commitTag, theChannel);            
    return 0;
}

int 
DomainDecompositionAnalysis::recvSelf(int commitTag, 
				      Channel &theChannel, 
				      FEM_ObjectBroker &theBroker)
{
    // receive the data identifyng the objects in the aggregation
    ID data(14);
    int dataTag = this->getDbTag();
    theChannel.recvID(dataTag, commitTag, data);

    //
    // now ask the object broker an object of each type
    // and invoke recvSelf() on the object to init it.
    //

    theHandler = theBroker.getNewConstraintHandler(data(0));
    if (theHandler != 0) {
	theHandler->setDbTag(data(7));
	theHandler->recvSelf(commitTag, theChannel,theBroker);
    }
    else {
	opserr << "DomainDecompositionAnalysis::recvSelf";
	opserr << " - failed to get the ConstraintHandler\n";
	return -1;
    }



    theNumberer = theBroker.getNewNumberer(data(1));
    if (theNumberer != 0) {
	theNumberer->setDbTag(data(8));	
	theNumberer->recvSelf(commitTag, theChannel,theBroker);
    }
    else {
	opserr << "DomainDecompositionAnalysis::recvSelf";
	opserr << " - failed to get the DOF Numberer\n";
	return -1;
    }    


    theModel = theBroker.getNewAnalysisModel(data(2));
    if (theModel != 0) {
	theModel->setDbTag(data(9));
	theModel->recvSelf(commitTag, theChannel,theBroker);
    }
    else {
	opserr << "DomainDecompositionAnalysis::recvSelf";
	opserr << " - failed to get the AnalysisModel\n";
	return -1;
    }        




    theAlgorithm = theBroker.getNewDomainDecompAlgo(data(3));
    if (theAlgorithm != 0) {
	theAlgorithm->setDbTag(data(10));
	theAlgorithm->recvSelf(commitTag, theChannel,theBroker);
    }
    else {
	opserr << "DomainDecompositionAnalysis::recvSelf";
	opserr << " - failed to get the Domain Decomp Algo\n";
	return -1;
    }            

    theIntegrator = theBroker.getNewIncrementalIntegrator(data(4));
    if (theIntegrator != 0) {
	theIntegrator->setDbTag(data(11));
	theIntegrator->recvSelf(commitTag, theChannel,theBroker);
    }
    else {
	opserr << "DomainDecompositionAnalysis::recvSelf";
	opserr << " - failed to get the IncrementalIntegrator\n";
	return -1;
    }        	

    theSOE = theBroker.getPtrNewDDLinearSOE(data(5),data(6));
    theSolver = theBroker.getNewDomainSolver();

    if (theSOE == 0 || theSolver == 0) {
	opserr << "DomainDecompositionAnalysis::recvSelf";
	opserr << " - failed to get the LinearSOE and the DomainSolver \n";
	return -1;
    }  else {
	theSOE->setDbTag(data(12));
	theSolver->setDbTag(data(13));
	theSOE->recvSelf(commitTag, theChannel,theBroker);
	theSolver->recvSelf(commitTag, theChannel,theBroker);
    }

    // set the links in all the objects

    theModel->setLinks(*theSubdomain);
    theHandler->setLinks(*theSubdomain,*theModel,*theIntegrator);
    theNumberer->setLinks(*theModel);
    theIntegrator->setLinks(*theModel,*theSOE);
    theAlgorithm->setLinks(*theModel,*theIntegrator,*theSOE,
			   *theSolver,*theSubdomain);

    theSubdomain->setDomainDecompAnalysis(*this);

    return 0;
}




