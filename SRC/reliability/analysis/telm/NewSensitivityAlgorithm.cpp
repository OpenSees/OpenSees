/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 2001, The Regents of the University of California    **
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
** Reliability module developed by:                                   **
**   Terje Haukaas (haukaas@ce.berkeley.edu)                          **
**   Armen Der Kiureghian (adk@ce.berkeley.edu)                       **
**                                                                    **
** ****************************************************************** */

// $Revision: 1.3 $
// $Date: 2008-08-26 17:34:50 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/analysis/telm/NewSensitivityAlgorithm.cpp,v $
                                                                        
#include <NewSensitivityAlgorithm.h>
#include <LinearSOE.h>
#include <Vector.h>
#include <EquiSolnAlgo.h>
#include <ReliabilityDomain.h>
//#include <RandomVariablePositioner.h>
//#include <RandomVariablePositionerIter.h>
//#include <ParameterPositionerIter.h>

#include <fstream>
#include <iomanip>
#include <iostream>

using std::ifstream;
using std::ios;

using std::setw;
using std::setprecision;


NewSensitivityAlgorithm::NewSensitivityAlgorithm(ReliabilityDomain *passedReliabilityDomain,
						 Domain* passedFEDomain,
						 EquiSolnAlgo *passedAlgorithm,
						 Integrator *passedSensitivityIntegrator,					
						 int passedAnalysisTypeTag)
  
:SensitivityAlgorithm(passedFEDomain, passedAlgorithm,
										   passedSensitivityIntegrator, passedAnalysisTypeTag)
{
	// The reliability domain is needed to get hold 
	// of the random variable positioners:
	theReliabilityDomain = passedReliabilityDomain;

	// The finite element equation solution algorithm is 
	// needed to get hold of the system of equations (SOE):
	theAlgorithm = passedAlgorithm;

	// The sensitivity integrator is needed to assemble the 
	// new right-hand side of the system of equations:
	theSensitivityIntegrator = passedSensitivityIntegrator;

	// Tag to tell whether grads should be computed at each step
	// and whether they should be computed wrt. random variables
	analysisTypeTag = passedAnalysisTypeTag;

	//gradRVPositioner=0;
	//gradParaPositioner=0;
	idGradPositioner=0;
	numGradPositioner=0;

	output.open("sensAlgo.txt", ios::out);

}




NewSensitivityAlgorithm::~NewSensitivityAlgorithm()
{
}



int 
NewSensitivityAlgorithm::computeSensitivities(bool fromFEM)
{
	// Meaning of analysisTypeTag:
	// 1: compute at each step wrt. random variables
	// 2: compute at each step wrt. parameters
	// 3: compute by command wrt. random variables
	// 4: compute by command wrt. parameters

	
	// Initial declarations
	int gradNumber;

//	int numGrads, numPos, i, gradNumber, posNumber;
//	RandomVariablePositioner *theRandomVariablePositioner;
//	ParameterPositioner *theParameterPositioner;




	// Get pointer to the system of equations (SOE)
//	LinearSOE *theSOE = theAlgorithm->getLinearSOEptr();


	// Get pointer to incremental integrator
//	IncrementalIntegrator *theIncInt = theAlgorithm->getIncrementalIntegratorPtr();


	// Form current tangent at converged state
	// (would be nice with an if-statement here in case
	// the current tangent is already formed)
	if (!fromFEM){
		if (theIncInt->formTangent(CURRENT_TANGENT) < 0){
			opserr << "WARNING SensitivityAlgorithm::computeGradients() -";
			opserr << "the Integrator failed in formTangent()\n";
			return -1;
		}
	}
//	if (theIncInt->formTangent(CURRENT_TANGENT) < 0){
//		opserr << "WARNING SensitivityAlgorithm::computeGradients() -";
//		opserr << "the Integrator failed in formTangent()\n";
//		return -1;
//	}


	// Get number of random variables and random variable positioners
//	if (analysisTypeTag==1 || analysisTypeTag==3) {
//		numGrads = theReliabilityDomain->getNumberOfRandomVariables();
//		numPos = theReliabilityDomain->getNumberOfRandomVariablePositioners();
//	}
//	else {
//		numPos = theReliabilityDomain->getNumberOfParameterPositioners();
//		numGrads = numPos;
//	}
	

	// Zero out the old right-hand side of the SOE
	theSOE->zeroB();
		

	// Form the part of the RHS which are indepent of parameter
	theSensitivityIntegrator->formIndependentSensitivityRHS();

	if(analysisTypeTag==1 || analysisTypeTag==3){
		int act;
		bool done=false;
		for (gradNumber=1; gradNumber<=numGrads; gradNumber++ )  {
			int idtemp=gradNumber-1;
			bool activeExist=false;
			for(int i=0;i<numGradPositioner[idtemp];i++){
				int ii=idGradPositioner[idtemp][i];
				//act=gradRVPositioner[ii]->activate(true);
				if(act!=0)activeExist=true; 
			}
			if(activeExist){
				theSOE->zeroB();
				// Form new right-hand side
				theSensitivityIntegrator->formSensitivityRHS(gradNumber);
				// Solve the system of equation with the new right-hand side
	
                /*
				Vector outv=theSOE->getA();
				double a1=outv(0);
				double a2=outv(1);
				double a3=outv(2);
				double a4=outv(3);
				double a5=outv(4);
				double a6=outv(5);
				double a7=outv(6);
				double a8=outv(7);
				double a9=outv(8);
				double a10=outv(9);
				double a11=outv(10);
				double a12=outv(11);
				double a13=outv(12);
				double a14=outv(13);
				double a15=outv(14);
				double a16=outv(15);
				double a17=outv(16);
				double a18=outv(17);
				double a19=outv(18);
				double a20=outv(19);
				double a21=outv(20);
				double a22=outv(21);
				double a23=outv(22);
				double a24=outv(23);
//				output.setf(ios::right);
//				output.setf(ios::scientific, ios::floatfield);
//				output<<" amat "<< "\n";
//				for(int ijk=0; ijk<outv.Size(); ijk++){
//				output << setw(30) << setprecision(10) <<outv(ijk);
//				output << "\n";
//				}
//				output.flush();
				Vector bbb=theSOE->getB();
				double b1=bbb(0);
				double b2=bbb(1);
				double b3=bbb(2);
				double b4=bbb(3);
				double b5=bbb(4);
				double b6=bbb(5);
//				output<<" bvec "<< "\n";
//				for(int ijk=0; ijk<bbb.Size(); ijk++){
//				output << setw(30) << setprecision(10) <<bbb(ijk);
//				output << "\n";
//				}


				//Vector bbb=theSOE->getB();
				theSOE->solve();
				// Save 'v' to the nodes for a "sensNodeDisp node? dof?" command
				Vector ccc=theSOE->getX();
				double c1=ccc(0);
				double c2=ccc(1);
				double c3=ccc(2);
				double c4=ccc(3);
				double c5=ccc(4);
				double c6=ccc(5);


//				output<<" bvec "<< "\n";
//				for(int ijk=0; ijk<ccc.Size(); ijk++){
//				output << setw(30) << setprecision(10) <<ccc(ijk);
//				output << "\n";
//				}
//				output.flush();

                */

				theSensitivityIntegrator->saveSensitivity( theSOE->getX(), gradNumber, numGrads );
			}else{
				done=true;
//				theSensitivityIntegrator->updateGradNumber(gradNumber);
//				theSensitivityIntegrator->saveSensitivity( *zeroVector, gradNumber, numGrads );
			}
			for(int i=0;i<numGradPositioner[idtemp];i++){
				int ii=idGradPositioner[idtemp][i];
				//gradRVPositioner[ii]->activate(false);
			}
			// Commit unconditional history variables (also for elastic problems; strain sens may be needed anyway)
			theSensitivityIntegrator->commitSensitivity(gradNumber, numGrads);
			if(done) break;
		} // End loop for each gradient
	}else{
		int act;
		bool done=false;
		for (gradNumber=1; gradNumber<=numGrads; gradNumber++ )  {
			int idtemp=gradNumber-1;
			bool activeExist=false;
			for(int i=0;i<numGradPositioner[idtemp];i++){
				int ii=idGradPositioner[idtemp][i];
				//act=gradParaPositioner[ii]->activate(true);
				if(act!=0)activeExist=true; 
			}
			if(activeExist){
				theSOE->zeroB();
				// Form new right-hand side
				theSensitivityIntegrator->formSensitivityRHS(gradNumber);
				// Solve the system of equation with the new right-hand side
				theSOE->solve();
				// Save 'v' to the nodes for a "sensNodeDisp node? dof?" command
				theSensitivityIntegrator->saveSensitivity( theSOE->getX(), gradNumber, numGrads );
			}else{
				done=true;
//				theSensitivityIntegrator->updateGradNumber(gradNumber);
//				theSensitivityIntegrator->saveSensitivity( *zeroVector, gradNumber, numGrads );
			}
			for(int i=0;i<numGradPositioner[idtemp];i++){
				int ii=idGradPositioner[idtemp][i];
				//gradParaPositioner[ii]->activate(false);
			}
			theSensitivityIntegrator->commitSensitivity(gradNumber, numGrads);
			if(done) break;
		}
	} // End loop for each gradient



    return 0;
}

bool 
NewSensitivityAlgorithm::shouldComputeAtEachStep(void)
{
	if (analysisTypeTag==1 || analysisTypeTag==2) {
		return true;
	}
	else {
		return false;
	}
}


////////////////////////////////////////////////////////////////////////////////
/////////////////  ADDED by K FUJIMURA /////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
int 
NewSensitivityAlgorithm::sensitivityDomainChanged(void)
{
  /*
	if(gradRVPositioner!=0){
		for(int i=0;i<numPos;i++) gradRVPositioner[i]=0;
		delete [] gradRVPositioner;
		gradRVPositioner=0;
	}
	if(gradParaPositioner!=0){
		for(int i=0;i<numPos;i++) gradParaPositioner[i]=0;
		delete [] gradParaPositioner;
		gradParaPositioner=0;
	}
  */
	if(idGradPositioner!=0){
		for(int i=0;i<numGrads;i++) delete [] idGradPositioner[i];
		delete [] idGradPositioner;
		idGradPositioner=0;
	}
	if(numGradPositioner!=0){
		delete [] numGradPositioner;
		numGradPositioner=0;
	}
	//RandomVariablePositioner *theRandomVariablePositioner;
	//ParameterPositioner *theParameterPositioner;

	theSOE = theAlgorithm->getLinearSOEptr();
	theIncInt = theAlgorithm->getIncrementalIntegratorPtr();
	// Get number of random variables and random variable positioners
	if (analysisTypeTag==1 || analysisTypeTag==3) {
		numGrads = theReliabilityDomain->getNumberOfRandomVariables();
		//numPos = theReliabilityDomain->getNumberOfRandomVariablePositioners();
		numPos = 0;
		//gradRVPositioner = new RandomVariablePositioner*[numPos];
		for(int i=0;i<numPos;i++){
		  //theRandomVariablePositioner = theReliabilityDomain->getRandomVariablePositionerPtr(i+1);
		  //gradRVPositioner[i]=theRandomVariablePositioner;
		}
	}
	else {
	  //numPos = theReliabilityDomain->getNumberOfParameterPositioners();
	  numPos = 0;
		numGrads = numPos;
		//gradParaPositioner = new ParameterPositioner*[numPos];
		for(int i=0;i<numPos;i++){
		  //theParameterPositioner = theReliabilityDomain->getParameterPositionerPtr(i+1);
		  //gradParaPositioner[i]=theParameterPositioner;
		}
	}

	opserr << "FATAL NewSensitivityAlgorithm::TALK TO MHS\n";
	exit(0);
//	theSensitivityIntegrator->sensitivityDomainChanged(numGrads);//Abbas...............................

	if (analysisTypeTag==1 || analysisTypeTag==3) {
		// inactivate all positioner //
		for (int i=1; i<=numPos; i++ ) {
		  //theRandomVariablePositioner = theReliabilityDomain->getRandomVariablePositionerPtr(i);
		  //theRandomVariablePositioner->activate(false);
		}
		int* itemp=new int[numGrads];
		numGradPositioner = new int[numGrads];
		idGradPositioner = new int*[numGrads];
		int numEach;
		for (int gradNumber=1; gradNumber<=numGrads; gradNumber++ )  {
			numEach=0;
			/*
			RandomVariablePositionerIter rvpIter = theReliabilityDomain->getRandomVariablePositioners();
			//for (int i=1; i<=numPos; i++ ) {
			
			while ((theRandomVariablePositioner = rvpIter()) != 0) {
				//theRandomVariablePositioner = theReliabilityDomain->getRandomVariablePositionerPtr(i);
				int rvNumber = theRandomVariablePositioner->getRvIndex();
				if ( rvNumber==gradNumber ) {
					itemp[numEach]=rvNumber;
					numEach++;
				}
			}
			*/
			numGradPositioner[gradNumber-1]=numEach;
			idGradPositioner[gradNumber-1]=new int[numEach];
			for(int i=0; i<numEach; i++) 
				idGradPositioner[gradNumber-1][i]=itemp[i];
		}
		delete [] itemp;
		itemp=0;
	}else{
	  /*
		ParameterPositionerIter pIter = theReliabilityDomain->getParameterPositioners();
		while ((theParameterPositioner = pIter()) != 0) {
		//for (int i=1; i<=numPos; i++ ) {
			//theParameterPositioner = theReliabilityDomain->getParameterPositionerPtr(i);
			theParameterPositioner->activate(false);
		}
		int* itemp=new int[numGrads];
		numGradPositioner = new int[numGrads];
		idGradPositioner = new int*[numGrads];
		int numEach;
		for (int gradNumber=1; gradNumber<=numGrads; gradNumber++ )  {
			numEach=0;
			for (int i=1; i<=numPos; i++ ) {
				itemp[numEach]=i-1;
				numEach++;
			}
			numGradPositioner[gradNumber-1]=numEach;
			idGradPositioner[gradNumber-1]=new int[numEach];
			for(int i=0; i<numEach; i++) 
				idGradPositioner[gradNumber-1][i]=itemp[i];
		}
		delete [] itemp;
		itemp=0;
	  */
	}
	return 0;
}

//int 
//NewSensitivityAlgorithm::computeSensitivities()
//{
//		opserr << "Fatal NewSensitivityAlgorithm::computeSensitivities()";
//		opserr << "This function should not be called from this object\n";
//		opserr << "This function should be called from OrigSensitivityAlgorithm\n";
//		return -1;
//}/

bool 
NewSensitivityAlgorithm::newAlgorithm(void)
{
	return true;
}
