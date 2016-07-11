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

// $Revision: 1.1 $
// $Date: 2008-02-29 19:43:52 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/analysis/telm/Analyzer.cpp,v $
                                                                     


#include <Analyzer.h>

Analyzer::Analyzer
		 (ReliabilityDomain* passedReliabilityDomain,
		  Domain* passedStructuralDomain,
		  InitialStaticAnalysis* passedInitialStaticAnalysis,
		 // SensitivityAlgorithm* passedSensitivityAlgorithm,
		 Integrator* passedSensitivityAlgorithm,
	        // SensitivityIntegrator *passedSensitivityIntegrator,
            	Integrator *passedSensitivityIntegrator,

	       	 int passedNumstep,
		  double passeddelta,
		  int passedNumLoadPatterns,
		  int* passedLoadPatterns,
		  bool passedprint)
{
	//// initialize all /////
	numOrgPatterns=0;
	theOrgPatterns=0;

	theReliabilityDomain = passedReliabilityDomain;
	theDomain = passedStructuralDomain;
	theInitialStaticAnalysis = passedInitialStaticAnalysis;
	theSensitivityAlgorithm = passedSensitivityAlgorithm;
	theSensitivityIntegrator = passedSensitivityIntegrator;

	if(theSensitivityAlgorithm == 0){
		activeSensitivity=false;
	}else{
		activeSensitivity=true;
	}
	print = passedprint;
	if(print){
		output.open("Analyzer.txt", ios::out);
	}
	Numstep = passedNumstep;
	delta = passeddelta;

	NumLoadPatterns = passedNumLoadPatterns;
	if( NumLoadPatterns == 0 ){
		LoadPatterns = 0;
	}else{
		if(LoadPatterns!=0){ delete [] LoadPatterns; LoadPatterns=0;}
		LoadPatterns = new int[NumLoadPatterns];
		if(LoadPatterns==0){
			opserr << "ERROR\n";
			opserr << "in sufficient memory Analyzer::Analyzer\n";
			opserr << "allocate LoadPatterns\n";
			exit(-1);
		}
		for( int i=0; i<NumLoadPatterns; i++){
			LoadPatterns[i] = passedLoadPatterns[i];
		}
	}

	if(print){
		output << "=====SelectLoadAnalyzer=====\n";
		output << "\n";
		output << " Number of Selected Loads" << NumLoadPatterns << "\n";
		output << "\n";
		output << " Load Pattern IDs \n";
		for( int i=0; i<NumLoadPatterns; i++)
			output << LoadPatterns[i] << "\n";
	}
	output.flush();
	saveLoads();
}
Analyzer::~Analyzer()
{
	if(numOrgPatterns!=0){
		for(int i=0; i<numOrgPatterns; i++)	theOrgPatterns[i]=0;
		delete [] theOrgPatterns;
		theOrgPatterns = 0;
	}
	if( LoadPatterns != 0 ) {delete LoadPatterns;LoadPatterns = 0;}
}
void Analyzer::saveLoads()
{
	if(numOrgPatterns!=0){
		for(int i=0; i<numOrgPatterns; i++)	theOrgPatterns[i]=0;
		delete [] theOrgPatterns;
		theOrgPatterns=0;
	}

	numOrgPatterns=0;
	LoadPattern* thePattern;
	LoadPatternIter& thePatterns = theDomain->getLoadPatterns();	
	while((thePattern = thePatterns()) != 0){
		numOrgPatterns++;
	}
	if(theOrgPatterns!=0){
		for( int i=0; i<numOrgPatterns; i++){
//			delete theOrgPatterns[i];
			theOrgPatterns[i]=0;
		}
		delete [] theOrgPatterns;
		theOrgPatterns=0;
	}
	theOrgPatterns = new LoadPattern*[numOrgPatterns];
	if(theOrgPatterns==0){
		opserr << "ERROR\n";
		opserr << "in sufficient memory Analyzer::saveLoads\n";
		opserr << "allocate theOrgPatterns\n";
		exit(-1);
	}
	int itemp=0;
	thePatterns.reset();
	while((thePattern = thePatterns()) != 0){
		theOrgPatterns[itemp]=thePattern;
		itemp++;
	}
}
void Analyzer::modifyLoads()
{
	LoadPattern* thePattern;
	LoadPattern* thePat;
	if(NumLoadPatterns != 0){
		LoadPatternIter& thePatterns = theDomain->getLoadPatterns();	
		///// remove all /////
		while((thePattern = thePatterns()) != 0){
			int tag=thePattern->getTag();
			thePat=theDomain->removeLoadPattern(tag);
		}
		for( int i=0; i<numOrgPatterns; i++){
			thePattern=theOrgPatterns[i];
			int tag=thePattern->getTag();
			bool found=false;
			for(int i=0; i<NumLoadPatterns; i++){
				if(tag == LoadPatterns[i]){
					found=true;
					break;
				}
			}
			if(found) theDomain->addLoadPattern(thePattern);
		}
	}

	if(print){
		output << "\n";
		output << " after modify load  in theDomain \n";
		output << "\n";
		LoadPatternIter& thePatterns = theDomain->getLoadPatterns();	
		while((thePattern = thePatterns()) != 0){
			int tag=thePattern->getTag();
			output << " load pattern " << tag <<"\n";
		}
		output << "\n";
		output << " after modify load in theOrgPatterns \n";
		output << "\n";
		for( int i=0; i<numOrgPatterns; i++){
			thePattern=theOrgPatterns[i];
			int tag=thePattern->getTag();
			output << " load pattern " << tag <<"\n";
		}
	}
}	

void Analyzer::recoverLoads()
{
	LoadPattern* thePattern;
	LoadPattern* thePat;
	LoadPatternIter& thePatterns = theDomain->getLoadPatterns();	
	///// remove all currently in the domain  /////
	while((thePattern = thePatterns()) != 0){
		int tag=thePattern->getTag();
		thePat=theDomain->removeLoadPattern(tag);
	}
	//// add all original /////
	for( int i=0; i<numOrgPatterns; i++){
		thePattern=theOrgPatterns[i];
		theDomain->addLoadPattern(thePattern);
	}

	if(print){
		output << "\n";
		output << " after recover load  in theDomain \n";
		output << "\n";
		LoadPatternIter& thePatterns = theDomain->getLoadPatterns();	
		while((thePattern = thePatterns()) != 0){
			int tag=thePattern->getTag();
			output << " load pattern " << tag <<"\n";
		}
	}
}
void Analyzer::printresults()
{
	output<<"\n";
	output<<" ---- result of analysis -----\n";
	output<<"\n";
	output.setf(ios::right);
	output.setf(ios::scientific, ios::floatfield);
	NodeIter& theNodes = theDomain->getNodes();
	Node* theNode;
	while((theNode = theNodes()) != 0){
		int tag=theNode->getTag();
		int ndof=theNode->getNumberDOF();
		Vector Disp=theNode->getDisp();
		output << setw(10) << tag;
		for(int i=0; i<ndof; i++)
			output << setw(15) << setprecision(5) << Disp(i);
		output<<"\n";
	}
	output.flush();
}
void Analyzer::setNsteps(int passedNsteps)
{
	Numstep = passedNsteps;
}
double Analyzer::getDt()
{
	return delta;
}	
int Analyzer::getNstep()
{
 	return Numstep;
}	
