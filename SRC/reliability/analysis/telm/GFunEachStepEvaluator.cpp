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
// $Date: 2008-03-13 22:26:49 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/analysis/telm/GFunEachStepEvaluator.cpp,v $

#include <GFunEachStepEvaluator.h>
#include <ReliabilityDomain.h>

GFunEachStepEvaluator::GFunEachStepEvaluator(Tcl_Interp *passedTclInterp,
											 ReliabilityDomain *passedReliabilityDomain,
											 Domain* passedDomain,
											 int passednumSteps,
											 bool passedprint)
{
	theTclInterp=passedTclInterp;
	theDomain=passedDomain;
	print=passedprint;
	numRVs=passedReliabilityDomain->getNumberOfRandomVariables();
	numLSF=passedReliabilityDomain->getNumberOfLimitStateFunctions();
	numSteps=passednumSteps;

	if(print){
		output.open("GFunEachStepEvaluator.txt", ios::out);
		output << "\n";
		output << "GFunEachStepEvaluator\n";
		output << "\n";
		output << "numLSF  "<<numLSF<<"\n";
		output << "numSteps"<<numSteps<<"\n";
		output << "\n";
		output.flush();
	}

	theLSFValues=new Matrix(numLSF, numSteps);
	if(theLSFValues == 0) {
		 opserr << " memory allocation error 5 in GFunEachStepEvaluator \n";
		 exit(-1);
	}
	theConvFlag=new Matrix(numLSF, numSteps);
	if(theConvFlag == 0) {
		 opserr << " memory allocation error 5 in GFunEachStepEvaluator \n";
		 exit(-1);
	}


	PFthershold= new Vector(numLSF);

	thePFCoeffs=0;
	thePFCoeffs= new TaggedObjectStorage*[numLSF];
	if(thePFCoeffs==0){
		opserr << "GFunEachStepEvaluator::setLimitState - out of memory "
			   << "for PerfomanceFunctionCoeff \n";
		exit(-1);}
	for(int i=0; i<numLSF; i++){
		thePFCoeffs[i]=new ArrayOfTaggedObjects(32);
		if(thePFCoeffs[i]==0){
		opserr << "GFunEachStepEvaluator::setLimitState - out of memory "
			   << "for PerfomanceFunctionCoeff \n";
		exit(-1);}
	}

	thePFIters=0;
	thePFIters = new PerformanceFunctionCoefficientIter*[numLSF];
	if(thePFIters==0){
		opserr << "GFunEachStepEvaluator::setLimitState - out of memory "
			   << "for thePFIters\n";
		exit(-1);}
	for(int i=0; i<numLSF; i++){
		thePFIters[i]=0;
	}

	analyzeLSF(passedReliabilityDomain);

}
GFunEachStepEvaluator::~GFunEachStepEvaluator()
{
	if(PFthershold!=0){ delete PFthershold; PFthershold=0; }
	if(thePFCoeffs!=0){
		for(int i=0;i<numLSF; i++){
			delete thePFCoeffs[i]; thePFCoeffs[i]=0;
		}
		delete [] thePFCoeffs;
		thePFCoeffs=0;
	}
	if(thePFIters!=0){
		for(int i=0;i<numLSF; i++){
			delete thePFIters[i]; thePFIters[i]=0;
		}
		delete [] thePFIters;
		thePFIters=0;
	}
	if(theLSFValues!=0){
		delete theLSFValues; theLSFValues=0;
	}
}
void GFunEachStepEvaluator::analyzeLSF(ReliabilityDomain* theRelib)
{
	for( int lsf=0; lsf<numLSF; lsf++){
		LimitStateFunction* theLSF = theRelib->getLimitStateFunctionPtr(lsf+1);
		// Get the limit-state function expression
		TaggedObjectStorage* TempPerformFuncCoeffs=0;
		PerformanceFunctionCoefficientIter* TempPfCoeffIter=0;
		(*PFthershold)(lsf)=setLimitState(theLSF,lsf);
//			TempPfCoeffIter,TempPfCoeffIter);
//			thePFCoeffs[lsf]=TempPerformFuncCoeffs;
//			thePFIters[lsf]=TempPfCoeffIter;
		if(print){
			output.setf( ios::scientific );
			output << "Performance Function Constant.."<<lsf<<"\n";
			output << setw(15) << setprecision(5) << (*PFthershold)(lsf) << "\n";
			output << "\n";
			output << "      Node";
			output << " Direction";
			output << "     Coefficient\n";
			thePFIters[lsf]->reset();
			PerformanceFunctionCoeff* thePfCoeff;
			while((thePfCoeff = (*thePFIters[lsf])()) != 0){
				int N=thePfCoeff->getNodeID();
				int D=thePfCoeff->getDirection();
				double c=thePfCoeff->getCoefficient();
				output<< setw(10) << N;
				output<< setw(10) << D;
				output<< setw(15) << c;
			}
			output.flush();
		}
	}
}
double GFunEachStepEvaluator::setLimitState
							 (LimitStateFunction* theLSF, int lsf)
//							  TaggedObjectStorage*& thePerformFuncCoeffs,
//							  PerformanceFunctionCoefficientIter*& thePfCoeffIter)
{
    PerformanceFunctionCoeff* thePFCoeff;

	const char *theExpression = theLSF->getExpression();
	
	// This parsing should go away -- MHS 10/7/2011
	//char *theTokExpression = theLSF->getTokenizedExpression();
	char *theTokExpression = ""; // So that it compiles

	char separators[5] = "}{";
	char *dollarSign = "$";
	char *underscore = "_";
	char tempchar[100]="";

	char lsf_forTokenizing[500];
	strcpy(lsf_forTokenizing,theExpression);
	char tclAssignment[100];

//  1st set all disp to zero and evaluate to have the constant term //
	char *tokenPtr = strtok( lsf_forTokenizing, separators);
	while ( tokenPtr != NULL ) {
		strcpy(tempchar,tokenPtr);
		if ( strncmp(tokenPtr, "ud", 2) == 0) {
			opserr << "Performance Function need to be \n";
			opserr << "a linear function of disp to apply \n";
			opserr << "initial shape analysys \n";
			opserr << "skip initial shape for this lsf \n";
			exit(-1);
		}
		// If a nodal displacement is detected
		else if ( strncmp(tokenPtr, "u", 1) == 0) {
			int nodeNumber, direction;
			sscanf(tempchar,"u_%i_%i", &nodeNumber, &direction);
			sprintf(tclAssignment,"set u_%d_%d 0.0",nodeNumber,direction);
			Tcl_Eval( theTclInterp, tclAssignment);
		}
		else if ( strncmp(tokenPtr, "rec",3) == 0) {
			opserr << "Performance Function need to be \n";
			opserr << "a linear function of disp to apply \n";
			opserr << "initial shape analysys \n";
			opserr << "skip initial shape for this lsf \n";
			exit(-1);
		}
		tokenPtr = strtok( NULL, separators);
	}
	double pfthreshold = 0.0;
	double dthreshold = 0.0;
	Tcl_ExprDouble( theTclInterp, theTokExpression, &pfthreshold );

	// Re-create possible recorders for subsequent analyses
	strcpy(lsf_forTokenizing,theExpression);
	tokenPtr = strtok( lsf_forTokenizing, separators);
	while ( tokenPtr != NULL ) {
		strcpy(tempchar,tokenPtr);

        if ( strncmp(tokenPtr, "u", 1) == 0) {
			int nodeNumber, direction;
			sscanf(tempchar,"u_%i_%i", &nodeNumber, &direction);
			sprintf(tclAssignment,"set u_%d_%d 1.0",nodeNumber,direction);
			Tcl_Eval( theTclInterp, tclAssignment);
			Tcl_ExprDouble( theTclInterp, theTokExpression, &dthreshold );
			double coeff=dthreshold-pfthreshold;
			int Tag=thePFCoeffs[lsf]->getNumComponents()+1;
            thePFCoeff=new PerformanceFunctionCoeff
		     (Tag, nodeNumber, direction, coeff);
			bool result = thePFCoeffs[lsf]->addComponent(thePFCoeff);
			if(!result){
		      opserr << "OutCrossingAnalysis::AnalyzeGfun - out of memory3\n";
			exit(-1);}
			sprintf(tclAssignment,"set u_%d_%d 0.0",nodeNumber,direction);
			Tcl_Eval( theTclInterp, tclAssignment);
		}
		
		tokenPtr = strtok( NULL, separators);
	}
	thePFIters[lsf]=new PerformanceFunctionCoefficientIter(thePFCoeffs[lsf]);
	if(thePFIters[lsf]==0){
       opserr << "GFunEachStepEvaluator::setLimitState memory3\n";
	   exit(1);
	}
	return pfthreshold;
}
void GFunEachStepEvaluator::initialize()
{
	for(int i=0; i<numLSF; i++){
		for( int j=0; j<numSteps; j++) {
			(*theLSFValues)(i,j)=0.0;
			(*theConvFlag)(i,j)=0;
		}
	}
}
void GFunEachStepEvaluator::evaluateG(int istep)
{
	Node* theNode =0;
	PerformanceFunctionCoeff* thePfCoeff =0 ;
	for( int lsf=0; lsf<numLSF; lsf++){
		thePFIters[lsf]->reset();
		double sum=0.0;
		while((thePfCoeff = (*thePFIters[lsf])()) != 0){
		    int N=thePfCoeff->getNodeID();
			theNode = theDomain->getNode(N);
		    int D=thePfCoeff->getDirection();
		    double c=thePfCoeff->getCoefficient();
			sum+=(theNode->getDisp())(D-1)*c;
		}
		(*theLSFValues)(lsf,istep)=sum;
		(*theConvFlag)(lsf,istep)=1.0;

	}
}
void GFunEachStepEvaluator::evaluateTrialG(int istep)
{
	Node* theNode =0;
	PerformanceFunctionCoeff* thePfCoeff =0 ;
	for( int lsf=0; lsf<numLSF; lsf++){
		thePFIters[lsf]->reset();
		double sum=0.0;
		while((thePfCoeff = (*thePFIters[lsf])()) != 0){
		    int N=thePfCoeff->getNodeID();
			theNode = theDomain->getNode(N);
		    int D=thePfCoeff->getDirection();
		    double c=thePfCoeff->getCoefficient();
			sum+=(theNode->getTrialDisp())(D-1)*c;
		}
		(*theLSFValues)(lsf,istep)=sum;
		(*theConvFlag)(lsf,istep)=0.0;

	}
}

