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
                                                                        
// $Revision: 1.2 $
// $Date: 2008-05-27 20:04:30 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/analysis/telm/AnalyzerGradGEvaluator.cpp,v $

//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu)
//

#include <AnalyzerGradGEvaluator.h>
#include <Vector.h>
#include <Matrix.h>
#include <GradGEvaluator.h>
#include <GFunEvaluator.h>
#include <ReliabilityDomain.h>
#include <LimitStateFunction.h>
#include <RandomVariable.h>
#include <tcl.h>
#include <string.h>

#include <fstream>
#include <iomanip>
#include <iostream>
using std::ifstream;
using std::ios;
using std::setw;
using std::setprecision;
using std::setiosflags;


AnalyzerGradGEvaluator::AnalyzerGradGEvaluator(
					Tcl_Interp *passedTclInterp,
					ReliabilityDomain *passedReliabilityDomain,
                    Domain* passedDomain,
					FunctionEvaluator* passedGFunEvaluator,
					bool PdoGradientCheck)
:GradientEvaluator(passedReliabilityDomain, passedGFunEvaluator)
{
	
	doGradientCheck = PdoGradientCheck;
	theDomain=passedDomain;

	bool fdf=false;
	this->setfinitedifference(fdf);


	int nrv = passedReliabilityDomain->getNumberOfRandomVariables();
	grad_g = new Vector(nrv);
	grad_g_matrix = 0;

	DgDdispl = 0;
}

AnalyzerGradGEvaluator::~AnalyzerGradGEvaluator()
{
	if (grad_g != 0) {delete grad_g;grad_g;}

	if (DgDdispl != 0){	delete DgDdispl;DgDdispl=0;}

	if (grad_g_matrix != 0){delete grad_g_matrix;grad_g_matrix=0;}
}




Vector
AnalyzerGradGEvaluator::getGradG()
{
	return (*grad_g);
}


Matrix
AnalyzerGradGEvaluator::getAllGradG()
{
	if (grad_g_matrix==0) {
		Matrix dummy(1,1);
		return dummy;
	}
	else {
		return (*grad_g_matrix);
	}
}

int
AnalyzerGradGEvaluator::computeGradG(double g, const Vector &passed_x)
{
	numberOfEvalIncSens++;
	// Zero out the previous result matrix
	if (DgDdispl != 0) {
		delete DgDdispl;
		DgDdispl = 0;
	}

	// Call base class method
	computeParameterDerivatives(g);
	int nrv = theReliabilityDomain->getNumberOfRandomVariables();

	Node* theNode =0;
	PerformanceFunctionCoeff* thePfCoeff =0 ;
	thePfCoeffIter->reset();
	int numCoeff=0;
    while((thePfCoeff = (*thePfCoeffIter)()) != 0)numCoeff++;
	DgDdispl = new Matrix(numCoeff, 3);
	double onedudx;
	Vector dudx(nrv);

	if(grad_g !=0){
		delete grad_g;
		grad_g = 0;
	}
	grad_g = new Vector(nrv);

	int idnumber=0;
	thePfCoeffIter->reset();
    while((thePfCoeff = (*thePfCoeffIter)()) != 0){
	    int N=thePfCoeff->getNodeID();
		int nodeNumber=N;
		theNode = theDomain->getNode(N);
	    int D=thePfCoeff->getDirection();
		int direction=D;
	    double c=thePfCoeff->getCoefficient();
		double onedgdu = c;
		(*DgDdispl)(idnumber,0) = (double)nodeNumber;
		(*DgDdispl)(idnumber,1) = (double)direction;
		(*DgDdispl)(idnumber,2) = onedgdu;
		for (int i=1; i<=nrv; i++) {
			onedudx=theNode->getDispSensitivity(direction,i);
			dudx( (i-1) ) = onedudx;
		}
		(*grad_g) += onedgdu*dudx;
	}

	if (doGradientCheck) {
		char myString[100];
		ofstream outputFile( "DDMgradients.out", ios::out );
		opserr << endln;
		for (int ddm=0; ddm<grad_g->Size(); ddm++) {
			opserr << "DDM("<< (ddm+1) << ") = " << (*grad_g)(ddm) << endln;
			sprintf(myString,"%20.16e ",(*grad_g)(ddm));
			outputFile << myString << endln;
		}
		outputFile.close();
		opserr << "PRESS Ctrl+C TO TERMINATE APPLICATION!" << endln;
		while(true) {
		}
	}

	return 0;

}



int
AnalyzerGradGEvaluator::computeAllGradG(const Vector &gFunValues, const Vector &passed_x)
{
	numberOfEvalIncSens++;
	// Allocate result matrix
	Vector gradG(passed_x.Size());

	if (grad_g_matrix != 0) {
		delete grad_g_matrix;
		grad_g_matrix = 0;
	}
	grad_g_matrix = new Matrix(passed_x.Size(), gFunValues.Size());


//	if (grad_g_matrix == 0) {
//		grad_g_matrix = new Matrix(passed_x.Size(), gFunValues.Size());
//	}
//	else {
//		grad_g_matrix->Zero();
//	}


	// Loop over performance functions
	for (int j=1; j<=gFunValues.Size(); j++) {

		// Set tag of active limit-state function
		theReliabilityDomain->setTagOfActiveLimitStateFunction(j);

		this->computeGradG(gFunValues(j-1),passed_x);
		gradG = this->getGradG();

		for (int i=1; i<=passed_x.Size(); i++) {
	
			(*grad_g_matrix)(i-1,j-1) = gradG(i-1);
		}
	}

	return 0;
}


Matrix 
AnalyzerGradGEvaluator::getDgDdispl()
{
	return (*DgDdispl);
}



void AnalyzerGradGEvaluator::setReliabilityDomain(ReliabilityDomain *pdomain)
{
	theReliabilityDomain=pdomain;
	int nrv = theReliabilityDomain->getNumberOfRandomVariables();
	if(	grad_g != 0){
		delete grad_g;
		grad_g = 0;
		grad_g = new Vector(nrv);
	}
	if(grad_g_matrix != 0){
		delete grad_g_matrix;
		grad_g_matrix = 0;
	}

	if( DgDdispl != 0 ){ 
		delete DgDdispl;
		DgDdispl = 0;
	}
}

