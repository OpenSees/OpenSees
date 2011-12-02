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
                                                                        
// $Revision: 1.14 $
// $Date: 2008-04-10 16:26:27 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/analysis/sensitivity/OpenSeesGradGEvaluator.cpp,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu)
//

#include <OpenSeesGradGEvaluator.h>
#include <Vector.h>
#include <Matrix.h>
#include <GradGEvaluator.h>
#include <ReliabilityDomain.h>
#include <LimitStateFunction.h>
#include <RandomVariable.h>
#include <RandomVariableIter.h>
#include <SensitivityAlgorithm.h>
#include <tcl.h>
#include <string.h>

#include <RandomVariableIter.h>
#include <LimitStateFunctionIter.h>

#include <fstream>
#include <iomanip>
#include <iostream>
using std::ifstream;
using std::ios;
using std::setw;
using std::setprecision;
using std::setiosflags;


OpenSeesGradGEvaluator::OpenSeesGradGEvaluator(Tcl_Interp *passedTclInterp,
					       ReliabilityDomain *passedReliabilityDomain,
					       SensitivityAlgorithm *theAlgo,
					       bool PdoGradientCheck)
:GradGEvaluator(passedReliabilityDomain, passedTclInterp)
{
	theReliabilityDomain = passedReliabilityDomain;
	doGradientCheck = PdoGradientCheck;
	theSensAlgo = theAlgo;

	int nrv = passedReliabilityDomain->getNumberOfRandomVariables();
	grad_g = new Vector(nrv);
	grad_g_matrix = 0;

	DgDdispl = 0;

	/////S added by K Fujimura /////
	bool fdf=true;
	this->setfinitedifference(fdf);
	/////E added by K Fujimura /////
}

OpenSeesGradGEvaluator::~OpenSeesGradGEvaluator()
{
	if (grad_g != 0) 
		delete grad_g;

	if (DgDdispl != 0)
		delete DgDdispl;

	if (grad_g_matrix != 0)
		delete grad_g_matrix;
}




Vector
OpenSeesGradGEvaluator::getGradG()
{
	return (*grad_g);
}


Matrix
OpenSeesGradGEvaluator::getAllGradG()
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
OpenSeesGradGEvaluator::computeGradG(double g, const Vector &passed_x)
{
	numberOfEvalIncSens++;///// added by K Fujimura /////

	// Zero out the previous result matrix
	if (DgDdispl != 0) {
		delete DgDdispl;
		DgDdispl = 0;
	}

	// Call base class method
	computeParameterDerivatives(g);

	// Initial declaractions
	double perturbationFactor = 0.001; // (is multiplied by stdv and added to others...)
	char tclAssignment[1000];
	char *dollarSign = "$";
	char *underscore = "_";
	char lsf_expression[1000];
	char separators[5] = "}{";
	int nrv = theReliabilityDomain->getNumberOfRandomVariables();
	RandomVariable *theRandomVariable;
	char tempchar[1000];
	char newSeparators[5] = "_";
	double g_perturbed;
	int i;
	double onedudx;
	Vector dudx(nrv);
	

	// Compute gradients if this is a path-INdependent analysis
	// (This command only has effect if it IS path-independent.)
	if (theSensAlgo != 0 && !(theSensAlgo->shouldComputeAtEachStep()) )
	  theSensAlgo->computeSensitivities();

	// Initialize gradient vector
	grad_g->Zero();


	// "Download" limit-state function from reliability domain
	int lsf = theReliabilityDomain->getTagOfActiveLimitStateFunction();
	LimitStateFunction *theLimitStateFunction = 
		theReliabilityDomain->getLimitStateFunctionPtr(lsf);
	char *theExpression = theLimitStateFunction->getExpression();
	char lsf_copy[1000];
	strcpy(lsf_copy,theExpression);


	// Tokenize the limit-state function and COMPUTE GRADIENTS
	char *tokenPtr = strtok( lsf_copy, separators); 
	while ( tokenPtr != NULL ) {

		strcpy(tempchar,tokenPtr);

		if ( strncmp(tokenPtr, "x",1) == 0) {

			// Get random variable tag
			int rvTag;
			sscanf(tempchar,"x_%i",&rvTag);

			// Perturb its value according to its standard deviation
			theRandomVariable = theReliabilityDomain->getRandomVariablePtr(rvTag);

			double stdv = theRandomVariable->getStdv();
			//int rvIndex = theRandomVariable->getIndex();
			int rvIndex = theReliabilityDomain->getRandomVariableIndex(rvTag);

			sprintf(tclAssignment , "set x_%d  %35.20f", rvTag, (passed_x(rvIndex)+perturbationFactor*stdv) );
			if (Tcl_Eval( theTclInterp, tclAssignment) == TCL_ERROR) {
			  opserr << "ERROR OpenSeesGradGEvaluator -- Tcl_Eval returned error" << endln;
			  opserr << theTclInterp->result << endln;
			  return -1;
			}

			// Evaluate limit-state function again
			char *theTokenizedExpression = theLimitStateFunction->getTokenizedExpression();
			if (Tcl_ExprDouble( theTclInterp, theTokenizedExpression, &g_perturbed ) == TCL_ERROR) {
			  opserr << "ERROR OpenSeesGradGEvaluator -- Tcl_ExprDouble returned error" << endln;
			  opserr << theTclInterp->result << endln;
			  return -1;
			}

			// Make assignment back to its original value
			sprintf(tclAssignment , "set x_%d  %35.20f", rvTag, passed_x(rvIndex) );
			if (Tcl_Eval( theTclInterp, tclAssignment) == TCL_ERROR) {
			  opserr << "ERROR OpenSeesGradGEvaluator -- Tcl_Eval returned error" << endln;
			  opserr << theTclInterp->result << endln;
			  return -1;
			}

			// Add gradient contribution
			(*grad_g)(rvIndex) += (g_perturbed-g)/(perturbationFactor*stdv);
		}
		// If a nodal velocity is detected
		else if ( strncmp(tokenPtr, "ud", 2) == 0) {

			// Get node number and dof number
			int nodeNumber, direction;
			sscanf(tempchar,"ud_%i_%i", &nodeNumber, &direction);

			// Keep the original value
			double originalValue;
			sprintf(tclAssignment,"$ud_%d_%d", nodeNumber, direction);
			if (Tcl_ExprDouble( theTclInterp, tclAssignment, &originalValue) == TCL_ERROR) {
			  opserr << "ERROR OpenSeesGradGEvaluator -- Tcl_ExprDouble returned error" << endln;
			  opserr << theTclInterp->result << endln;
			  return -1;
			}

			// Set perturbed value in the Tcl workspace

			// ---------------------- Quan and michele 
			double newValue;
			
			//if (originalValue ==0 ) newValue = 0.0001;
			if (fabs(originalValue)<1.0e-14 ) newValue = 0.0001;  //---  Quan Jan 2007 ---

			else newValue = originalValue*(1.0+perturbationFactor);

			sprintf(tclAssignment,"set ud_%d_%d %35.20f", nodeNumber, direction, newValue);
			if (Tcl_Eval( theTclInterp, tclAssignment) == TCL_ERROR) {
			  opserr << "ERROR OpenSeesGradGEvaluator -- Tcl_Eval returned error" << endln;
			  opserr << theTclInterp->result << endln;
			  return -1;
			}
			// Evaluate the limit-state function again
			char *theTokenizedExpression = theLimitStateFunction->getTokenizedExpression();
			if (Tcl_ExprDouble( theTclInterp, theTokenizedExpression, &g_perturbed ) == TCL_ERROR) {
			  opserr << "ERROR OpenSeesGradGEvaluator -- Tcl_ExprDouble returned error" << endln;
			  opserr << theTclInterp->result << endln;
			  return -1;
			}

			// Compute gradient
			double onedgdu = (g_perturbed-g)/(newValue-originalValue);

			//double onedgdu = (g_perturbed-g)/(originalValue*perturbationFactor);

			// Make assignment back to its original value
			sprintf(tclAssignment,"set ud_%d_%d %35.20f", nodeNumber, direction, originalValue);
			if (Tcl_Eval( theTclInterp, tclAssignment) == TCL_ERROR) {
			  opserr << "ERROR OpenSeesGradGEvaluator -- Tcl_Eval returned error" << endln;
			  opserr << theTclInterp->result << endln;
			  return -1;
			}
			// Obtain DDM gradient vector
			RandomVariableIter &rvIter = 
			  theReliabilityDomain->getRandomVariables();
			RandomVariable *theRV;
			//for (int i=1; i<=nrv; i++) {
			while ((theRV = rvIter()) != 0) {
			  int rvTag = theRV->getTag();
			  //int i = theRV->getIndex();
			  int i = theReliabilityDomain->getRandomVariableIndex(rvTag);
			  sprintf(tclAssignment , "set sens [sensNodeVel %d %d %d ]",nodeNumber,direction,i);
			  if (Tcl_Eval( theTclInterp, tclAssignment) == TCL_ERROR) {
			    opserr << "ERROR OpenSeesGradGEvaluator -- Tcl_Eval returned error" << endln;
			    opserr << theTclInterp->result << endln;
			    return -1;
			  }
			  sprintf(tclAssignment , "$sens ");
			  if (Tcl_ExprDouble( theTclInterp, tclAssignment, &onedudx ) == TCL_ERROR) {
			    opserr << "ERROR OpenSeesGradGEvaluator -- Tcl_ExprDouble returned error" << endln;
			    opserr << theTclInterp->result << endln;
			    return -1;
			  }
			  dudx(i) = onedudx;
			}

			// Add gradient contribution
			(*grad_g) += onedgdu*dudx;

		}
		// If a nodal displacement is detected
		else if ( strncmp(tokenPtr, "u", 1) == 0) {

			// Get node number and dof number
			int nodeNumber, direction;
			sscanf(tempchar,"u_%i_%i", &nodeNumber, &direction);

			// Keep the original value


			double originalValue;
			sprintf(tclAssignment,"$u_%d_%d", nodeNumber, direction);
			if (Tcl_ExprDouble( theTclInterp, tclAssignment, &originalValue) == TCL_ERROR) {
			  opserr << "ERROR OpenSeesGradGEvaluator -- Tcl_ExprDouble returned error" << endln;
			  opserr << theTclInterp->result << endln;
			  return -1;
			}

			// Set perturbed value in the Tcl workspace

			// ---------------------- Quan and michele   

			double newValue;
			//if (originalValue ==0 ) newValue = 0.0001;
			if (fabs(originalValue)<1.0e-14 ) newValue = 0.0001;  //---  Quan Jan 2006 ---
			else 	newValue = originalValue*(1.0+perturbationFactor);

			sprintf(tclAssignment,"set u_%d_%d %35.20f", nodeNumber, direction, newValue);
			if (Tcl_Eval( theTclInterp, tclAssignment) == TCL_ERROR) {
			  opserr << "ERROR OpenSeesGradGEvaluator -- Tcl_Eval returned error" << endln;
			  opserr << theTclInterp->result << endln;
			  return -1;
			}

			// Evaluate the limit-state function again
			char *theTokenizedExpression = theLimitStateFunction->getTokenizedExpression();
			if (Tcl_ExprDouble( theTclInterp, theTokenizedExpression, &g_perturbed ) == TCL_ERROR) {
			  opserr << "ERROR OpenSeesGradGEvaluator -- Tcl_ExprDouble returned error" << endln;
			  opserr << theTclInterp->result << endln;
			  return -1;
			}

			// Compute gradient


			double onedgdu = (g_perturbed-g)/(newValue - originalValue);

			// Store the DgDdispl in a matrix
			if (DgDdispl == 0) {
				DgDdispl = new Matrix(1, 3);
				(*DgDdispl)(0,0) = (double)nodeNumber;
				(*DgDdispl)(0,1) = (double)direction;
				(*DgDdispl)(0,2) = onedgdu;
			}
			else {
				int oldSize = DgDdispl->noRows();
				Matrix tempMatrix = *DgDdispl;
				delete DgDdispl;
				DgDdispl = new Matrix(oldSize+1, 3);
				for (i=0; i<oldSize; i++) {
					(*DgDdispl)(i,0) = tempMatrix(i,0);
					(*DgDdispl)(i,1) = tempMatrix(i,1);
					(*DgDdispl)(i,2) = tempMatrix(i,2);
				}
				(*DgDdispl)(oldSize,0) = (double)nodeNumber;
				(*DgDdispl)(oldSize,1) = (double)direction;
				(*DgDdispl)(oldSize,2) = onedgdu;
			}


			// Make assignment back to its original value
			sprintf(tclAssignment,"set u_%d_%d %35.20f", nodeNumber, direction, originalValue);
			if (Tcl_Eval( theTclInterp, tclAssignment) == TCL_ERROR) {
			  opserr << "ERROR OpenSeesGradGEvaluator -- Tcl_Eval returned error" << endln;
			  opserr << theTclInterp->result << endln;
			  return -1;
			}
			// Obtain DDM gradient vector
			RandomVariableIter &rvIter = 
			  theReliabilityDomain->getRandomVariables();
			RandomVariable *theRV;
			//for (int i=1; i<=nrv; i++) {
			while ((theRV = rvIter()) != 0) {
			  int rvTag = theRV->getTag();
			  //int i = theRV->getIndex();
			  int i = theReliabilityDomain->getRandomVariableIndex(rvTag);
			  sprintf(tclAssignment , "set sens [sensNodeDisp %d %d %d ]",nodeNumber,direction,i);
			  if (Tcl_Eval( theTclInterp, tclAssignment) == TCL_ERROR) {
			    opserr << theTclInterp->result << endln;
			    opserr << "ERROR OpenSeesGradGEvaluator -- Tcl_Eval returned error" << endln;
			    return -1;
			  }
			  sprintf(tclAssignment , "$sens ");
			  if (Tcl_ExprDouble( theTclInterp, tclAssignment, &onedudx ) == TCL_ERROR) {
			    opserr << "ERROR OpenSeesGradGEvaluator -- Tcl_ExprDouble returned error" << endln;
			    opserr << theTclInterp->result << endln;
			    return -1;
			  }
			  dudx(i) = onedudx;
			}

			//opserr <<"OSGGE: " << nodeNumber << ' ' << direction << endln;;
			//opserr <<"OSGGE u_: " << onedgdu << ' ' << nrv << ' ' << dudx ;

			// Add gradient contribution
			(*grad_g) += onedgdu*dudx;

		}
		else if ( strncmp(tokenPtr, "rec_element",11) == 0) {

		  // Initial declarations
		  char restString[100];
		  
		  // Start obtaining information about element number etc. 
		  int eleNumber;
		  sscanf(tempchar,"rec_element_%i_%s", &eleNumber, restString);
		  
		  if ( strncmp(restString, "section",7) == 0) {
		    int sectionNumber;
		    int rowNumber;
		    sscanf(restString,"section_%i_%s", &sectionNumber, restString);
		    if ( strncmp(restString, "force",5) == 0) {
		      sscanf(restString,"force_%i", &rowNumber);

		      // Keep the original value
		      double originalValue;
		      sprintf(tclAssignment,"$rec_element_%d_section_%d_force_%d", eleNumber, sectionNumber, rowNumber);
		      if (Tcl_ExprDouble( theTclInterp, tclAssignment, &originalValue) == TCL_ERROR) {
			opserr << "ERROR OpenSeesGradGEvaluator -- Tcl_ExprDouble returned error" << endln;
			opserr << theTclInterp->result << endln;
			return -1;
		      }
		      
		      // Set perturbed value in the Tcl workspace
		      double newValue = originalValue*(1.0+perturbationFactor);
		      sprintf(tclAssignment,"set rec_element_%d_section_%d_force_%d %35.20f", eleNumber, sectionNumber, rowNumber, newValue);
		      if (Tcl_Eval( theTclInterp, tclAssignment) == TCL_ERROR) {
			opserr << "ERROR OpenSeesGradGEvaluator -- Tcl_Eval returned error" << endln;
			opserr << theTclInterp->result << endln;
			return -1;
		      }

		      // Evaluate the limit-state function again
		      char *theTokenizedExpression = theLimitStateFunction->getTokenizedExpression();
		      if (Tcl_ExprDouble( theTclInterp, theTokenizedExpression, &g_perturbed ) == TCL_ERROR) {
			opserr << "ERROR OpenSeesGradGEvaluator -- Tcl_ExprDouble returned error" << endln;
			opserr << theTclInterp->result << endln;
			return -1;
		      }

		      // Compute gradient
		      double onedgdu = (g_perturbed-g)/(originalValue*perturbationFactor);
			
		      //opserr << "OSGGE " << g_perturbed << ' ' << g << endln;




		      // Make assignment back to its original value
		      sprintf(tclAssignment,"set rec_element_%d_section_%d_force_%d %35.20f", eleNumber, sectionNumber, rowNumber, originalValue);
		      if (Tcl_Eval( theTclInterp, tclAssignment) == TCL_ERROR) {
			opserr << "ERROR OpenSeesGradGEvaluator -- Tcl_Eval returned error" << endln;
			opserr << theTclInterp->result << endln;
			return -1;
		      }
		      
		      // Obtain DDM gradient vector
		      RandomVariableIter &rvIter = 
			theReliabilityDomain->getRandomVariables();
		      RandomVariable *theRV;
		      //for (int i=1; i<=nrv; i++) {
		      while ((theRV = rvIter()) != 0) {
			int rvTag = theRV->getTag();
			//int i = theRV->getIndex();
			int i = theReliabilityDomain->getRandomVariableIndex(rvTag);
			sprintf(tclAssignment , "set sens [sensSectionForce %d %d %d %d]",eleNumber,sectionNumber,rowNumber,i);
			if (Tcl_Eval( theTclInterp, tclAssignment) == TCL_ERROR) {
			  opserr << "ERROR OpenSeesGradGEvaluator -- Tcl_Eval returned error" << endln;
			  opserr << theTclInterp->result << endln;
			  return -1;
			}
			sprintf(tclAssignment , "$sens ");
			if (Tcl_ExprDouble( theTclInterp, tclAssignment, &onedudx ) == TCL_ERROR) {
			  opserr << "ERROR OpenSeesGradGEvaluator -- Tcl_ExprDouble returned error" << endln;
			  opserr << theTclInterp->result << endln;
			    return -1;
			}
			dudx(i) = onedudx;
		      }
		      //opserr <<"OSGGE rec_element: " << onedgdu << ' ' << nrv << ' ' << dudx ;

		      // Add gradient contribution
		      (*grad_g) += onedgdu*dudx;
		    }
		    //opserr << "OSGGE: " << (*grad_g) << endln;
		  }
		  




		}

		tokenPtr = strtok( NULL, separators);  // read next token and go up and check the while condition again
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
OpenSeesGradGEvaluator::computeAllGradG(const Vector &gFunValues,
					const Vector &passed_x)
{
	
	numberOfEvalIncSens++;///// added by K Fujimura /////

  int sizeX = passed_x.Size();
	// Allocate result matrix
	Vector gradG(sizeX);
	if (grad_g_matrix == 0) {
		grad_g_matrix = new Matrix(sizeX, gFunValues.Size());
	}
	else {
		grad_g_matrix->Zero();
	}

	LimitStateFunctionIter lsfIter = theReliabilityDomain->getLimitStateFunctions();
	LimitStateFunction *theLSF;
	// Loop over performance functions
	//for (int j=1; j<=gFunValues.Size(); j++) {
	while ((theLSF = lsfIter()) != 0) {
	  int lsfTag = theLSF->getTag();
	  int j = theReliabilityDomain->getLimitStateFunctionIndex(lsfTag);

	  // Set tag of active limit-state function
	  theReliabilityDomain->setTagOfActiveLimitStateFunction(lsfTag);

	  this->computeGradG(gFunValues(j),passed_x);
	  const Vector &gradG = *grad_g;
	  
	  for (int i = 0; i < sizeX; i++) {
	    //opserr << "OSFFE: (i,j) = " << i << ' ' << j << endln;
	    (*grad_g_matrix)(i,j) = gradG(i);
	  }
	}

	return 0;
}


Matrix 
OpenSeesGradGEvaluator::getDgDdispl()
{
	return (*DgDdispl);
}

