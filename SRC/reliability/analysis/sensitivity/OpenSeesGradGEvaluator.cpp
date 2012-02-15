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
                                                                        
// $Revision: 1.21 $
// $Date: 2010-09-13 21:38:31 $
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
#include <Node.h>
#include <Element.h>
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
					       GFunEvaluator *passedGFunEvaluator,
						   ReliabilityDomain *passedReliabilityDomain,
						   Domain *passedOpenSeesDomain,
					       SensitivityAlgorithm *theAlgo,
					       bool PdoGradientCheck)
:GradGEvaluator(passedReliabilityDomain, passedGFunEvaluator, passedTclInterp), 
	theOpenSeesDomain(passedOpenSeesDomain)
{
	
	doGradientCheck = PdoGradientCheck;
	theSensAlgo = theAlgo;

	int nrv = passedReliabilityDomain->getNumberOfRandomVariables();
	grad_g = new Vector(nrv);
	grad_g_matrix = 0;

	DgDdispl = 0;

	/////S added by K Fujimura /////
	bool fdf = true;
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


double
OpenSeesGradGEvaluator::nodeGradient(int nodeNumber, int direction, int indx, char* dispOrWhat)
{
	
	// now obtain the information directly from domain
	Node *theNode = theOpenSeesDomain->getNode(nodeNumber);
	if (theNode == 0) {
		opserr << "GFunEvaluator::nodeTclVariable -- node with tag " << nodeNumber
			<< " not found in OpenSees Domain" << endln;
		return 0;
	}

	double gFunValue = 0.0;
	 
	if ( strncmp(dispOrWhat, "disp", 4) == 0) {
		gFunValue = theNode->getDispSensitivity(direction,indx);
	}
	else if ( strncmp(dispOrWhat, "vel", 3) == 0) {
		gFunValue = theNode->getVelSensitivity(direction,indx);
	}
	else if ( strncmp(dispOrWhat, "accel", 5) == 0) {
		gFunValue = theNode->getAccSensitivity(direction,indx);
	}
	else {
		opserr << "OpenSeesGradGEvaluator::nodeGradient in syntax (" << dispOrWhat << ") not supported" << endln;
	}

	return gFunValue;
}

double
OpenSeesGradGEvaluator::elementGradient(int eleNumber, int indx, char* inString)
{
	
	// now obtain the information directly from domain 
	Element *theElement = theOpenSeesDomain->getElement(eleNumber);
	if (theElement == 0) {
		opserr << "OpenSeesGradGEvaluator::elementGradient -- element with tag " << eleNumber
			<< " not found in OpenSees Domain" << endln;
		return 0;
	}
	
	double gFunValue = 0.0;
	int rowNumber; // index into vector containing element response
	char restString[100];
	strcpy(restString, inString);
	
	if ( strncmp(restString, "section",7) == 0) {
		int sectionNumber;
		sscanf(restString,"section_%i_%s", &sectionNumber, restString);
		
		if ( strncmp(restString, "force",5) == 0) {
			sscanf(restString,"force_%i", &rowNumber);
			//[sensSectionForce %d %d %d %d]",eleNumber,sectionNumber,rowNumber,i
			//gFunValue = theElement->getBlahBlahSectionForceSensitivity(sectionNumber,indx,rowNumber);
		}
		else {
			opserr << "OpenSeesGradGEvaluator::elementGradient in syntax (" << inString << ") not supported" << endln;
		}
	}
	else {
		// should work for any user input response quantity type (As long as element supports it)
		// No DDMs exist for these yet
		opserr << "OpenSeesGradGEvaluator::elementGradient in syntax (" << inString << ") not supported" << endln;
	}
	
	return gFunValue;
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
  // Perb Factor needs to be passed in by the user as well, should be implemented
  double perturbationFactor = 0.01; // (is multiplied by stdv and added to others...)
  int nrv = theReliabilityDomain->getNumberOfRandomVariables();
  RandomVariable *theRandomVariable;
  double g_perturbed;
  Vector dudx(nrv);
  char varName[20] = "";
  char arrName[40] = "";
  
  // Compute gradients if this is a path-INdependent analysis
  // (This command only has effect if it IS path-independent.)
  //if (theSensAlgo != 0 && !(theSensAlgo->shouldComputeAtEachStep()) ) {
  if (theSensAlgo != 0)
    theSensAlgo->computeSensitivities();
  //}
  
  // Initialize gradient vector
  grad_g->Zero();
  
  // "Download" limit-state function from reliability domain
  int lsf = theReliabilityDomain->getTagOfActiveLimitStateFunction();
  LimitStateFunction *theLimitStateFunction =
    theReliabilityDomain->getLimitStateFunctionPtr(lsf);
  
  // Loop through random variables
  for (int i = 0; i < nrv; i++) {
    double result = 0;

    RandomVariable *theRV =
      theReliabilityDomain->getRandomVariablePtrFromIndex(i);
    int tag = theRV->getTag();
    const char *gradExpression =
      theLimitStateFunction->getGradientExpression(tag);

    // There is an analytic expression for dg/dx
    if (gradExpression != 0) {
      if (Tcl_ExprDouble( theTclInterp, gradExpression, &result) == TCL_ERROR) {
	opserr << "ERROR OpenSeesGradGEvaluator -- error in Tcl_ExprDouble for the analytic gradient command" << gradExpression << endln;
	opserr << " caused error:" << endln << theTclInterp->result <<endln;
	return -1;
      }
      (*grad_g)(i) += result;
      //opserr << tag << ' ' << (*grad_g)(i) << endln;
    }

    // If no analytic expression, then determine dg/du by FD and du/dx by DDM
    else {
      char *theExpression = theLimitStateFunction->getExpression();
      opserr << theExpression << endln;
      Tcl_Obj *passedList = theLimitStateFunction->getParameters();
	
      // Cycle over all the LSF parameters and compute gradients
      Tcl_Obj *paramList = Tcl_DuplicateObj(passedList);
      int llength;
      Tcl_Obj *objPtr;
      if (Tcl_ListObjLength(theTclInterp,paramList,&llength) != TCL_OK) {
	opserr << "ERROR OpenSeesGradGEvaluator getting list length. " << endln;
	opserr << theTclInterp->result << endln;
      }
      
      for (int jk = 0; jk < llength; jk++) {
	Tcl_ListObjIndex(theTclInterp,paramList,jk,&objPtr);
	char *listStr = Tcl_GetStringFromObj(objPtr,NULL);
	opserr << jk << ' ' << listStr << endln;

	// If a RV is detected
	if ( strncmp(listStr, "x_", 2) == 0 || strncmp(listStr, "xrv", 3) == 0 ) {
	  
	  
	  // Get random variable tag
	  int rvTag = 0;
	  int oldNew = 0;
	  Tcl_Obj *outp;
	  int args = sscanf(listStr,"xrv(%i)",&rvTag);
	  if (args != 1) {
	    args = sscanf(listStr,"x_%i",&rvTag);
	    oldNew = 1;
	  }

	  // Perturb its value according to its standard deviation
	  theRandomVariable = theReliabilityDomain->getRandomVariablePtr(rvTag);
	  int rvIndex = theReliabilityDomain->getRandomVariableIndex(rvTag);
	  double stdv = theRandomVariable->getStdv();
	  double dx = perturbationFactor*stdv;
	  double x0 = passed_x(rvIndex);
	  double xval = x0 + dx;
	  
	  char theIndex[20];
	  sprintf(theIndex,"%d",rvTag);
	  if (oldNew == 0)
	    outp = Tcl_SetVar2Ex(theTclInterp,"xrv",theIndex,Tcl_NewDoubleObj(xval),TCL_LEAVE_ERR_MSG);
	  else
	    outp = Tcl_SetVar2Ex(theTclInterp,listStr,NULL,Tcl_NewDoubleObj(xval),TCL_LEAVE_ERR_MSG);
	  
	  if (outp == NULL) {
	    opserr << "ERROR OpenSeesGradGEvaluator -- error setting RV with tag " << rvTag << endln;
	    opserr << theTclInterp->result << endln;
	    return -1;
	  }
	  
	  //passed_x(rvIndex) = xval;
	  Vector xNew = passed_x;
	  xNew(rvIndex) = xval;
	  // Evaluate limit-state function again
	  if (theGFunEvaluator->evaluateG(xNew) < 0) {
	    opserr << "ERROR OpenSeesGradGEvaluator -- error evaluating LSF" << endln;
	    opserr << theTclInterp->result << endln;
	    return -1;
	  }
	  g_perturbed = theGFunEvaluator->getG();
	  
	  opserr << g << ' ' <<  g_perturbed << endln;
	  // Make assignment back to its original value
	  xval = x0;
	  if (oldNew == 0)
	    outp = Tcl_SetVar2Ex(theTclInterp,"xrv",theIndex,Tcl_NewDoubleObj(xval),TCL_LEAVE_ERR_MSG);
	  else
	    outp = Tcl_SetVar2Ex(theTclInterp,listStr,NULL,Tcl_NewDoubleObj(xval),TCL_LEAVE_ERR_MSG);
	  
	  if (outp == NULL) {
	    opserr << "ERROR OpenSeesGradGEvaluator -- error setting RV with tag " << rvTag << endln;
	    opserr << theTclInterp->result << endln;
	    return -1;
	  }
	  
	  // Add gradient contribution
	  (*grad_g)(rvIndex) += (g_perturbed-g)/dx;
	}
	
	/*
	// If a node quantity is detected
	else if ( strncmp(listStr, "u", 1) == 0 || 
		  strncmp(listStr, "node", 4) == 0 ||
		  strncmp(listStr, "rec_node", 8) == 0 ) {
	  
	  // Get node number, dof number, and type of response
	  // take advantage of GFunEvaluator parser, do NOT recode everything here
	  int nodeNumber, direction;
	  char dispOrWhat[10] = "";
	  if ( strncmp(listStr, "u", 1) == 0 )
	    theGFunEvaluator->uParse(listStr, &nodeNumber, &direction, dispOrWhat, varName, arrName);
	  else if ( strncmp(listStr, "node", 4) == 0 || strncmp(listStr, "rec_node", 8) == 0 )
	    theGFunEvaluator->nodeParse(listStr, &nodeNumber, &direction, dispOrWhat, varName, arrName);
	  
	  // retain original value
	  double originalValue = 0.0;
	  Tcl_Obj *tempDouble = Tcl_GetVar2Ex(theTclInterp, varName, arrName, TCL_LEAVE_ERR_MSG);
	  if (tempDouble == NULL) {
	    opserr << "ERROR OpenSeesGradGEvaluator -- GetVar error" << endln;
	    opserr << theTclInterp->result << endln;
	    return -1;
	  }
	  Tcl_GetDoubleFromObj(theTclInterp, tempDouble, &originalValue);
	  
	  // Set perturbed value in the Tcl workspace
	  // ---------------------- Quan and michele 
	  double newValue;
	  if (fabs(originalValue) < 1.0e-14 )
	    newValue = 0.0001;
	  else
	    newValue = originalValue*(1.0+perturbationFactor);
	  
	  // set Tcl value
	  if (Tcl_SetVar2Ex(theTclInterp,varName,arrName,Tcl_NewDoubleObj(newValue),TCL_LEAVE_ERR_MSG) == NULL) {
	    opserr << "ERROR OpenSeesGradGEvaluator -- SetVar error" << endln;
	    opserr << theTclInterp->result << endln;
	    return -1;
	  }
	  
	  // Evaluate the limit-state function again
	  if (theGFunEvaluator->evaluateGnoRecompute(theExpression) < 0) {
	    opserr << "ERROR OpenSeesGradGEvaluator -- error evaluating LSF" << endln;
	    opserr << theTclInterp->result << endln;
	    return -1;
	  }
	  g_perturbed = theGFunEvaluator->getG();
	  
	  // Compute gradient
	  double onedgdu = (g_perturbed-g)/(newValue-originalValue);
	  //opserr << "node dgdu = " << onedgdu << endln;
	  
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
	    for (int i=0; i<oldSize; i++) {
	      (*DgDdispl)(i,0) = tempMatrix(i,0);
	      (*DgDdispl)(i,1) = tempMatrix(i,1);
	      (*DgDdispl)(i,2) = tempMatrix(i,2);
	    }
	    (*DgDdispl)(oldSize,0) = (double)nodeNumber;
	    (*DgDdispl)(oldSize,1) = (double)direction;
	    (*DgDdispl)(oldSize,2) = onedgdu;
	  }
	  
	  // Make assignment back to its original value
	  theGFunEvaluator->nodeTclVariable(nodeNumber, direction, dispOrWhat, varName, arrName);
	  
	  // Obtain DDM gradient vector
	  // KRM -- note I changed this back to for-loop  because this rvIter is nested 
	  // above the rvIter going on inside GFunEvaluator
	  //RandomVariableIter &rvIter = theReliabilityDomain->getRandomVariables();
	  //RandomVariable *theRV;
	  for (int ig = 0; ig < nrv; ig++) {
	    //while ((theRV = rvIter()) != 0) {
	    //int rvTag = theRV->getTag();
	    //int i = theRV->getIndex();
	    //int idx = theReliabilityDomain->getRandomVariableIndex(rvTag);
	    dudx(ig) = nodeGradient(nodeNumber, direction, ig, dispOrWhat);
	  }
	  
	  // Add gradient contribution
	  (*grad_g) += onedgdu*dudx;
	  
	}     
	*/

	/*
	// If an element quantity is detected
	else if ( strncmp(listStr, "element", 7) == 0 || strncmp(listStr, "rec_element", 11) == 0 ) {
	  
	  //opserr << "DDM with element argument " << listStr << endln;
	  // Get element number and arguments
	  // take advantage of GFunEvaluator parser, do NOT recode everything here
	  int eleNumber = 0;
	  char restString[100] = "";
	  theGFunEvaluator->elementParse(listStr, &eleNumber, varName, restString);
	  
	  // retain original value
	  char arrName[100];
	  sprintf(arrName,"%d,%s",eleNumber,restString);
	  double originalValue = 0.0;
	  Tcl_Obj *tempDouble = Tcl_GetVar2Ex(theTclInterp, varName, arrName, TCL_LEAVE_ERR_MSG);
	  if (tempDouble == NULL) {
	    opserr << "ERROR OpenSeesGradGEvaluator -- GetVar error" << endln;
	    opserr << theTclInterp->result << endln;
	    return -1;
	  }
	  Tcl_GetDoubleFromObj(theTclInterp, tempDouble, &originalValue);
	  
	  // Set perturbed value in the Tcl workspace
	  double newValue = originalValue*(1.0+perturbationFactor);
	  
	  // set Tcl value
	  if (Tcl_SetVar2Ex(theTclInterp,varName,arrName,Tcl_NewDoubleObj(newValue),TCL_LEAVE_ERR_MSG) == NULL) {
	    opserr << "ERROR OpenSeesGradGEvaluator -- SetVar error" << endln;
	    opserr << theTclInterp->result << endln;
	    return -1;
	  }
	  
	  // Evaluate the limit-state function again
	  if (theGFunEvaluator->evaluateGnoRecompute(theExpression) < 0) {
	    opserr << "ERROR OpenSeesGradGEvaluator -- error evaluating LSF" << endln;
	    opserr << theTclInterp->result << endln;
	    return -1;
	  }
	  g_perturbed = theGFunEvaluator->getG();
	  
	  // Compute gradient
	  double onedgdu = (g_perturbed-g)/(newValue-originalValue);
	  
	  // Make assignment back to its original value
	  theGFunEvaluator->elementTclVariable(eleNumber, varName, restString);
	  
	  // Obtain DDM gradient vector
	  // KRM -- note I changed this back to for-loop  because this rvIter is nested 
	  // above the rvIter going on inside GFunEvaluator
	  //RandomVariableIter &rvIter = theReliabilityDomain->getRandomVariables();
	  //RandomVariable *theRV;
	  for (int ig = 0; ig < nrv; ig++) {
	    //while ((theRV = rvIter()) != 0) {
	    //int rvTag = theRV->getTag();
	    //int i = theRV->getIndex();
	    //int idx = theReliabilityDomain->getRandomVariableIndex(rvTag);
	    dudx(ig) = elementGradient(eleNumber, ig, restString);
	  }
	  //opserr << dudx;
	  
	  // Add gradient contribution
	  (*grad_g) += onedgdu*dudx;
	  
	}
*/
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
	// should never have an infinite loop
	while(true) {
	}
      }

      // reclaim Tcl object space
      Tcl_DecrRefCount(paramList);

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

	// not sure if this causes a conflict with lsfIter resetting somewhere else....TBD
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

